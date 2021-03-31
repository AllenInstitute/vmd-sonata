import os
import h5py
import pandas as pd
import numpy as np
from itertools import combinations, permutations
from neuron import h

from bmtk.simulator.bionet import nrn
from bmtk.utils import sonata
from bmtk.builder import NetworkBuilder
from bmtk.simulator.bionet.morphology import Morphology


pd.set_option('display.max_columns', None)


def sample_df():
    sonata_file = sonata.File(
        data_files=[
            '/data/work_files/V1_network_update/Biophysical_network/network_fixed_synapses/v1_nodes.h5',
            '/data/work_files/V1_network_update/Biophysical_network/network_fixed_synapses/v1_v1_edges.h5'
        ],
        data_type_files=[
            '/data/work_files/V1_network_update/Biophysical_network/network_fixed_synapses/v1_node_types.csv',
            '/data/work_files/V1_network_update/Biophysical_network/network_fixed_synapses/v1_v1_edge_types.csv'
        ]
    )

    nodes_df = sonata_file.nodes['v1'].to_dataframe(index_by_id=False)
    nodes_df = nodes_df[['node_id', 'location', 'ei', 'model_type']]

    h5 = h5py.File('/data/work_files/V1_network_update/Biophysical_network/network_fixed_synapses/v1_v1_edges.h5', 'r')
    edges_h5 = h5['/edges/v1_to_v1']
    edges_df = pd.DataFrame({
        'target_node_id': edges_h5['target_node_id'],
        'source_node_id': edges_h5['source_node_id']
    })
    # edges_df = sonata_file.edges['v1_to_v1'].get_group(0).to_dataframe()
    edges_df = edges_df.merge(nodes_df, how='left', left_on='target_node_id', right_on='node_id')
    edges_df = edges_df[edges_df['model_type'] == 'biophysical']
    edges_df = edges_df.drop(columns=['node_id', 'model_type'])
    edges_df = edges_df.rename(columns={'location': 'trg_location', 'ei': 'trg_ei', 'morphology': 'trg_morphology'})
    print(edges_df)

    edges_df = edges_df.merge(nodes_df, how='left', left_on='source_node_id', right_on='node_id')
    edges_df = edges_df[edges_df['model_type'] == 'biophysical']
    edges_df = edges_df.drop(columns=['node_id', 'model_type'])
    edges_df = edges_df.rename(columns={'location': 'src_location', 'ei': 'src_ei', 'morphology': 'src_morphology'})
    print(edges_df)

    e2e_edges_df = edges_df[(edges_df['src_ei'] == 'e') & (edges_df['trg_ei'] == 'e')]
    l4_l2_edges_df = e2e_edges_df[
        (e2e_edges_df['src_location'] == 'VisL4') & (e2e_edges_df['trg_location'] == 'VisL23')
    ]

    l2_l4_edges_df = e2e_edges_df[
        (e2e_edges_df['src_location'] == 'VisL23') & (e2e_edges_df['trg_location'] == 'VisL4')
    ]

    combined_edges_df = pd.concat([l2_l4_edges_df, l4_l2_edges_df])
    combined_edges_df.to_csv('metadata/combined_edges.csv', sep=' ', index=False)


def select_nodes():
    sonata_file = sonata.File(
        data_files=['/data/work_files/V1_network_update/Biophysical_network/network_fixed_synapses/v1_nodes.h5'],
        data_type_files=['/data/work_files/V1_network_update/Biophysical_network/network_fixed_synapses/v1_node_types.csv']
    )
    nodes_df = sonata_file.nodes['v1'].to_dataframe(index_by_id=False)[['node_id', 'pop_name']]

    edges_df = pd.read_csv('metadata/combined_edges.csv', sep=' ')

    l4_l2_edges_df = edges_df[edges_df['src_location'] == 'VisL4'][['source_node_id', 'target_node_id']]
    l4_l2_edges_df = l4_l2_edges_df.rename(columns={'source_node_id': 'l4_node_id', 'target_node_id': 'l2_node_id'})
    l4_l2_counts_df = l4_l2_edges_df.groupby(['l4_node_id', 'l2_node_id']).size().to_frame('->').reset_index()
    l4_l2_counts_df = l4_l2_counts_df[l4_l2_counts_df['->'] > 3]

    l2_l4_edges_df = edges_df[edges_df['src_location'] == 'VisL23']
    l2_l4_edges_df = l2_l4_edges_df.rename(columns={'source_node_id': 'l2_node_id', 'target_node_id': 'l4_node_id'})
    l2_l4_counts_df = l2_l4_edges_df.groupby(['l4_node_id', 'l2_node_id']).size().to_frame('<-').reset_index()

    sym_counts_df = l4_l2_counts_df.merge(l2_l4_counts_df, how='inner', on=['l4_node_id', 'l2_node_id'])
    sym_counts_df = sym_counts_df.merge(nodes_df, how='left', left_on='l4_node_id', right_on='node_id')
    sym_counts_df = sym_counts_df.rename(columns={'pop_name': 'src_pop_name'})

    sym_counts_df = sym_counts_df.merge(nodes_df, how='left', left_on='l2_node_id', right_on='node_id')
    sym_counts_df = sym_counts_df.rename(columns={'pop_name': 'trg_pop_name'})

    sym_counts_df = sym_counts_df.sort_values('->')
    sym_counts_df.to_csv('metadata/symmetric_edge_counts.csv', sep=' ', index=False)
    # print(l4_l2_counts_df.merge(l2_l4_counts_df, how='inner', on=['l4_node_id', 'l2_node_id']).sort_values('->'))
    print(sym_counts_df)


class SWCParser(object):
    def __init__(self, swc_path, fix_axon=True):
        self._swc_path = swc_path
        self._swc_df = pd.read_csv(swc_path, sep=' ', names=['id', 'type', 'x', 'y', 'z', 'r', 'pid'], comment='#')

        nrn.load_neuron_modules(None, None)
        self._hobj = h.Biophys1(swc_path)
        self._sections = []
        self._swc_ids = []
        section_list = list(self._hobj.all)
        for ln in range(len(self._swc_df)):
            sec_id = int(self._hobj.nl.point2sec[ln])
            swc_id = int(self._hobj.nl.pt2id(ln))
            sec = section_list[sec_id]
            self._sections.append(sec.name())
            self._swc_ids.append(swc_id)

        self._swc_df['section_name'] = self._sections

        # self._swc_table = pd.DataFrame({
        #     'section_name': self._sections,
        #     'swc_id': self._swc_ids
        # }).

        if fix_axon:
            self._fix_axon()
            axon_indices = [i for i, sec_name in enumerate(self._sections) if 'axon' in sec_name]
            self._sections = [self._sections[i] for i in range(len(self._sections)) if i not in axon_indices]
            self._swc_ids = [self._swc_ids[i] for i in range(len(self._swc_ids)) if i not in axon_indices]

    def get_swc_id(self, sec_id, sec_x):
        section_name = self._sections[sec_id]
        swc_rows = self._swc_df[self._swc_df['section_name'] == section_name]
        swc_ids = swc_rows['id'].values
        swc_place = np.max((0.0, sec_x*len(swc_ids) - 1.0))
        swc_indx = int(np.ceil(swc_place))
        swc_id = swc_ids[swc_indx]
        swc_dist = swc_indx - swc_place

        print(self._swc_path, section_name, swc_id)
        return swc_id, swc_dist

    def _fix_axon(self):
        """Removes and refixes axon"""
        axon_diams = [self._hobj.axon[0].diam, self._hobj.axon[0].diam]
        axon_indices = []
        for i, sec in enumerate(self._hobj.all):
            section_name = sec.name().split(".")[1][:4]
            if section_name == 'axon':
                axon_diams[1] = sec.diam
                axon_indices.append(i)

        for sec in self._hobj.axon:
            h.delete_section(sec=sec)

        h.execute('create axon[2]', self._hobj)
        for index, sec in enumerate(self._hobj.axon):
            sec.L = 30
            sec.diam = 1

            self._hobj.axonal.append(sec=sec)
            self._hobj.all.append(sec=sec)  # need to remove this comment

        self._hobj.axon[0].connect(self._hobj.soma[0], 1.0, 0)
        self._hobj.axon[1].connect(self._hobj.axon[0], 1.0, 0)

        h.define_shape()


def find_branches(swc_file):
    swc_df = pd.read_csv(swc_file, sep=' ', names=['id', 'type', 'x', 'y', 'z', 'r', 'pid'], comment='#')
    dends_df = swc_df[swc_df['type'] == 3]
    dend_ids = dends_df['id'].values
    branch_ends = np.sort(dends_df[(dends_df['id'] - dends_df['pid'] > 1) & (dends_df['pid'] > 1)]['pid'].values)
    # branch_ids = np.sort(branch_ids)
    branch_begs = np.sort(dends_df[dends_df['pid'].isin(branch_ends)]['id'].values)
    branch_begs = np.concatenate(([dend_ids.min()], branch_begs))
    print(np.diff(branch_begs))
    print(branch_begs)
    print(branch_ends)


def build_edges(source_node_id, target_node_id):
    sonata_file = sonata.File(
        data_files=[
            '/data/work_files/V1_network_update/Biophysical_network/network_fixed_synapses/v1_nodes.h5',
            '/data/work_files/V1_network_update/Biophysical_network/network_fixed_synapses/v1_v1_edges.h5'
        ],
        data_type_files=[
            '/data/work_files/V1_network_update/Biophysical_network/network_fixed_synapses/v1_node_types.csv',
            '/data/work_files/V1_network_update/Biophysical_network/network_fixed_synapses/v1_v1_edge_types.csv'
        ]
    )

    v1_edges = sonata_file.edges['v1_to_v1']
    v1_nodes = sonata_file.nodes['v1']

    l4_node = v1_nodes.get_node_id(source_node_id)
    l4_node_morphology_path = os.path.join('biophys_components/morphologies', l4_node['morphology'])
    l4_node_props = l4_node.node_type_properties
    l4_node_props.update({str(k): [v] for k, v in l4_node.group_props.items()})
    l4_node_props['rotation_angle_zaxis'] *= -1

    l2_node = v1_nodes.get_node_id(target_node_id)
    l2_node_morphology_path = os.path.join('biophys_components/morphologies', l2_node['morphology'])
    l2_node_props = l2_node.node_type_properties
    l2_node_props.update({str(k): [v] for k, v in l2_node.group_props.items()})
    l2_node_props['rotation_angle_zaxis'] *= -1

    swc_parsers = {
        source_node_id: SWCParser(l4_node_morphology_path),
        target_node_id: SWCParser(l2_node_morphology_path)
    }

    v1 = NetworkBuilder('v1')
    v1.add_nodes(**l4_node_props)
    v1.add_nodes(**l2_node_props)

    def query_str2dict(qry_str):
        qry_dict = {}
        for cond in qry_str.split('&'):
            k, v = cond.split('==')
            v = v.strip("\'")
            try:
                v = int(v)
            except ValueError:
                pass

            qry_dict[k] = v
        return qry_dict

    for src_id, trg_id in permutations([source_node_id, target_node_id]):
        sec_ids = []
        sec_xs = []
        syn_weights = []
        edge_props = {}
        trg_swc = swc_parsers[trg_id]
        trg_swc_ids = []
        trg_swc_dists = []

        for e in v1_edges.get_source(src_id):
            if e.target_node_id == trg_id:
                sec_id = e._group_props['sec_id']
                sec_x = e._group_props['sec_x']
                swc_id, swc_dist = trg_swc.get_swc_id(sec_id, sec_x)

                trg_swc_ids.append(swc_id)
                trg_swc_dists.append(swc_dist)
                edge_props[e._edge_type_props['edge_type_id']] = e._edge_type_props
                sec_ids.append(sec_id)
                sec_xs.append(sec_x)
                syn_weights.append(e._group_props['syn_weight'])

        for edge_type_id, edge_type_props in edge_props.items():
            trg_query = query_str2dict(edge_type_props['target_query'])
            del edge_type_props['target_query']
            src_query = query_str2dict(edge_type_props['source_query'])
            del edge_type_props['source_query']

            cm = v1.add_edges(
                source=src_query,
                target=trg_query,
                connection_rule=len(syn_weights),
                **edge_type_props
            )
            # cm.add_properties('syn_weight', rule=lambda *_: np.random.uniform(0.0, 100.0), dtypes=np.float)
            cm.add_properties('afferent_section_id', rule=sec_ids, dtypes=np.int)
            cm.add_properties('afferent_section_pos', rule=sec_xs, dtypes=np.float)
            cm.add_properties('syn_weight', rule=syn_weights, dtypes=np.float)
            cm.add_properties('afferent_swc_id', rule=trg_swc_ids, dtypes=np.int)
            cm.add_properties('afferent_swc_pos', rule=trg_swc_dists, dtypes=np.float)

    v1.build()
    v1.save(output_dir='network')

def sonata2csv():
    sonata_file = sonata.File(
        data_files=['network/v1_nodes.h5', 'network/v1_v1_edges.h5'],
        data_type_files=['network/v1_node_types.csv', 'network/v1_v1_edge_types.csv']
    )
    v1_nodes = sonata_file.nodes['v1']
    nodes_df = v1_nodes.to_dataframe(index_by_id=False)
    nodes_df.to_csv('network/csv/nodes.csv', sep=' ', index=False)

    v1_edges = sonata_file.edges['v1_to_v1']
    edges_df = v1_edges.get_group(0).to_dataframe()
    edges_df.to_csv('network/csv/edges.csv', sep=' ', index=False)




if __name__ == '__main__':
    # sample_df()
    # select_nodes()
    # l4_node_id = 122801
    # l2_node_id = 71546
    l4_node_id = 110435
    l2_node_id = 72492

    build_edges(source_node_id=l4_node_id, target_node_id=l2_node_id)
    sonata2csv()
