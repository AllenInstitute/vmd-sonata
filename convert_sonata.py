import os
from pathlib import PurePath
import pandas as pd
import numpy as np
from neuron import h
from six import string_types

from bmtk.utils import sonata
from bmtk.utils.sonata.config import SonataConfig
from bmtk.simulator.bionet import nrn
from bmtk.simulator.bionet.morphology import Morphology


class SWCReader(object):
    """A class for pulling out section id, section locations, coordinates from a SWC file. Useful when building a
    network that requires exact locations of pre- or post-synaptic locations. Requires NEURON.

    Attributes
    ==========
    swc_file - path to a SWC morphology file
    fix_axon - If set to true, the axon will be removed and replaced with a 30 um stub, as defined for all Allen
        Cell-Type models (default: True).
    random_seed - integer value to seed the random genator, used by choose_sections method.
    """

    def __init__(self, swc_file, random_seed=10, fix_axon=True):
        nrn.load_neuron_modules(None, None)
        self._swc_file = swc_file
        self._hobj = h.Biophys1(swc_file)
        if fix_axon:
            self._fix_axon()

        self._morphology = Morphology(self._hobj)
        self._morphology.set_seg_props()
        self._morphology.calc_seg_coords()
        self._prng = np.random.RandomState(random_seed)

        self._secs = []
        self._save_sections()

    def _save_sections(self):
        for sec in self._hobj.all:
            for _ in sec:
                self._secs.append(sec)

    def _fix_axon(self):
        """Removes and refixes axon"""
        axon_diams = [self._hobj.axon[0].diam, self._hobj.axon[0].diam]
        for sec in self._hobj.all:
            section_name = sec.name().split(".")[1][:4]
            if section_name == 'axon':
                axon_diams[1] = sec.diam

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

    def get_coord(self, sec_ids, sec_xs, soma_center=(0.0, 0.0, 0.0), rotation_matrix=None):
        """Takes in a list of section_ids and section_x values and returns a list of coordinates, assuming the soma
        is at the center of the system.

        :param sec_ids: [float]: list of N section_ids
        :param sec_xs: [float]: list of N cooresponding section_x's
        :param soma_center: location of soma in respect to the coordinate system. (default (0, 0, 0)).
        :param rotation_matrix: List of rotations (not yet implemented)
        :return: [(float, float, float)]: for seach sec_ids/sec_xs returna the x,y,z coordinates as a tuple
        """
        adjusted = self._morphology.get_soma_pos() - np.array(soma_center)
        absolute_coords = []
        for sec_id, sec_x in zip(sec_ids, sec_xs):
            sec = self._secs[sec_id]
            n_coords = int(h.n3d(sec=sec))
            coord_indx = int(sec_x*(n_coords - 1))
            swc_coords = np.array([h.x3d(coord_indx, sec=sec), h.y3d(coord_indx, sec=sec), h.x3d(coord_indx, sec=sec)])
            adjusted_coords = swc_coords - adjusted

            if rotation_matrix is not None:
                adjusted_coords = np.dot(rotation_matrix, adjusted_coords)

            absolute_coords.append(adjusted_coords)

        return absolute_coords

    def get_dist(self, sec_ids):
        """Returns arc-length distance from soma for a list of section_ids"""
        return [self._morphology.seg_prop['dist'][sec_id] for sec_id in sec_ids]


def rotation_matrix(axis, theta):
    """Return the rotation matrix associated with counterclockwise rotation about the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis * np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d

    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


class CellSections(object):
    def __init__(self, node_id, model_type, x, y, z, rot_x, rot_y, rot_z, swc_file, model_processing=None):
        self.node_id = node_id
        self.model_type = model_type
        self.x = x
        self.y = y
        self.z = z
        self.soma_center = (x, y, z)
        self.rot_x = rot_x
        self.rot_y = rot_y
        self.rot_z = rot_z
        self.swc_file = swc_file
        self.model_processing = model_processing
        self.is_biophysical = self.model_type == 'biophysical'

        rotx_mat = rotation_matrix([1, 0, 0], rot_x)
        roty_mat = rotation_matrix([0, 1, 0], rot_y)  # rotate segments around yaxis normal to pia
        rotz_mat = rotation_matrix([0, 0, 1], -rot_z)  # rotate segments around zaxis to get a proper orientation
        self.rotation_matrix = np.dot(rotx_mat, roty_mat.dot(rotz_mat))

        self.swc_reader = SWCReader(swc_file=swc_file, fix_axon=self.is_biophysical) if self.model_type == 'biophysical' else None

    def get_coords(self, sec_id, sec_x):
        return self.swc_reader.get_coord([sec_id], [sec_x], soma_center=self.soma_center,
                                         rotation_matrix=self.rotation_matrix)[0]

    @classmethod
    def load_row(cls, node_id, row, config_dir):
        morphology_path = os.path.join(config_dir, row['morphology'])
        return CellSections(
            node_id,
            row['model_type'],
            row['x'], row['y'], row['z'],
            row['rotation_angle_xaxis'], row['rotation_angle_yaxis'], row['rotation_angle_zaxis'],
            morphology_path,
            row['model_processing']
        )


def save_nodes_csv(circuit_config, population):
    config = SonataConfig.from_json(circuit_config)
    morphology_dir = config['components']['morphologies_dir']
    config_dir = config['manifest']['configdir']

    nodes_h5 = [n['nodes_file'] for n in config['networks']['nodes']]
    node_types_csv = [n['node_types_file'] for n in config['networks']['nodes']]
    l4_net = sonata.File(
        data_files=nodes_h5,
        data_type_files=node_types_csv
    )

    net_df = l4_net.nodes[population].to_dataframe()
    for rot_axis in ['rotation_angle_xaxis', 'rotation_angle_yaxis', 'rotation_angle_zaxis']:
        if rot_axis not in net_df.columns:
            net_df[rot_axis] = 0.0

        net_df[rot_axis] = net_df[rot_axis].fillna(0.0)

    net_df = net_df[['x', 'y', 'z', 'rotation_angle_xaxis', 'rotation_angle_yaxis', 'rotation_angle_zaxis',
                     'morphology', 'model_processing', 'model_name', 'model_type']]

    p = PurePath(morphology_dir)
    morp_rel_path = p.relative_to(config_dir)
    net_df['morphology'] = net_df.apply(
        lambda r: os.path.join(morp_rel_path, r['morphology']) if isinstance(r['morphology'], string_types) else None,
        axis=1
    )

    net_df.to_csv(os.path.join(config_dir, 'network_cells.csv'), sep=' ', na_rep="None")


def save_synapses_csv(circuit_config, population):
    config = SonataConfig.from_json(circuit_config)
    config_dir = config['manifest']['configdir']

    edges_h5 = [n['edges_file'] for n in config['networks']['edges']]
    edge_types_csv = [n['edge_types_file'] for n in config['networks']['edges']]

    l4_net = sonata.File(
        data_files=edges_h5,
        data_type_files=edge_types_csv
    )
    l4_edges = l4_net.edges['{0}_to_{0}'.format(population)]
    l4_nodes_df = pd.read_csv(os.path.join(config_dir, 'network_cells.csv'), sep=' ')
    cells = {}
    for nid, r in l4_nodes_df.iterrows():
        cells[nid] = CellSections.load_row(node_id=nid, row=r, config_dir=config_dir)

    src_node_ids = []
    trg_node_ids = []
    sec_ids = []
    sec_xs = []
    syn_x = []
    syn_y = []
    syn_z = []
    for e in l4_edges:
        src_node_ids.append(e.source_node_id)
        trg_node_ids.append(e.target_node_id)
        cell = cells[e.target_node_id]
        print(cell.node_id, cell.model_type, cell.is_biophysical)
        if cell.is_biophysical:
            cell_secs = cells[e.target_node_id]
            sec_id = e['sec_id']
            sec_x = e['sec_x']
            syn_coords = cell_secs.get_coords(sec_id, sec_x)

            sec_ids.append(sec_id)
            sec_xs.append(sec_x)
            syn_x.append(syn_coords[0])
            syn_y.append(syn_coords[1])
            syn_z.append(syn_coords[2])
        else:
            sec_ids.append(-1)
            sec_xs.append(-1)
            syn_x.append(np.nan)
            syn_y.append(np.nan)
            syn_z.append(np.nan)

    edges_df = pd.DataFrame({
        'target_node_id': trg_node_ids,
        'source_node_id': src_node_ids,
        'section_id': sec_ids,
        'section_x': sec_xs,
        'afferent_x': syn_x,
        'afferent_y': syn_y,
        'afferent_z': syn_z
    })
    edges_df.to_csv(os.path.join(config_dir, 'network_synapses.csv'), sep=' ', index=False, na_rep=np.nan)


if __name__ == '__main__':
    save_nodes_csv(circuit_config='layer4_450cells/circuit_config.json', population='internal')
    save_synapses_csv(circuit_config='layer4_450cells/circuit_config.json', population='internal')
