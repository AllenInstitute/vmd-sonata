"""Simulates an example network of 450 cell receiving two kinds of exernal input as defined in the configuration file"""
import sys

from bmtk.simulator import bionet


def run(config_file):
    conf = bionet.Config.from_json(config_file, validate=True)
    conf.build_env()

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    sim.run()


if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        # Make sure to run only one at a time
        run('simulation_config.iclamp.json')  # Current clamp stimulation
