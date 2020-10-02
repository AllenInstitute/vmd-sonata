# Network layer4_450cells

This is a small network sampled from the 45,000 mouse V1 Layer 4 model: https://portal.brain-map.org/explore/models/l4-mv1. Not
including the feedforward inputs.

## Directory Structure:
 * _circuit_config.json_ - configuration file for initializing links including paths to components and network directories.

 * __components/__ - external files for initializing the network. Most importantly is the __components/morphologies/__ 
 which contain [swc files](http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html) for describing
 the morphology before rotation and translation.
 
 * __network/__ - SONATA format files for the network format.
 
 * _network_cells.csv_ and _network_synapses.csv_ - conversion of network's morphologically relevant into space-separated 
    csv file.
 
 
 ## network_cells.csv
 
 Space-separated csv file with the necessacary properties to display the network. Each row represents a different cell with
 a unique soma (cell body) location in cartesian space and rotation, although subsets of cells will share the same morphology.
 
 There are two different "cell type" models. Biophysically detail cells are a standard model of the cell with a full morphology
 including soma, axon and dendrities as described in the morphology file. Point_process cells are models of cells that do
 not contain an axon or an dendrities - and can be effectivily be thought of as single points or even just separated somas 
 (and thus they don't really have an euler rotation)
 
 #### columns
 
 |       |     |
 |-------|-----|
 | node_id| unique id given to each cell |
 | model_type| Either "biophysical" indicating a cell has an associated swc with full morphology, or "point_process" cell lacking axon and dendritic branches |
 | x, y, z | Cartesian coordinates at the center of a cell's soma |
 | rotation_angle_[x,y,z]_axis| rotation to be applied to biophysical cell. Assumes rotation is to be done first by z-axis, then a y-axis rotation, followed by x-axis rotation|
 | morphology| For biophysical cells a link to the (relative to this directory) location of the swc files containing cell morphologyies |
 


## network_synapses.csv

Contains location of where two cells synapse, usually represented by a "button" on part of the dendritic branch. This can 
be a nice feature to have but not as important a visualization compared to being able to show the cells and their relative
locations and rotations alone. 

In theory a synapse has a location on the pre-synaptic (source) neuron and post-synaptic (target) neuron. But in practice
most models tend to only be concerned with the location at the target-neuron and the source comes from part of the axon.

Doesn't apply for connections with post-synaptic "point_process" cells since their is only one place to put the synapse.

 #### columns

|     |      |
|-----|------|
|target_node_id| node_id of the post-synaptic node |
|source_node_id| node_id of the pre-synaptic node |
|section_id| identifier of cell section used by NEURON |
|section_x| Also used by NEURON, a value 0 - 1 to indicate where along a given section the synapse is placed|
|afferent_[x,y,z]| The location of the synpase in the cartesian space|


