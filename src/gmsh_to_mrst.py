'''
Run this program as

python gmsh_to_mrst.py input.msh output.mat

'''

import gmsh
import sys
import numpy
from scipy.io import savemat

gmsh.initialize()
model_name = sys.argv[1]

# Utilities
def reshape(x, cols):
    return x.reshape(len(x)//cols, cols)
def transpose(x):
    return x.reshape(len(x), 1)

# Load model
gmsh.open(model_name);
gdim = gmsh.model.getDimension()
print('Model ' + gmsh.model.getCurrent() + ' (' + str(gdim) + 'D)')
mrst_type = numpy.asarray(["Primordial Soup", "gmsh"], dtype=object)

# Nodes
node_tags, coords, _ = gmsh.model.mesh.getNodes()
num_nodes = len(numpy.unique(node_tags))
# Coords seem to be always 3, regardless of gdim
coords.resize(num_nodes, 3)
coords = coords[:, 0:gdim]
print(f"{num_nodes=}")

# Save to MRST data structure
mrst_nodes = {"num": float(num_nodes), 
              "coords": coords}

# Setup edges for all dimtags
dimtags_gdim = gmsh.model.getEntities(gdim)
gmsh.model.mesh.createEdges(dimtags_gdim)

# Setup all edges with default label
edge_labels = {}
default_edge_label = 0
num_nodes_per_edge = 2
for dimtag in dimtags_gdim:
    element_types, _, _ = gmsh.model.mesh.getElements(dimtag[0], dimtag[1])

    for element_type in element_types:
        edge_node_tags = gmsh.model.mesh.getElementEdgeNodes(element_type, dimtag[1])
        edge_nodes = numpy.sort(reshape(edge_node_tags, num_nodes_per_edge), axis=1)
        for en in edge_nodes:
            edge_labels[tuple(en)] = default_edge_label

# Setup edges with labels
dimtags_edges = gmsh.model.getEntities(gdim-1)
for dimtag in dimtags_edges:
    assert(dimtag[0] == 1)
    element_type, _, element_node_tags = gmsh.model.mesh.getElements(dimtag[0], dimtag[1])
    assert(numpy.all(element_type == 1)) # lines are of type 1
    for k in range(len(element_type)):
        edge_nodes = numpy.sort(reshape(element_node_tags[k], num_nodes_per_edge), axis=1)
        for en in edge_nodes:
            edge_labels[tuple(en)] = dimtag[1]

# For each element (only triangle or quad in 2D), collect the edges
# (as given by the nodes) and the elements (as given by the edges)
all_edges = []
all_edge_tags = []
all_elements_by_edges = []
all_element_labels = []
num_elements = 0

for dimtag in dimtags_gdim:
    # Let dimtag[1] denote the label
    label = dimtag[1]
    element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements(dimtag[0], dimtag[1])
    
    for eti, element_type in enumerate(element_types):
        prop = gmsh.model.mesh.getElementProperties(element_type)
        name = prop[0]
        print(f"{name=}", f"{label=}")
        num_nodes_per_element = prop[3]
        num_edges_per_element = num_nodes_per_element

        edge_node_tags = gmsh.model.mesh.getElementEdgeNodes(element_type, dimtag[1])
        edge_tags, _ = gmsh.model.mesh.getEdges(edge_node_tags)
        edges = numpy.sort(reshape(edge_node_tags, num_nodes_per_edge), axis=1)

        # Make these unique
        _, ii = numpy.unique(edge_tags, return_index=True)
        all_edge_tags.append(edge_tags[ii])
        all_edges.append(edges[ii,:])

        # Find representation of the elements as given by edges
        elements_by_edges = reshape(edge_tags, num_edges_per_element)
        all_elements_by_edges.append(elements_by_edges)
        print("num elements with this label", elements_by_edges.shape[0])
        num_elements += elements_by_edges.shape[0]
        element_labels = numpy.full(elements_by_edges.shape[0], label, dtype=int)
        all_element_labels.append(element_labels)
        
# Merge edges as given by the nodes
all_edges = numpy.concatenate(all_edges)
all_edge_tags = numpy.concatenate(all_edge_tags)
all_edge_tags, ii = numpy.unique(all_edge_tags, return_index=True)
all_edges = all_edges[ii,:]
num_edges = all_edges.shape[0]

# Store the labels
all_edge_labels = numpy.empty(num_edges, dtype=numpy.uint64)
for k, en in enumerate(all_edges):
    t = tuple(en)
    all_edge_labels[k] = edge_labels[t]

# Save to MRST data structure (only include tags if there are more than one unique tag)
nodePos = numpy.arange(1, num_nodes_per_edge*(num_edges+1), num_nodes_per_edge).astype(numpy.uint64)
mrst_faces = {"num": float(num_edges),
              "nodePos": transpose(nodePos), 
              "nodes": transpose(all_edges.flatten())}
if len(numpy.unique(all_edge_labels)) > 1:
    mrst_faces["tags"] = transpose(all_edge_labels)

# Save elements as given by the edges
diff = []
faces = []
for ee in all_elements_by_edges:
    diff.append(numpy.tile(ee.shape[1], (1, ee.shape[0]))[0])
    faces.append(ee.flatten())
diff = numpy.concatenate(diff)
facePos = numpy.append([1], numpy.cumsum(diff)+1).astype(numpy.uint64)
all_element_labels = numpy.concatenate(all_element_labels).astype(numpy.uint64)

# Save to MRST data structure (only include tags if there are more
# than one unique tag)
mrst_cells = {"num": float(num_elements),
              "facePos": transpose(facePos), 
              "faces": transpose(numpy.concatenate(faces))}
if len(numpy.unique(all_element_labels)) > 1:
    mrst_cells["tags"] = transpose(all_element_labels)

# Grid
G =  {"cells": mrst_cells,
      "faces": mrst_faces,
      "nodes": mrst_nodes,
      "type": mrst_type, 
      "griddim": float(gdim)}

if len(sys.argv) == 3:
    matfile = sys.argv[2]
else:
    matfile = "Gdata.mat"

savemat(matfile, G)
