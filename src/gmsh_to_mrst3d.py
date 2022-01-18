"""
Run this program as

python gmsh_to_mrst.py input.msh output.mat

"""

import gmsh
import sys
import numpy
from scipy.io import savemat

gmsh.initialize()
# model_name = sys.argv[1]
model_name = "vwell25.msh"

# Utilities
def reshape(x, cols):
    return x.reshape(len(x) // cols, cols)


def transpose(x):
    return x.reshape(len(x), 1)


# Load model
gmsh.open(model_name)
gdim = gmsh.model.getDimension()
assert gdim == 3, "This script is currently only for 3D"
print("Model " + gmsh.model.getCurrent() + " (" + str(gdim) + "D)")
mrst_type = numpy.asarray(["Primordial Soup", "gmsh"], dtype=object)

# Nodes
node_tags, coords, _ = gmsh.model.mesh.getNodes()
num_nodes = len(numpy.unique(node_tags))
# Coords seem to be always 3, regardless of gdim
coords.resize(num_nodes, 3)
coords = coords[:, 0:gdim]
print(f"{num_nodes=}")

# Save to MRST data structure
mrst_nodes = {"num": float(num_nodes), "coords": coords}

# Setup faces for all dimtags
dimtags_gdim = gmsh.model.getEntities(gdim)
# gmsh.model.mesh.createEdges(dimtags_gdim)
gmsh.model.mesh.createFaces(dimtags_gdim)

# Setup all faces with default label
face_labels = {}
default_face_label = 0

# Face types are triangle or quad
face_types = [3, 4]

# Mmapping of the face type (see gmsh manual) to the number of nodes
# per face (3 for triangular faces, 4 for quadrangular faces). This is
# for use in getElementFaceNodes only.
num_nodes_per_face_type = {3: 3, 4: 4}

for dimtag in dimtags_gdim:
    element_types, _, _ = gmsh.model.mesh.getElements(dimtag[0], dimtag[1])

    for element_type in element_types:
        for face_type in face_types:
            face_node_tags = gmsh.model.mesh.getElementFaceNodes(
                element_type, face_type
            )
            if face_node_tags.size:
                face_nodes = numpy.sort(
                    reshape(face_node_tags, num_nodes_per_face_type[face_type]), axis=1
                )
                for fn in face_nodes:
                    face_labels[tuple(fn)] = default_face_label

# Setup faces with labels
dimtags_faces = gmsh.model.getEntities(gdim - 1)
for dimtag in dimtags_faces:
    assert dimtag[0] == 2
    element_type, _, element_node_tags = gmsh.model.mesh.getElements(
        dimtag[0], dimtag[1]
    )
    # Element type is typically 3 for quadrilateral, 2 for triangle
    assert numpy.all(element_type == 3 or element_type == 2)
    for k in range(len(element_type)):
        prop = gmsh.model.mesh.getElementProperties(element_type[k])
        num_nodes_per_element = prop[3]
        face_nodes = numpy.sort(
            reshape(element_node_tags[k], num_nodes_per_element), axis=1
        )
        for fn in face_nodes:
            face_labels[tuple(fn)] = dimtag[1]

# For each element, collect the faces (as given by the nodes) and the
# elements (as given by the faces)
all_faces = []
all_face_tags = []
all_elements_by_faces = []
all_element_labels = []
num_elements = 0

# Mapping from element type to the number of faces per element (eg
# hexahedron is type 5, has 6 faces)
num_faces_per_element = {5: 6}

for dimtag in dimtags_gdim:
    # Let dimtag[1] denote the label
    label = dimtag[1]
    element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements(
        dimtag[0], dimtag[1]
    )

    for eti, element_type in enumerate(element_types):
        prop = gmsh.model.mesh.getElementProperties(element_type)
        name = prop[0]
        print(f"{name=}", f"{label=}")
        num_nodes_per_element = prop[3]

        for face_type in face_types:
            face_node_tags = gmsh.model.mesh.getElementFaceNodes(
                element_type, face_type
            )
            if face_node_tags.size:
                face_tags, _ = gmsh.model.mesh.getFaces(face_type, face_node_tags)
                faces = numpy.sort(
                    reshape(face_node_tags, num_nodes_per_face_type[face_type]), axis=1
                )

                # Make these unique
                _, ii = numpy.unique(face_tags, return_index=True)
                all_face_tags.append(face_tags[ii])
                all_faces.append(faces[ii, :])

                # Find representation of the elements as given by the faces
                # FIXME we assume that all elements have only one face type
                elements_by_faces = reshape(
                    face_tags, num_faces_per_element[element_type]
                )
                all_elements_by_faces.append(elements_by_faces)
                print("num elements with this label", elements_by_faces.shape[0])
                num_elements += elements_by_faces.shape[0]
                element_labels = numpy.full(
                    elements_by_faces.shape[0], label, dtype=int
                )
                all_element_labels.append(element_labels)

# Merge faces as given by the nodes
all_faces = numpy.concatenate(all_faces)
all_face_tags = numpy.concatenate(all_face_tags)
all_face_tags, ii = numpy.unique(all_face_tags, return_index=True)
all_faces = all_faces[ii, :]
num_faces = all_faces.shape[0]

# Store the labels
all_face_labels = numpy.empty(num_faces, dtype=numpy.uint64)
for k, fn in enumerate(all_faces):
    t = tuple(fn)
    all_face_labels[k] = face_labels[t]

# Save to MRST data structure (only include tags if there are more than one unique tag)
num_nodes_per_face = 4 # FIXME 
nodePos = numpy.arange(
    1, num_nodes_per_face * (num_faces + 1), num_nodes_per_face
).astype(numpy.uint64)
mrst_faces = {
    "num": float(num_faces),
    "nodePos": transpose(nodePos).astype(float),
    "nodes": transpose(all_faces.flatten()).astype(float),
}
if len(numpy.unique(all_face_labels)) > 1:
    mrst_faces["tags"] = transpose(all_face_labels).astype(float)

# Save elements as given by the faces
# FIXME uniform elements
diff = []
faces = []
for ee in all_elements_by_faces:
    diff.append(numpy.tile(ee.shape[1], (1, ee.shape[0]))[0])
    faces.append(ee.flatten())
diff = numpy.concatenate(diff)
facePos = numpy.append([1], numpy.cumsum(diff) + 1).astype(numpy.uint64)
all_element_labels = numpy.concatenate(all_element_labels).astype(numpy.uint64)

# Save to MRST data structure (only include tags if there are more
# than one unique tag)
mrst_cells = {
    "num": float(num_elements),
    "facePos": transpose(facePos).astype(float),
    "faces": transpose(numpy.concatenate(faces)).astype(float),
}
if len(numpy.unique(all_element_labels)) > 1:
    mrst_cells["tags"] = transpose(all_element_labels).astype(float)

# Grid
G = {
    "cells": mrst_cells,
    "faces": mrst_faces,
    "nodes": mrst_nodes,
    "type": mrst_type,
    "griddim": float(gdim),
}

if len(sys.argv) == 3:
    matfile = sys.argv[2]
else:
    matfile = "Gdata.mat"

savemat(matfile, G)
