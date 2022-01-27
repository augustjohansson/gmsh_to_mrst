"""
Run this program as

python gmsh_to_mrst.py input.msh output.mat

"""

import sys
import numpy as np
import gmsh
from scipy.io import savemat

# Set default markers (we don't use the work tag since it's used in
# gmsh for something different)
default_edge_marker = 0
default_face_marker = 0
default_cell_marker = 0

gmsh.initialize()
model_name = sys.argv[1]


def reshape(array, cols):
    return array.reshape(len(array) // cols, cols)


def transpose(array):
    return array.reshape(len(array), 1)


def flatten(lst):
    return [elem for sublist in lst for elem in sublist]


# Load model
gmsh.open(model_name)
gdim = gmsh.model.getDimension()
print("Model " + gmsh.model.getCurrent() + " (" + str(gdim) + "D)")

# Nodes
node_tags, coords, _ = gmsh.model.mesh.getNodes()
num_nodes = len(np.unique(node_tags))

# Coords are always in 3 columns, regardless of gdim
coords.resize(num_nodes, 3)
coords = coords[:, 0:gdim]
print(f"{num_nodes=}")

# MRST grid nodes
mrst_nodes = {"num": float(num_nodes), "coords": coords}

# Separate functions for gdim=2 and gdim=3
if gdim == 2:
    # Setup edges for all dimtags
    dimtags_gdim = gmsh.model.getEntities(gdim)
    gmsh.model.mesh.createEdges(dimtags_gdim)

    # Setup all edges with default label
    edge_markers = {}
    num_nodes_per_edge = 2
    for dimtag in dimtags_gdim:
        element_types, _, _ = gmsh.model.mesh.getElements(dimtag[0], dimtag[1])

        for element_type in element_types:
            edge_node_tags = gmsh.model.mesh.getElementEdgeNodes(
                element_type, dimtag[1]
            )
            edge_nodes = np.sort(reshape(edge_node_tags, num_nodes_per_edge), axis=1)
            for en in edge_nodes:
                edge_markers[tuple(en)] = default_edge_marker

    # Setup edges with labels
    dimtags_edges = gmsh.model.getEntities(gdim - 1)
    for dimtag in dimtags_edges:
        assert dimtag[0] == 1
        element_type, _, element_node_tags = gmsh.model.mesh.getElements(
            dimtag[0], dimtag[1]
        )
        assert np.all(element_type == 1)  # lines are of type 1
        for k in range(len(element_type)):
            edge_nodes = np.sort(
                reshape(element_node_tags[k], num_nodes_per_edge), axis=1
            )
            for en in edge_nodes:
                edge_markers[tuple(en)] = dimtag[1]

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
        element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements(
            dimtag[0], dimtag[1]
        )

        for eti, element_type in enumerate(element_types):
            prop = gmsh.model.mesh.getElementProperties(element_type)
            name = prop[0]
            print(f"{name=}", f"{label=}")
            num_nodes_per_element = prop[3]
            num_edges_per_element = num_nodes_per_element

            edge_node_tags = gmsh.model.mesh.getElementEdgeNodes(
                element_type, dimtag[1]
            )
            edge_tags, _ = gmsh.model.mesh.getEdges(edge_node_tags)
            edges = np.sort(reshape(edge_node_tags, num_nodes_per_edge), axis=1)

            # Make these unique
            _, ii = np.unique(edge_tags, return_index=True)
            all_edge_tags.append(edge_tags[ii])
            all_edges.append(edges[ii, :])

            # Find representation of the elements as given by edges
            elements_by_edges = reshape(edge_tags, num_edges_per_element)
            all_elements_by_edges.append(elements_by_edges)
            print("num elements with this label", elements_by_edges.shape[0])
            num_elements += elements_by_edges.shape[0]
            element_labels = np.full(elements_by_edges.shape[0], label, dtype=int)
            all_element_labels.append(element_labels)

    # Merge edges as given by the nodes
    all_edges = np.concatenate(all_edges)
    all_edge_tags = np.concatenate(all_edge_tags)
    all_edge_tags, ii = np.unique(all_edge_tags, return_index=True)
    all_edges = all_edges[ii, :]
    num_edges = all_edges.shape[0]

    # Store the labels
    all_edge_markers = np.empty(num_edges, dtype=np.uint64)
    for k, en in enumerate(all_edges):
        t = tuple(en)
        all_edge_markers[k] = edge_markers[t]

    # Save to MRST data structure (only include tags if there are more than one unique tag)
    nodePos = np.arange(
        1, num_nodes_per_edge * (num_edges + 1), num_nodes_per_edge
    ).astype(np.uint64)
    mrst_faces = {
        "num": float(num_edges),
        "nodePos": transpose(nodePos).astype(float),
        "nodes": transpose(all_edges.flatten()).astype(float),
    }
    if len(np.unique(all_edge_markers)) > 1:
        mrst_faces["tags"] = transpose(all_edge_markers).astype(float)

    # Save elements as given by the edges
    diff = []
    faces = []
    for ee in all_elements_by_edges:
        diff.append(np.tile(ee.shape[1], (1, ee.shape[0]))[0])
        faces.append(ee.flatten())
    diff = np.concatenate(diff)
    facePos = np.append([1], np.cumsum(diff) + 1).astype(np.uint64)
    all_element_labels = np.concatenate(all_element_labels).astype(np.uint64)

    # Save to MRST data structure (only include tags if there are more
    # than one unique tag)
    mrst_cells = {
        "num": float(num_elements),
        "facePos": transpose(facePos).astype(float),
        "faces": transpose(np.concatenate(faces)).astype(float),
    }
    if len(np.unique(all_element_labels)) > 1:
        mrst_cells["tags"] = transpose(all_element_labels).astype(float)

elif gdim == 3:

    # Get total number of faces (tris or quads)
    gmsh.model.mesh.createFaces()
    face_types = [3, 4]
    face_types_names = ["Triangle", "Quadrilateral"]
    num_faces = 0
    for face_type, name in zip(face_types, face_types_names):
        ff = gmsh.model.mesh.getAllFaces(face_type)[0]
        num_faces += len(ff)
        print("\tnumber of face type '", name, "':", len(ff))
    print(f"{num_faces=}")

    # Get total number of cells
    num_cells = 0
    element_types, all_element_tags, _ = gmsh.model.mesh.getElements(gdim)
    for element_type, element_tags in zip(element_types, all_element_tags):
        num_cells += len(element_tags)
        prop = gmsh.model.mesh.getElementProperties(element_type)
        print("\tnumber of element type '", prop[0], "':", len(element_tags))
    print(f"{num_cells=}")

    # Setup entities of codimension 1 with the dimtag as marker. Note that
    # not all the faces are labeled, only those surfaces manually created
    # in the modeling (for example, if a volume is created by extruding a
    # surface, only the faces of this surface are labeled).

    # Identify face markers by the nodes representing the face
    face_markers_by_nodes = {}

    for dimtag in gmsh.model.getEntities(gdim - 1):
        element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements(
            dimtag[0], dimtag[1]
        )

        # Element type for the faces is either 2 or 3 (tris or quads)
        assert np.all(np.logical_or(element_types == 2, element_types == 3))

        # The face labels are identified by the set of nodes
        for element_type, element_node_tag in zip(element_types, element_node_tags):
            prop = gmsh.model.mesh.getElementProperties(element_type)
            num_nodes_per_element = prop[3]
            for fn in reshape(element_node_tag, num_nodes_per_element):
                face_markers_by_nodes[frozenset(fn)] = dimtag[1]

    # Create the cells to faces map and the faces to nodes map
    print("Create cell to face and face to node maps")
    cells_2_faces = {}
    faces_2_nodes = {}
    nodes_2_faces = {}

    # Save cell markers
    cell_markers_by_tag = {}

    # Mapping from element type to the number of faces per element (eg
    # hexahedron is type 5, has 6 faces)
    # num_faces_per_element = {5: 6, 4: 4, 7: 5, 6: 5}
    num_faces_per_element = {5: 6, 4: 4}

    # For the prism and pyramid we need how many faces of each kind (tri
    # or quads)
    num_faces_prism = {3: 2, 4: 3}
    num_faces_pyramid = {3: 4, 4: 1}

    # Mapping of the face type (not the element type; cf gmsh manual) to
    # the number of nodes per face
    num_nodes_per_face_type = {3: 3, 4: 4}

    for dimtag in gmsh.model.getEntities(gdim):
        dim = dimtag[0]
        label = dimtag[1]

        # Get all element types in this domain
        element_types, _, _ = gmsh.model.mesh.getElements(dim, label)

        for element_type in element_types:

            # Face tags are retrieved in same order as elements
            element_tags, element_node_tags = gmsh.model.mesh.getElementsByType(
                element_type
            )

            # Set the label for all these elements
            for cell in element_tags:
                cell_markers_by_tag[cell] = label

            # Loop over all face types (tris or quads)
            for face_type in face_types:
                face_node_tags = gmsh.model.mesh.getElementFaceNodes(
                    element_type, face_type
                )

                # Some element types may not have any faces of this type
                if face_node_tags.size:
                    face_tags, face_orientations = gmsh.model.mesh.getFaces(
                        face_type, face_node_tags
                    )

                    # How many faces of this face_type do we have? Prisms
                    # and pyramids are special since they have faces of
                    # different types.
                    if element_type == 6:
                        num_local_faces = num_faces_prism[face_type]
                    elif element_type == 7:
                        num_local_faces = num_faces_pyramid[face_type]
                    else:
                        num_local_faces = num_faces_per_element[element_type]

                    # Store the face tags (face numbers) for each cell
                    # number (element tag)
                    for i in range(len(face_tags)):
                        element_tag = element_tags[i // num_local_faces]
                        if not element_tag in cells_2_faces:
                            cells_2_faces[element_tag] = [face_tags[i]]
                        else:
                            cells_2_faces[element_tag].append(face_tags[i])

                    # How many nodes do we have on this face type?
                    num_local_nodes = num_nodes_per_face_type[face_type]

                    # Store the node tags for the faces and vice versa
                    face_node_tags_table = reshape(face_node_tags, num_local_nodes)
                    for i, face_tag in enumerate(face_tags):
                        if not face_tag in faces_2_nodes:
                            faces_2_nodes[face_tag] = list(face_node_tags_table[i, :])
                            nodes_2_faces[
                                frozenset(face_node_tags_table[i, :])
                            ] = face_tag

    # Represent face labels by face tag instead of nodes (this could be a
    # list)
    face_markers_by_tags = {}
    for nodes, tag in nodes_2_faces.items():
        if nodes in face_markers_by_nodes:
            face_markers_by_tags[tag] = face_markers_by_nodes[nodes]
        else:
            face_markers_by_tags[tag] = default_face_marker

    # Check cell data
    print("Check data structures")
    assert len(cells_2_faces.keys()) == num_cells
    assert max(cells_2_faces.keys()) - min(cells_2_faces.keys()) + 1 == num_cells

    found_faces = []
    for f in cells_2_faces.values():
        found_faces.append(f)
    found_faces = flatten(found_faces)
    found_faces = np.unique(found_faces)
    assert min(found_faces) == 1
    assert max(found_faces) == num_faces

    assert set(cell_markers_by_tag.keys()).issubset(cells_2_faces)

    # Check face data
    assert len(faces_2_nodes.keys()) == num_faces
    assert max(faces_2_nodes.keys()) - min(faces_2_nodes.keys()) + 1 == num_faces

    found_nodes = []
    for n in faces_2_nodes.values():
        found_nodes.append(n)
    found_nodes = flatten(found_nodes)
    found_nodes = np.unique(found_nodes)
    assert min(found_nodes) == 1
    assert max(found_nodes) == num_nodes

    assert face_markers_by_tags.keys() == faces_2_nodes.keys()

    # Create MRST grid faces structure. Match order of faces_2_nodes.
    print("Create MRST data structures")
    nodes = []
    nodePos = [0]
    face_markers = []

    for k in sorted(faces_2_nodes.keys()):
        nodes.append(faces_2_nodes[k])
        nodePos.append(len(faces_2_nodes[k]))
        face_markers.append(face_markers_by_tags[k])

    nodes = np.array(flatten(nodes))
    nodePos = np.cumsum(nodePos)
    nodePos += 1
    face_markers = np.array(face_markers)

    mrst_faces = {
        "num": float(num_faces),
        "nodePos": transpose(nodePos).astype(float),
        "nodes": transpose(nodes).astype(float),
        "tags": transpose(face_markers).astype(float),
    }

    # Create MRST grid cells structure. Match order of cells_2_faces.
    faces = []
    facePos = [0]
    for k in sorted(cells_2_faces.keys()):
        faces.append(cells_2_faces[k])
        facePos.append(len(cells_2_faces[k]))
    faces = np.array(flatten(faces))
    facePos = np.cumsum(facePos)
    facePos += 1

    cell_markers = np.full(num_cells, default_cell_marker)
    for k, label in enumerate(cell_markers_by_tag.values()):
        cell_markers[k] = label

    mrst_cells = {
        "num": float(num_cells),
        "facePos": transpose(facePos).astype(float),
        "faces": transpose(faces).astype(float),
        "tags": transpose(cell_markers).astype(float),
    }

else:
    RuntimeError("Cannot handle dimension =", gdim)


# Grid
mrst_type = np.asarray(["Primordial Soup", "gmsh"], dtype=object)
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

print("Save to file", matfile)
savemat(matfile, G)

gmsh.finalize()
