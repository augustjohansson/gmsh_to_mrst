'''
https://gitlab.onelab.info/gmsh/gmsh/-/issues/1530
'''

import gmsh

gmsh.initialize()
model_name = "testbl3"
gmsh.model.add(model_name)
geo = gmsh.model.geo
mesh = gmsh.model.mesh
field = gmsh.model.mesh.field

# Corner points
lc = 0.2
p0 = geo.addPoint(0, 0, 0, lc)
p1 = geo.addPoint(1, 0, 0, lc)
p2 = geo.addPoint(1, 1, 0, lc)
p3 = geo.addPoint(0, 1, 0, lc)

# Fracture points
lcf = lc/10
pa = geo.addPoint(0.33, 0, 0, lcf)
pb = geo.addPoint(0.33, 0.33, 0, lcf)

# Higher dim objects
l00 = geo.addLine(p0, pa)
l01 = geo.addLine(pa, p1)
l1 = geo.addLine(p1, p2)
l2 = geo.addLine(p2, p3)
l3 = geo.addLine(p3, p0)
cl = geo.addCurveLoop([l00, l01, l1, l2, l3])
surf = geo.addPlaneSurface([cl])

# Fracture
lfrac = geo.addLine(pa, pb)
geo.synchronize()
mesh.embed(1, [lfrac], 2, surf)

# # Boundary layer on boundary edge
# '''
# Without mesh.embed, the boundary layer around [l00, l01] works.

# With mesh.embed, the following error occurs
# Error   : The 1D mesh seems not to be forming a closed loop (2 boundary nodes are considered once)
# '''
# bledge = field.add("BoundaryLayer")
# field.setNumbers(bledge, "CurvesList", [l00, l01])
# field.setNumbers(bledge, "PointsList", [p0, p1])
# field.setNumber(bledge, "SizeFar", 0.1)
# field.setNumber(bledge, "Size", lc/100)
# field.setNumber(bledge, "Thickness", 0.4)
# field.setNumber(bledge, "Ratio", 1.4)
# field.setNumber(bledge, "BetaLaw", 0.5)
# field.setAsBoundaryLayer(bledge)

# # Boundary layer on fracture
# blfracture = field.add("BoundaryLayer")
# field.setNumbers(blfracture, "CurvesList", [lfrac])
# field.setNumbers(blfracture, "PointsList", [pa, pb])
# field.setNumber(blfracture, "SizeFar", 0.1)
# field.setNumber(blfracture, "Size", lcf/10)
# field.setNumber(blfracture, "Thickness", 0.2)
# field.setNumber(blfracture, "Ratio", 1.4)
# field.setNumber(blfracture, "BetaLaw", 0.5)
# field.setAsBoundaryLayer(blfracture)


geo.synchronize()
mesh.setRecombine(2, surf)
mesh.generate(2)
gmsh.write(model_name + ".msh")
gmsh.finalize()
