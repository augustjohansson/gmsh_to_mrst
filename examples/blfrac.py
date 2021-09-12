import gmsh

gmsh.initialize()
model_name = "testbl4"
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
lcf = lc/100
pa = geo.addPoint(0.25, 0.48, 0, lcf)
pb = geo.addPoint(0.75, 0.48, 0, lcf)
pc = geo.addPoint(0.25, 0.52, 0, lcf)
pd = geo.addPoint(0.75, 0.52, 0, lcf)

# Higher dim objects
l0 = geo.addLine(p0, p1)
l1 = geo.addLine(p1, p2)
l2 = geo.addLine(p2, p3)
l3 = geo.addLine(p3, p0)
cl = geo.addCurveLoop([l0, l1, l2, l3])
surf = geo.addPlaneSurface([cl])

# Fractures
lfrac = geo.addLine(pa, pb)
lfrac2 = geo.addLine(pc, pd)
geo.synchronize()
mesh.embed(1, [lfrac, lfrac2], 2, surf)

# Boundary layer on fractures
blfracture = field.add("BoundaryLayer")
field.setNumbers(blfracture, "CurvesList", [lfrac, lfrac2])
field.setNumber(blfracture, "SizeFar", 0.1)
field.setNumber(blfracture, "Size", lcf/2)
field.setNumber(blfracture, "Thickness", 0.01)
field.setNumber(blfracture, "Ratio", 1.2)
field.setNumbers(blfracture, "FanPointsList", [pa, pb, pc, pd])
field.setAsBoundaryLayer(blfracture)

geo.synchronize()
mesh.setRecombine(2, surf)
mesh.generate(2)
gmsh.write(model_name + ".mesh")
gmsh.write(model_name + ".msh")
gmsh.finalize()

from subprocess import call
call(["python3", "gmsh_to_mrst.py", model_name + ".msh"])
