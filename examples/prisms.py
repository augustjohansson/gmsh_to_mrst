from pathlib import Path
import gmsh

gmsh.initialize()
model_name = Path(__file__).stem
gmsh.model.add(model_name)

L = 1000
h = L
x0 = 0
y0 = 0
x1 = L
y1 = L
p0 = gmsh.model.geo.addPoint(x0, y0, 0, h)
p1 = gmsh.model.geo.addPoint(x1, y0, 0, h)
p2 = gmsh.model.geo.addPoint(x1, y1, 0, h)
p3 = gmsh.model.geo.addPoint(x0, y1, 0, h)
l0 = gmsh.model.geo.addLine(p0, p1)
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p0)
pp = [p0, p1, p2, p3]
ll = [l0, l1, l2, l3]
cl = gmsh.model.geo.addCurveLoop(ll)
domain = gmsh.model.geo.addPlaneSurface([cl])
gmsh.model.geo.extrude([(2, domain)], 0, 0, 100, numElements=[1], recombine=True)
gmsh.model.geo.synchronize()

gdim = 3
gmsh.model.mesh.generate(3)
for dim in range(0, gdim + 1):
    for dt in gmsh.model.getEntities(dim):
        gmsh.model.addPhysicalGroup(dim, [dt[1]], dt[1])

gmsh.write(model_name + ".msh")
gmsh.write(model_name + ".mesh")
gmsh.write(model_name + ".m")
gmsh.finalize()
