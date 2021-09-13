import gmsh

gmsh.initialize()
model_name = "fractures"
gmsh.model.add(model_name)
geo = gmsh.model.geo
mesh = gmsh.model.mesh
field = gmsh.model.mesh.field

# Corner points
lc = 0.05
p0 = geo.addPoint(0, 0, 0, lc)
p1 = geo.addPoint(1, 0, 0, lc)
p2 = geo.addPoint(1, 1, 0, lc)
p3 = geo.addPoint(0, 1, 0, lc)

# Fracture points
lcf = lc/2
x=[   0.751267059305653,
      0.255095115459269,
      0.505957051665142,
      0.699076722656686,
      0.890903252535798,
      0.959291425205444,
      0.547215529963803,
      0.138624442828679,
      0.149294005559057,
      0.257508254123736,
      0.840717255983663,
      0.254282178971531,
      0.814284826068816,
      0.243524968724989,
      0.929263623187228,
      0.349983765984809,
      0.196595250431208,
      0.251083857976031,
      0.616044676146639,
      0.473288848902729,
      0.351659507062997,
      0.830828627896291,
      0.585264091152724,
      0.549723608291140,
      0.917193663829810,
      0.285839018820374,
      0.757200229110721,
      0.753729094278495,
      0.380445846975357,
      0.567821640725221,
      0.075854289563064,
      0.053950118666607,
      0.530797553008973,
      0.779167230102011,
      0.934010684229183,
      0.129906208473730,
      0.568823660872193,
      0.469390641058206,
      0.011902069501241,
   0.337122644398882]
y=[
    0.162182308193243,
    0.794284540683907,
    0.311215042044805,
    0.528533135506213,
    0.165648729499781,
    0.601981941401637,
    0.262971284540144,
    0.654079098476782,
    0.689214503140008,
    0.748151592823709,
    0.450541598502498,
    0.083821377996933,
    0.228976968716819,
    0.913337361501670,
    0.152378018969223,
    0.825816977489547,
    0.538342435260057,
    0.996134716626885,
    0.078175528753184,
    0.442678269775446,
    0.106652770180584,
    0.961898080855054,
    0.004634224134067,
    0.774910464711502,
    0.817303220653433,
    0.868694705363510,
    0.084435845510910,
    0.399782649098896,
    0.259870402850654,
    0.800068480224308,
    0.431413827463545,
    0.910647594429523,
    0.181847028302852,
    0.263802916521990,
    0.145538980384717,
    0.136068558708664,
    0.869292207640089,
    0.579704587365570,
    0.549860201836332,
   0.144954798223727];

p = []
for k in range(len(x)):
    p.append(geo.addPoint(x[k],y[k],0,lcf))

# Higher dim objects
l0 = geo.addLine(p0, p1)
l1 = geo.addLine(p1, p2)
l2 = geo.addLine(p2, p3)
l3 = geo.addLine(p3, p0)
cl = geo.addCurveLoop([l0, l1, l2, l3])
surf = geo.addPlaneSurface([cl])

# Fractures
lfrac = []
for k in range(len(x)//2):
    lfrac.append(geo.addLine(p[2*k],p[2*k+1]))
geo.synchronize()
mesh.embed(1, lfrac, 2, surf)

# # Boundary layer on fractures
# blfracture = field.add("BoundaryLayer")
# field.setNumbers(blfracture, "CurvesList", [lfrac, lfrac2])
# field.setNumber(blfracture, "SizeFar", 0.1)
# field.setNumber(blfracture, "Size", lcf/2)
# field.setNumber(blfracture, "Thickness", 0.01)
# field.setNumber(blfracture, "Ratio", 1.2)
# field.setNumbers(blfracture, "FanPointsList", [pa, pb, pc, pd])
# field.setAsBoundaryLayer(blfracture)

#gmsh.option.setNumber("Mesh.Algorithm", 8)

geo.synchronize()
#mesh.setRecombine(2, surf)
mesh.generate(2)
gmsh.write(model_name + ".mesh")
gmsh.write(model_name + ".msh")
gmsh.finalize()

from subprocess import call
call(["python3", "../src/gmsh_to_mrst.py", model_name + ".msh"])
