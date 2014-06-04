# This requires CGAL mesher applied to series of surfaces
#     go-smart-mesher_cgal --centre_radius 0.0001  --output_prefix ire \
#       --nearfield 0.001 --farfield 0.01 --granularity .0010 \
#       --output_gmsh --centre -0.0108 0.0174 -0.004 --zones tumour1.vtp:1 \
#       --needles needle1.vtp:1 needle2.vtp:2 needle3.vtp:3 needle4.vtp:4 \
#       needle5.vtp:5 needle6.vtp:6 --vessels vessel1.vtp:7 vessel2.vtp:8 \
#       vessel3.vtp:9 --extent exterior.vtp:10 --tissueid 0

# Questions for Ljubljana marked QLJ
import dolfin as d
from dolfin import dx

import itertools as it

# Define simulation variables
needles = {
    "positive": [1, 2, 3, 4],
    "negative": [5, 6]
}
voltages = [1300, 1500, 1300, 1900, 1300, 1300, 1300, 1900, 1300]

needle_pairs = it.product(needles["positive"], needles["negative"])

#QLJ: Why is there this final pair?
needle_pairs.append([5, 6])

# Load mesh and boundaries
mesh = d.Mesh("ire.xml")
patches = d.MeshFunction("size_t", mesh, "ire_facet_region.xml")

V = d.FunctionSpace(mesh, "Lagrange", 1)

u0 = d.Constant(0.0)
u5000 = d.Constant(5000.0)

bc1 = d.DirichletBC(V, u0, patches, 8)
bc2 = d.DirichletBC(V, u5000, patches, 9)

u = d.TrialFunction(V)
v = d.TestFunction(V)
f = d.Expression("sin(x[0])*cos(x[1])")
#g = Source()
g = d.Constant(0)

a = d.inner(d.nabla_grad(u), d.nabla_grad(v))*dx
L = g*v*dx

u = d.Function(V)
d.solve(a == L, u, bcs=[bc1, bc2])

d.plot(u)

file = d.File("poisson.pvd")
file << u

F = d.project(f, V=V)
file = d.File("poisson-solution.pvd")
file << F

F = d.project(g, V=V)
file = d.File("poisson-source.pvd")
file << F

error = d.Function(V)
error.vector()[:] = u.vector() - F.vector()
file = d.File("poisson-error.pvd")
file << error
