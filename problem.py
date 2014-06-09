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
import variables as v

import numpy as N
import matplotlib.pyplot as P

P.ion()


class IREProblem:
    def __init__(self):
        self.load_mesh()
        self.setup_fe()

    def load_mesh(self):
        # Load mesh and boundaries
        mesh = d.Mesh("data/ire.xml")

        #QLJ: There are properties (e.g. cond/perm) set for the inside of the electrodes - why is that?
        self.patches = d.MeshFunction("size_t", mesh, "data/ire_facet_region.xml")
        self.subdomains = d.MeshFunction("size_t", mesh, "data/ire_physical_region.xml")
        self.dxs = d.dx[self.subdomains]

        self.subdomains_array = N.asarray(self.subdomains.array(), dtype=N.int32)
        self.tissues_by_subdomain = {}
        for t in v.tissues.values():
            for j in t["indices"]:
                self.tissues_by_subdomain[j] = t

        self.mesh = mesh

    def setup_fe(self):
        V = d.FunctionSpace(self.mesh, "Lagrange", 1)
        DV = d.FunctionSpace(self.mesh, "Discontinuous Lagrange", 0)
        self.V = V
        self.DV = DV

        # Assign variables constant on each tissue type
        self.per_tissue_constants = []

        self.sigma_start = self.per_tissue_constant(lambda l: self.tissues_by_subdomain[l]["sigma"][0])
        self.sigma_end = self.per_tissue_constant(lambda l: self.tissues_by_subdomain[l]["sigma"][0])
        self.relative_permittivity = self.per_tissue_constant(lambda l: self.tissues_by_subdomain[l]["relative permittivity"])

        self.z = d.TrialFunction(self.V)
        self.w = d.TestFunction(self.V)
        print self.get_tumour_volume()

    def per_tissue_constant(self, generator):
        fefunction = d.Function(self.DV)
        generated_values = [generator(l) for l in N.unique(self.subdomains_array)]
        fefunction.vector()[:] = N.choose(self.subdomains_array, generated_values)
        return fefunction

    def get_tumour_volume(self):
        one = d.Function(self.V)
        one.vector()[:] = 1
        return sum(d.assemble(one * self.dxs(i)) for i in v.tissues["tumour"]["indices"])

    def solve(self):
        z, w = (self.z, self.w)
        u0 = d.Constant(0.0)

        a = d.inner(d.nabla_grad(z), d.nabla_grad(w)) * dx
        L = u0 * w * dx

        cond = d.Function(self.DV)

        U = d.Function(self.V)

        max_e = d.Function(self.V)
        max_e.vector()[:] = 0.0
        max_e_file = d.File("output/max_e.pvd")
        es = {}
        self.max_es = {}
        for voltage, pair in zip(v.voltages, v.electrode_pairs):
            print pair
            uV = d.Constant(voltage)
            term1_bc = d.DirichletBC(self.V, uV, self.patches, pair[0])
            term2_bc = d.DirichletBC(self.V, u0, self.patches, pair[1])

            d.solve(a == L, U, bcs=[term1_bc, term2_bc])
            e = d.project(d.sqrt(d.dot(d.grad(U), d.grad(U))), self.V)
            self.increase_conductivity(cond, e)

            for j in range(v.max_restarts):
                d.solve(a == L, U, bcs=[term1_bc, term2_bc])
                e_new = d.project(d.sqrt(d.dot(d.grad(U), d.grad(U))), self.V)
                e.vector()[:] = [max(*X) for X in zip(e.vector(), e_new.vector())]
                self.increase_conductivity(cond, e)

            es[pair] = e

            max_e_array = N.array([max(*X) for X in zip(max_e.vector(), e.vector())])
            max_e.vector()[:] = max_e_array
            max_e_new = d.Function(self.V)
            max_e_new.vector()[:] = max_e_array  # d.assemble(max_e * dxs(tumour_tissue["indices"]))
            self.max_es[pair] = max_e_new
            max_e_file << max_e

    def increase_conductivity(self, cond, e):
        emax = max(e.vector())

        def cond_function(l):
            sigma_start = self.tissues_by_subdomain[l]["sigma"][0]
            sigma_end = self.tissues_by_subdomain[l]["sigma"][1]
            threshold_reversible = self.tissues_by_subdomain[l]["threshold reversible"]
            threshold_irreversible = self.tissues_by_subdomain[l]["threshold irreversible"]

            if emax < threshold_reversible:
                return sigma_start
            if emax > threshold_irreversible:
                return sigma_end

            k = (sigma_end - sigma_start) / (threshold_irreversible - threshold_reversible)
            h = sigma_start - k * threshold_reversible
            return emax * k + h

        generated_values = [cond_function(l) for l in N.unique(self.subdomains_array)]
        cond.vector()[:] = N.choose(self.subdomains_array, generated_values)

    def plot_result(self):
        dofmap = self.DV.dofmap()
        cell_dofs = N.array([dofmap.cell_dofs(c)[0] for c in N.arange(self.mesh.num_cells()) if (self.subdomains[c] in v.tissues["tumour"]["indices"])])
        volumes = N.array([d.Cell(self.mesh, c).volume() for c in N.arange(self.mesh.num_cells()) if (self.subdomains[c] in v.tissues["tumour"]["indices"])])
        cc_haxis = N.linspace(5000, 1e5, 200)

        tumour_volume = self.get_tumour_volume()

        output_arrays = [cc_haxis / 100]  # V / cm
        for pair in v.electrode_pairs:
            e_dg = d.project(self.max_es[pair], self.DV)
            contributor = N.vectorize(lambda c: e_dg.vector()[c])
            contributions = contributor(cell_dofs)
            elim = N.vectorize(lambda x: volumes[contributions > x].sum() / tumour_volume)
            output_arrays.append(elim(cc_haxis))
        output = N.array(zip(*output_arrays))

        N.savetxt('output/coverage_curves.csv', output)

        for pair, a in zip(v.electrode_pairs, output_arrays[1:]):
            P.plot(output_arrays[0], a, label="%d - %d" % pair)

        P.draw()
        P.legend()
        P.show(block=True)
