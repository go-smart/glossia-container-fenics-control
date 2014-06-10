"""This requires CGAL mesher applied to series of surfaces
     go-smart-mesher_cgal --centre_radius 0.0001  --output_prefix ire \
       --nearfield 0.001 --farfield 0.01 --granularity .0010 \
       --output_gmsh --centre -0.0108 0.0174 -0.004 --zones tumour1.vtp:1 \
       --needles needle1.vtp:1 needle2.vtp:2 needle3.vtp:3 needle4.vtp:4 \
       needle5.vtp:5 needle6.vtp:6 --vessels vessel1.vtp:7 vessel2.vtp:8 \
       vessel3.vtp:9 --extent exterior.vtp:10 --tissueid 0
"""

# Questions for Ljubljana marked QLJ

# Use FEniCS for Finite Element
import dolfin as d

# Useful to import the derivative separately
from dolfin import dx

# Useful numerical libraries
import numpy as N
import matplotlib.pyplot as P

# Set interactive plotting on
P.ion()

# Use a separate Python file to declare variables
import variables as v


class IREProblem:
    """class IREProblem()

    This represents a Finite Element IRE problem using a similar algorithm to that of ULJ

    """
    def __init__(self):
        self.load_mesh()
        self.setup_fe()

    def load_mesh(self):
        # Load mesh and boundaries
        mesh = d.Mesh("data/ire.xml")

        #QLJ: There are properties (e.g. cond/perm) set for the inside of the electrodes - why is that?
        self.patches = d.MeshFunction("size_t", mesh, "data/ire_facet_region.xml")
        self.subdomains = d.MeshFunction("size_t", mesh, "data/ire_physical_region.xml")

        # Define differential over subdomains
        self.dxs = d.dx[self.subdomains]

        # Turn subdomains into a Numpy array
        self.subdomains_array = N.asarray(self.subdomains.array(), dtype=N.int32)

        # Create a map from subdomain indices to tissues
        self.tissues_by_subdomain = {}
        for t in v.tissues.values():
            for j in t["indices"]:
                self.tissues_by_subdomain[j] = t

        self.mesh = mesh

    def setup_fe(self):
        # Define the relevant function spaces
        V = d.FunctionSpace(self.mesh, "Lagrange", 1)
        self.V = V

        # DG0 is useful for defining piecewise constant functions
        DV = d.FunctionSpace(self.mesh, "Discontinuous Lagrange", 0)
        self.DV = DV

        # Define test and trial functions for FE
        self.z = d.TrialFunction(self.V)
        self.w = d.TestFunction(self.V)

        # Define variables that are constant on each tissue type (over DG function space)
        self.sigma_start = self.per_tissue_constant(lambda l: self.tissues_by_subdomain[l]["sigma"][0])
        self.sigma_end = self.per_tissue_constant(lambda l: self.tissues_by_subdomain[l]["sigma"][0])
        self.relative_permittivity = self.per_tissue_constant(lambda l: self.tissues_by_subdomain[l]["relative permittivity"])

    def per_tissue_constant(self, generator):
        fefunction = d.Function(self.DV)
        generated_values = [generator(l) for l in N.unique(self.subdomains_array)]
        fefunction.vector()[:] = N.choose(self.subdomains_array, generated_values)
        return fefunction

    def get_tumour_volume(self):
        # Perhaps there is a prettier way, but integrate a unit function over the tumour tets
        one = d.Function(self.V)
        one.vector()[:] = 1
        return sum(d.assemble(one * self.dxs(i)) for i in v.tissues["tumour"]["indices"])

    def solve(self):
        z, w = (self.z, self.w)
        u0 = d.Constant(0.0)

        # Define the linear and bilinear forms
        a = d.inner(d.nabla_grad(z), d.nabla_grad(w)) * dx
        L = u0 * w * dx

        #QLJ/TODO: How do we/ do we need to implement 'floating potential' boundary condition?

        # Define useful functions
        cond = d.Function(self.DV)
        U = d.Function(self.V)

        # Initialize the max_e vector, that will store the cumulative max e values
        max_e = d.Function(self.V)
        max_e.vector()[:] = 0.0
        max_e_file = d.File("output/max_e.pvd")

        es = {}
        self.max_es = {}

        # Loop through the voltages and electrode combinations
        for voltage, pair in zip(v.voltages, v.electrode_pairs):
            print "Electrodes %d (%lf) -> %d (0)" % (pair[0], voltage, pair[1])

            # Define the Dirichlet boundary conditions on the active needles
            uV = d.Constant(voltage)
            term1_bc = d.DirichletBC(self.V, uV, self.patches, pair[0])
            term2_bc = d.DirichletBC(self.V, u0, self.patches, pair[1])

            # Solve the FE problem
            d.solve(a == L, U, bcs=[term1_bc, term2_bc])

            # Extract the electric field norm
            e = d.project(d.sqrt(d.dot(d.grad(U), d.grad(U))), self.V)

            # Re-evaluate conductivity
            self.increase_conductivity(cond, e)

            for j in range(v.max_restarts):
                # Solve again
                #QLJ: is there any reason this happens before the loop and here, rather than during another iteration?
                d.solve(a == L, U, bcs=[term1_bc, term2_bc])

                # Extract electric field norm
                e_new = d.project(d.sqrt(d.dot(d.grad(U), d.grad(U))), self.V)

                # Take the max of the new field and the established electric field
                e.vector()[:] = [max(*X) for X in zip(e.vector(), e_new.vector())]

                # Re-evaluate conductivity
                self.increase_conductivity(cond, e)

            # Store this electric field norm, for this pair, for later reference
            es[pair] = e

            # Store the max of this electric field norm and that for all previous pairs
            max_e_array = N.array([max(*X) for X in zip(max_e.vector(), e.vector())])
            max_e.vector()[:] = max_e_array

            # Create a new max_e function for storage, or it will be overwritten by the next iteration
            max_e_new = d.Function(self.V)
            max_e_new.vector()[:] = max_e_array

            # Store this max e function for the cumulative coverage curve calculation later
            self.max_es[pair] = max_e_new

            # Save the max e function to a VTU
            max_e_file << max_e

    def increase_conductivity(self, cond, e):
        #TODO: Nonlinear conductivity is not correctly implemented
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
        # Calculate preliminary relationships
        dofmap = self.DV.dofmap()
        cell_dofs = N.array([dofmap.cell_dofs(c)[0] for c in N.arange(self.mesh.num_cells()) if (self.subdomains[c] in v.tissues["tumour"]["indices"])])
        volumes = N.array([d.Cell(self.mesh, c).volume() for c in N.arange(self.mesh.num_cells()) if (self.subdomains[c] in v.tissues["tumour"]["indices"])])

        # Create a horizontal axis
        cc_haxis = N.linspace(5000, 1e5, 200)

        # Calculate the tumour volume; this is what we will compare against
        tumour_volume = self.get_tumour_volume()

        # Initialize the output_arrays vector a rescale the x to V/cm
        output_arrays = [cc_haxis / 100]

        # Loop through the electrode pairs
        for pair in v.electrode_pairs:
            # Project the max e values for this pair to DG0 - this forces an evaluation of the function at the mid-point of each tet, DG0's only DOF
            e_dg = d.project(self.max_es[pair], self.DV)

            # Calculate the "max e" contribution for each cell
            contributor = N.vectorize(lambda c: e_dg.vector()[c])
            contributions = contributor(cell_dofs)

            # Sum the tet volumes for tets with a midpoint value greater than x, looping over x as e-norm thresholds (also scale to tumour volume)
            elim = N.vectorize(lambda x: volumes[contributions > x].sum() / tumour_volume)
            output_arrays.append(elim(cc_haxis))

        # Compile into a convenient array
        output = N.array(zip(*output_arrays))

        # Output cumulative coverage curves as CSV
        N.savetxt('output/coverage_curves.csv', output)

        # Plot the coverage curves
        for pair, a in zip(v.electrode_pairs, output_arrays[1:]):
            P.plot(output_arrays[0], a, label="%d - %d" % pair)

        # Draw the plot
        P.draw()

        # Show a legend for the plot
        P.legend()

        # Display the plot
        P.show(block=True)
