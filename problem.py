"""This requires CGAL mesher applied to series of surfaces. See readme.txt for details.
"""

from __future__ import print_function

# Use FEniCS for Finite Element
import fenics as d

# Useful to import the derivative separately
from dolfin import dx

# Useful numerical libraries
import numpy as N
import matplotlib.pyplot as P

# General tools
import os

# UFL
from ufl import *

# Set interactive plotting on
P.ion()

# Use a separate Python file to declare variables
import variables as v

input_mesh = "ire-cuboid-7"


class IREProblem:
    """class IREProblem()

    This represents a Finite Element IRE problem using a similar algorithm to that of ULJ

    """
    def __init__(self):
        self.load_mesh()
        self.setup_fe()
        self.prepare_increase_conductivity()

    def load_mesh(self):
        # Load mesh and boundaries
        mesh = d.Mesh("input/%s.xml" % input_mesh)

        # Import the binary data indicating the location of structures
        self.load_patient_data()

        self.patches = d.MeshFunction("size_t", mesh, "input/%s_facet_region.xml" % input_mesh)
        self.subdomains = d.MeshFunction("size_t", mesh, "input/%s_physical_region.xml" % input_mesh)

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

    def load_patient_data(self):
        # width = v.dim_width * v.delta_width
        # height = v.dim_height * v.delta_height
        # depth = v.dim_depth * v.delta_depth

        # self.box = d.BoxMesh(
        #     0, 0, 0, height, width, depth,
        #     v.dim_height, v.dim_width, v.dim_depth
        # )
        # self.BV = d.FunctionSpace(self.box, "CG", 1)
        # self.box_map = self.BV.dofmap().vertex_to_dof_map(self.box)

        indicators = {}
        # print(self.box.coordinates())
        for subdomain in ("liver", "vessels", "tumour"):
            values = N.empty((v.dim_height, v.dim_width, v.dim_depth), dtype='uintp')
            # values = N.empty((v.dim_height, v.dim_width, v.dim_depth, 6), dtype='uintp')
            for i in range(0, v.dim_depth):
                slice = N.loadtxt(os.path.join(
                    v.patient_data_location,
                    "patient-%s.%d.txt" % (subdomain, i + 1))
                )
                # r = (6 * i * v.dim_width * v.dim_height, 6 * (i + 1) * v.dim_width * v.dim_height)
                values[:, :, i] = slice.astype('uintp')
                # for j in range(6):
                #     values[:, :, i, j] = slice.astype('uintp')
            indicators[subdomain] = values
        self.indicators = indicators
        #     values = N.swapaxes(values, 0, 2)
        #     indicators[subdomain].set_values(values.flatten())

        #     f = d.File(subdomain + ".pvd")
        #     f << indicators[subdomain]
        #     print("Output " + subdomain)

    def interpolate_to_patient_data(self, function, indicator):
        values = N.empty((v.dim_height, v.dim_width, v.dim_depth), dtype='float')
        it = N.nditer(values, flags=['multi_index'])

        u = N.empty((1,))
        x = N.empty((3,))
        delta = (v.delta_height, v.delta_width, v.delta_depth)
        offset = (v.offset_x, v.offset_y, v.offset_z)
        while not it.finished:
            if indicator[it.multi_index] != 1:
                it.iternext()
                continue

            x[0] = it.multi_index[1] * delta[1] - offset[0]
            x[1] = it.multi_index[0] * delta[0] - offset[1]
            x[2] = it.multi_index[2] * delta[2] - offset[2]
            #print(it.multi_index)
            #print(delta)
            #print(function.function_space().mesh().coordinates()[0])
            #print(x)
            function.eval(u, x)
            values[...] = u[0]
            it.iternext()

        return values

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

        regions = self.per_tissue_constant(lambda l: l)
        region_file = d.File("output/%s-regions.pvd" % input_mesh)
        region_file << regions

    def per_tissue_constant(self, generator):
        fefunction = d.Function(self.DV)
        generated_values = dict((l, generator(l)) for l in N.unique(self.subdomains_array))
        vector = N.vectorize(generated_values.get)
        fefunction.vector()[:] = vector(self.subdomains_array)
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
        L = u0 * w * dx

        # Define useful functions
        cond = d.Function(self.DV)
        U = d.Function(self.V)

        # Initialize the max_e vector, that will store the cumulative max e values
        max_e = d.Function(self.V)
        max_e.vector()[:] = 0.0
        max_e_file = d.File("output/%s-max_e.pvd" % input_mesh)
        max_e_per_step = d.Function(self.V)
        max_e_per_step_file = d.File("output/%s-max_e_per_step.pvd" % input_mesh)

        self.es = {}
        self.max_es = {}
        fi = d.File("output/%s-cond.pvd" % input_mesh)

        potential_file = d.File("output/%s-potential.pvd" % input_mesh)

        # Loop through the voltages and electrode combinations
        for voltage, pair in zip(v.voltages, v.electrode_pairs):
            print("Electrodes %d (%lf) -> %d (0)" % (pair[0], voltage, pair[1]))

            cond = d.project(self.sigma_start, V=self.DV)

            # Define the Dirichlet boundary conditions on the active needles
            uV = d.Constant(voltage)
            term1_bc = d.DirichletBC(self.V, uV, self.patches, pair[0])
            term2_bc = d.DirichletBC(self.V, u0, self.patches, pair[1])

            e = d.Function(self.V)
            e.vector()[:] = max_e.vector()

            # Re-evaluate conductivity
            self.increase_conductivity(cond, e, pair)

            for j in range(v.max_restarts):
                # Update the bilinear form
                a = d.inner(d.nabla_grad(z), cond * d.nabla_grad(w)) * dx

                # Solve again
                print(" [solving...")
                #d.solve(a == L, U, bcs=[term1_bc, term2_bc])
                d.solve(a == L, U, bcs=[term1_bc, term2_bc])
                #d.solve(a == L, U, bcs=[term1_bc, term2_bc], solver_parameters={"linear_solver": "bicgstab", "krylov_solver": {"relative_tolerance": 1e-8}})
                print("  ....solved]")

                # Extract electric field norm
                for i in range(len(U.vector())):
                    if N.isnan(U.vector()[i]):
                        U.vector()[i] = 1e5

                e_new = d.project(d.sqrt(d.dot(d.grad(U), d.grad(U))), self.V)

                # Take the max of the new field and the established electric field
                e.vector()[:] = N.array([max(*X) for X in zip(e.vector(), e_new.vector())])

                # Re-evaluate conductivity
                fi << cond
                self.increase_conductivity(cond, e, pair)

            potential_file << U

            # Save the max e function to a VTU
            max_e_per_step.vector()[:] = e.vector()[:]
            max_e_per_step_file << max_e_per_step

            # Store this electric field norm, for this pair, for later reference
            self.es[pair] = e

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

    def prepare_increase_conductivity(self):
        def sigma_function(l, i):
            s = self.tissues_by_subdomain[l]["sigma"]
            if isinstance(s, tuple):
                return s[i]
            else:
                return s

        def threshold_function(l, i):
            s = self.tissues_by_subdomain[l]["sigma"]
            if isinstance(s, tuple):
                return self.tissues_by_subdomain[l][i]
            else:
                return 1 if i == "threshold reversible" else 0

        self.sigma_start = self.per_tissue_constant(lambda l: sigma_function(l, 0))
        self.sigma_end = self.per_tissue_constant(lambda l: sigma_function(l, 1))
        self.threshold_reversible = self.per_tissue_constant(lambda l: threshold_function(l, "threshold reversible"))
        self.threshold_irreversible = self.per_tissue_constant(lambda l: threshold_function(l, "threshold irreversible"))
        self.k = (self.sigma_end - self.sigma_start) / (self.threshold_irreversible - self.threshold_reversible)
        self.h = self.sigma_start - self.k * self.threshold_reversible

    def increase_conductivity(self, cond, e, pair):
        # Set up the three way choice function
        intermediate = e * self.k + self.h
        not_less_than = conditional(gt(e, self.threshold_irreversible), self.sigma_end, intermediate)
        cond_expression = conditional(lt(e, self.threshold_reversible), self.sigma_start, not_less_than)

        # Project this onto the function space
        cond_function = d.project(Max(cond_expression, cond), cond.function_space())
        cond.assign(cond_function)

    def plot_bitmap_result(self):
        # Create a horizontal axis
        cc_haxis = N.linspace(5000, 1e5, 200)

        # Calculate the tumour volume; this is what we will compare against
        tumour_volume = (self.indicators["tumour"] == 1).sum()

        # Initialize the output_arrays vector a rescale the x to V/cm
        output_arrays = [cc_haxis / 100]

        # Loop through the electrode pairs
        for pair in v.electrode_pairs:
            # Project the max e values for this pair to DG0 - this forces an evaluation of the function at the mid-point of each tet, DG0's only DOF
            e_dg = self.interpolate_to_patient_data(self.max_es[pair], self.indicators["tumour"])

            # Sum the tet volumes for tets with a midpoint value greater than x, looping over x as e-norm thresholds (also scale to tumour volume)
            elim = N.vectorize(lambda x: (e_dg > x).sum() / tumour_volume)
            output_arrays.append(elim(cc_haxis))

        # Compile into a convenient array
        output = N.array(zip(*output_arrays))

        # Output cumulative coverage curves as CSV
        N.savetxt('output/%s-coverage_curves_bitmap.csv' % input_mesh, output)

        # Plot the coverage curves
        for pair, a in zip(v.electrode_pairs, output_arrays[1:]):
            P.plot(output_arrays[0], a, label="%d - %d" % pair)

        # Draw the plot
        P.draw()
        P.title(r"Bitmap-based")
        P.xlabel(r"Threshold level of $|E|$ ($\mathrm{J}$)")
        P.ylabel(r"Fraction of tumour beneath level")

        # Show a legend for the plot
        P.legend(loc=3)

        # Display the plot
        P.show(block=True)

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
        N.savetxt('output/%s-coverage_curves.csv' % input_mesh, output)

        # Plot the coverage curves
        for pair, a in zip(v.electrode_pairs, output_arrays[1:]):
            P.plot(output_arrays[0], a, label="%d - %d" % pair)

        # Draw the plot
        P.draw()
        P.xlabel(r"Threshold level of $|E|$ ($\mathrm{J}$)")
        P.ylabel(r"Fraction of tumour beneath level")

        # Show a legend for the plot
        P.legend(loc=3)

        # Display the plot
        P.show(block=True)
