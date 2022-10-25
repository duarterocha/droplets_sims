from pyoomph import *
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.advection_diffusion import *
from pyoomph.equations.stokes_stream_func import *
from pyoomph.equations.ALE import *
from pyoomph.equations.poisson import *
from stationary_model.problem_def.equations import *
from pyoomph.expressions.units import * # units
from pyoomph.utils.dropgeom import DropletGeometry
from stationary_model.problem_def.mesh import Mesh
from pyoomph.meshes.remesher import * # to remesh at large distortions
from stationary_model.problem_def.plot import Plotter

'''
Define and solve problem
'''


class DropletTempProblem(Problem):

    def __init__(self, surfactants = True):
        super(DropletTempProblem, self).__init__()

        # Droplet Properties
        self._param_Ra = self.get_global_parameter("Ra") # Rayleigh
        self._param_Ma = self.get_global_parameter("Ma") # Marangoni
        self.param_Strength = self.get_global_parameter("Strength") # Surfactancts Strength

        # Evaporation rate
        self.evap_rate = - dot(grad(var("c",domain="gas")),var("normal"))

        # Mesh Properties
        self.radius = 1  # radius of droplet that defines the mesh
        self.contact_angle = 120 * degree # contact angle
        self.domain_size = 1  # size of domain
        self.resolution = 0.025  # resolution of mesh

        # Slip length
        self.sliplength = 1

        # Pinned contact line or constant contact angle
        self.pinned_contact_line = True

        # Surfactants
        self.surfactants_bool = surfactants
        self.average_amount_surfactants = 1

        # add the plotter
        self.plotter = Plotter(self)

    def set_Ra(self,value):
        self._param_Ra.value=value

    def get_Ra(self,symbolic=False):
        if symbolic:
            return self._param_Ra.get_symbol() # Symbolic for expressions
        else:
            return self._param_Ra.value # current value

    def set_Ma(self,value):
        self._param_Ma.value=value

    def get_Ma(self,symbolic=False):
        if symbolic:
            return self._param_Ma.get_symbol() # Symbolic for expressions
        else:
            return self._param_Ma.value # current value

    def set_Strength(self,value):
        self.param_Strength.value=value

    def get_Strength(self,symbolic=False):
        if symbolic:
            return self.param_Strength.get_symbol() # Symbolic for expressions
        else:
            return self.param_Strength.value # current value

    def define_problem(self):
        # Changing to an axisymmetric coordinate system
        self.set_coordinate_system(axisymmetric)

        # Insert mesh
        geom = DropletGeometry(volume=1, contact_angle=self.contact_angle)
        mesh = Mesh(radius=geom.curv_radius, contact_angle=geom.contact_angle)
        mesh.default_resolution = self.resolution
        mesh.mesh_mode = "tris"
        mesh.remesher = Remesher2d(mesh)  # add remeshing possibility
        self.add_mesh_template(mesh)

        '''Droplet'''
        d_eqs = MeshFileOutput()

        # Mesh motion
        d_eqs += PseudoElasticMesh()

        # Set equations
        d_eqs += NavierStokesEquations(bulkforce=self.get_Ra(symbolic=True)*var("T")*vector(0,-1)*-1)  # Stokes Equation
        d_eqs += AdvectionDiffusionEquations(fieldnames="T", space="C2")  # Advection-Diffusion Equation
        d_eqs += NavierStokesFreeSurface(surface_tension=100-self.get_Ma(symbolic=True) * var("T") - self.get_Strength(symbolic=True)*var("Gamma"), mass_transfer_rate= self.evap_rate) @ "interface" # Free interface

        # Boundary Conditions
        d_eqs += DirichletBC(mesh_x=0, velocity_x=0) @ "droplet_axis"  # No flow through axis
        d_eqs += NeumannBC(T=self.evap_rate) @ "interface"  # Temperature flux through boundary
        d_eqs += ConnectMeshAtInterface() @ "interface" # Connect moving liquid and gas meshes at interface
        d_eqs += DirichletBC(T=0) @ "droplet_surface"  # Fixed temperature at contact line
        d_eqs += DirichletBC(mesh_y=0, velocity_y=0) @ "droplet_surface"  # Allow slip but fix mesh at substrate
        d_eqs += NavierStokesSlipLength(sliplength=self.sliplength)  @ "droplet_surface"# Navier-Slip condition

        # Different contact line dynamics
        if self.pinned_contact_line:  # if pinned
            # Pinned contact line means mesh_x is fixed.
            # We enforce partial_t(mesh_x)=0 by adjusting the radial velocity at the contact line
            cl_constraint = partial_t(var("mesh_x")) - 0
            d_eqs += EnforcedBC(velocity_x=cl_constraint) @ "interface/droplet_surface"
        else:
            d_eqs += NavierStokesContactAngle(contact_angle=self.contact_angle) @ "interface/droplet_surface"  # and constant contact angle

        # Add surfactants
        if self.surfactants_bool:
            lagr_mult_eqs = GlobalLagrangeMultiplier(lagrange_multiplier=0)  # Add lagrange multiplier for surfactants
            self.add_equations(lagr_mult_eqs @ "globals")
            lagr_mult = var("lagrange_multiplier", domain="globals")
            d_eqs += SurfactantTransportEquation(lagr_mult, self.average_amount_surfactants) @ "interface"

        # Plotting
        d_eqs += LocalExpressions(evap_rate=self.evap_rate) @ "interface"  # Output file
        d_eqs += StreamFunctionFromVelocity()  # Stream functions
        for boundary in ["droplet_surface", "droplet_axis", "interface"]:
            d_eqs += DirichletBC(streamfunc=0) @ boundary # No outwards flow at boundaries

        # Initial conditions
        d_eqs += InitialCondition(T=0)  # Droplet at substrate temperature

        '''Gas Phase'''
        # Gas phase is treated in quasi-stationary limit, only considered diffusive vapour transport

        g_eqs = MeshFileOutput()

        # Mesh motion
        g_eqs += PseudoElasticMesh()

        # Set equations
        g_eqs += PoissonEquation(name="c", space="C2", source=None)  # Laplace Equation, same as diffusion equation here
        # Empty dummy equations to allow for plotting the gas substrate. Without it, we cannot plot it
        g_eqs += Equations() @ "gas_surface"

        # Boundary Conditions
        g_eqs += PoissonFarFieldMonopoleCondition(far_value=0) @ "gas_infinity"
        g_eqs += DirichletBC(c=1) @ "interface" # Concentration at interface
        g_eqs += DirichletBC(mesh_x=0) @ "gas_axis"  # Fixed mesh coordinates at the boundaries
        g_eqs += DirichletBC(mesh_y=0) @ "gas_surface"
        g_eqs += DirichletBC(mesh_x=True, mesh_y=True) @ "gas_infinity" # Keep initial value for curved interface


        # Initial conditions
        g_eqs += InitialCondition(c=0) # Initial vapor concentration in the gas phase initially

        # Control remeshing
        d_eqs += RemeshWhen(RemeshingOptions(max_expansion=1.2, min_expansion=0.8))
        g_eqs += RemeshWhen(RemeshingOptions(max_expansion=1.2, min_expansion=0.8))

        '''Integral Observables'''
        # Integrate to find area and positive values of stream function
        d_eqs += IntegralObservables(_cross_area=1,_pos_stream_area=heaviside(var("streamfunc")),_coordinante_system=cartesian)
        # Find positive stream function areas
        d_eqs += IntegralObservables(positive_stream_area_fraction=lambda _pos_stream_area, _cross_area: _pos_stream_area / _cross_area)
        # Integrate the square of the velocity
        d_eqs += IntegralObservables(_v_sqr_int=dot(var("velocity"), var("velocity")), _volume=1)
        # Then divide it by the volume and take the square root => average velocity magnitude
        d_eqs += IntegralObservables(avg_velo=lambda _v_sqr_int, _volume: square_root(_v_sqr_int / _volume))
        # Volume as integral observable
        d_eqs += IntegralObservables(volume=1)
        d_eqs += IntegralObservableOutput(filename="observables",
                                         first_column=["time", "Strength", "Ma", "Ra"])  # output to write it to file
        # Add equations
        self.add_equations(d_eqs @ "droplet" + g_eqs @ "gas")