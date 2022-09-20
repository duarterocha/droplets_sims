from pyoomph import *
from pyoomph.equations.poisson import *
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.advection_diffusion import *
from pyoomph.equations.stokes_stream_func import *
from scripts.equations import *
from pyoomph.expressions.units import * # units
from pyoomph.utils.dropgeom import DropletGeometry
from scripts.mesh import Mesh
from scripts.plot import Plotter

'''
Define and solve problem
'''


class DropletTempProblem(Problem):

    def __init__(self):
        super(DropletTempProblem, self).__init__()

        # Droplet Properties
        self.Ra = self.get_global_parameter("Ra") # Rayleigh
        self.Ma = self.get_global_parameter("Ma") # Marangoni

        # Evaporation rate
        self.evap_rate = - dot(grad(var("c",domain="gas")),var("normal"))

        # Mesh Properties
        self.radius = 1  # radius of droplet that defines the mesh
        self.contact_angle = 120 * degree # contact angle
        self.domain_size = 1  # size of domain
        self.resolution = 0.05  # resolution of mesh

        # add the plotter
        self.plotter = Plotter(self)

    def define_problem(self):
        # Changing to an axisymmetric coordinate system
        self.set_coordinate_system(axisymmetric)

        # Insert mesh
        geom = DropletGeometry(volume=1, contact_angle=self.contact_angle)
        mesh = Mesh(radius=geom.curv_radius, contact_angle=geom.contact_angle)
        mesh.default_resolution = self.resolution
        self.add_mesh_template(mesh)

        '''Droplet'''
        d_eqs = MeshFileOutput()

        # Set equations
        d_eqs += StokesEquations(bulkforce=self.Ra.get_symbol()*var("T")*vector(0,-1)*-1)  # Stokes Equation
        d_eqs += AdvectionDiffusionEquations(diffusivity=1, fieldnames="T", space="C2")  # Advection-Diffusion Equation
        d_eqs += NavierStokesFreeSurface(surface_tension=-self.Ma.get_symbol() * var("T"), static_interface=True) @ "interface" # Free interface

        # Boundary Conditions
        d_eqs += DirichletBC(velocity_x=0, velocity_y=0) @ "droplet_surface"  # No slip at substrate
        d_eqs += DirichletBC(velocity_x=0) @ "droplet_axis"  # No flow through axis
        d_eqs += DirichletBC(T=0) @ "droplet_surface"  # Fixed temperature at contact line
        d_eqs += DirichletBC(pressure=0) @ "droplet_surface/interface"  # Offset pressure
        d_eqs += NeumannBC(T=self.evap_rate) @ "interface"  # Temperature flux through boundary

        # Plotting
        d_eqs += LocalExpressions(evap_rate=self.evap_rate) @ "interface"  # Output file
        d_eqs += StreamFunctionFromVelocity()  # Stream functions
        for boundary in ["droplet_surface", "droplet_axis", "interface"]:
            d_eqs += DirichletBC(streamfunc=0) @ boundary # No outwards flow at boundaries

        '''Gas Phase'''
        g_eqs = MeshFileOutput()

        # Set equations
        g_eqs += PoissonEquation(name="c", space="C2", source=None)  # Laplace Equation
        # Empty dummy equations to allow for plotting the gas substrate. Without it, we cannot plot it
        g_eqs += Equations() @ "gas_surface"

        # Boundary Conditions
        g_eqs += PoissonFarFieldMonopoleCondition(far_value=0) @ "gas_infinity"
        g_eqs += DirichletBC(c=1) @ "interface" # Concentration at interface

        '''Integral Observables'''
        # Integrate to find area and positive values of stream function
        d_eqs += IntegralObservables(_cross_area=1,_pos_stream_area=heaviside(var("streamfunc")),_coordinante_system=cartesian)
        # Find positive stream function areas
        d_eqs += IntegralObservables(positive_stream_area_fraction=lambda _pos_stream_area, _cross_area: _pos_stream_area / _cross_area)
        # Integrate the square of the velocity
        d_eqs += IntegralObservables(_v_sqr_int=dot(var("velocity"), var("velocity")), _volume=1)
        # Then divide it by the volume and take the square root => average velocity magnitude
        d_eqs += IntegralObservables(avg_velo=lambda _v_sqr_int, _volume: square_root(_v_sqr_int / _volume))
        d_eqs += IntegralObservableOutput(filename="observables",
                                         first_column=["Ma", "Ra"])  # output to write it to file
        # Add equations
        self.add_equations(d_eqs @ "droplet" + g_eqs @ "gas")