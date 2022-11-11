from pyoomph import *
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.advection_diffusion import *
from pyoomph.equations.stokes_stream_func import *
from pyoomph.equations.multi_component import *
from pyoomph.equations.ALE import *
from pyoomph.equations.poisson import *
from pyoomph.expressions.units import *  # units
from pyoomph.utils.dropgeom import DropletGeometry
from pyoomph.meshes.remesher import *  # to remesh at large distortions
from pyoomph.materials.default_materials import *
from temporal_model_dimensional.problem_def.equations import *
from temporal_model_dimensional.problem_def.mesh import Mesh
from temporal_model_dimensional.problem_def.plot import Plotter

'''
Define and solve problem
'''


class EvapWaterDropletProblem(Problem):

    def __init__(self):
        super(EvapWaterDropletProblem, self).__init__()

        # Droplet properties
        self.water = PureLiquidWater()
        self.water.evaluate_at_condition("dynamic_viscosity", temperature=21 * celsius)
        self.droplet_density = self.water.mass_density  # Density of droplet
        self.droplet_viscosity = self.water.dynamic_viscosity  # Viscosity of droplet
        self.specific_heat_capacity = self.water.specific_heat_capacity
        self.conductivity = self.water.thermal_conductivity  # Droplet thermal conductivity
        self.thermal_diffusivity = self.conductivity / (self.droplet_density * self.specific_heat_capacity)  # Droplet diffusivity
        self.gravity = 9.81 * meter / second**2
        self.g = self.gravity * vector(0, -1) * minimum(var("time") / (2*second), 1)  # Gravity acceleration
        self.volume = (3.5 * 10**(-9)) * meter ** 3  # Volume of droplet
        self.substrate_temperature = (21 + 273.15) * kelvin  # Isothermal substrate temperature

        # Gas properties
        self.vap_diffusivity = 25.55e-6 * meter ** 2 / second  # Gas diffusivity
        self.c_sat = 0.0173 * kilogram / meter ** 3  # Saturated partial density of vapor
        self.c_infty = 0.75 * self.c_sat  # Partial vapor density far away

        # Interface properties
        self.contact_angle = 155 * degree  # Initial contact angle
        self.surface_tension_ref = 72 * milli * newton / meter
        self.surface_tension = self.water.default_surface_tension["gas"]
        self.Ma_G = self.get_global_parameter("surfactant_Ma_G") # Surfactancts Strength
        self.initial_G = 1  # Fix average amount of surfactant concentration
        self.surfactants_diffusivity = 1e-9 * meter ** 2 / second  # Surfactants diffusivity
        self.latent_heat = self.water.latent_heat_of_evaporation
        self.evap_rate = - self.vap_diffusivity * dot(grad(var("c", domain="gas")), var("normal"))  # Evaporation rate
        self.pinned_contact_line = True  # Pinned contact line or constant contact angle

        # Mesh properties
        self.resolution = 0.025 * milli * meter  # Resolution of mesh
        self.sliplength = 1 * micro * meter  # Slip length
        self.droplet_geom = DropletGeometry(volume=self.volume, contact_angle=self.contact_angle)
        self.droplet_radius = self.droplet_geom.curv_radius

        # Add the plotter
        self.plotter = Plotter(self)

    def set_Ma_G(self,value):
        self.Ma_G.value=value

    def get_Ma_G(self,symbolic=False):
        if symbolic:
            return self.Ma_G.get_symbol() # Symbolic for expressions
        else:
            return self.Ma_G.value # current value

    def define_problem(self):
        # Settings: Axisymmetric and typical scales
        self.set_coordinate_system(axisymmetric)
        self.set_scaling(temporal=second, spatial=self.volume ** rational_num(1,3))
        self.set_scaling(velocity=scale_factor("spatial") / scale_factor("temporal"))
        self.set_scaling(temperature=self.substrate_temperature)
        self.set_scaling(pressure= self.surface_tension_ref / scale_factor("spatial"))
        self.set_scaling(c=self.c_sat)

        # Insert mesh
        mesh = Mesh(radius=self.droplet_geom.curv_radius, contact_angle=self.droplet_geom.contact_angle)
        mesh.default_resolution = self.resolution
        mesh.mesh_mode = "tris"
        mesh.remesher = Remesher2d(mesh)  # Add remeshing possibility
        self.add_mesh_template(mesh)

        '''Gas equations'''

        g_eqs = MeshFileOutput()
        g_eqs += PseudoElasticMesh()  # Mesh motion
        g_eqs += PoissonEquation(name='c', space="C2", source=None, coefficient=self.vap_diffusivity / self.c_sat * scale_factor('temporal'))  # Laplace equation for gas
        g_eqs += DirichletBC(mesh_x=0) @ "gas_axis"  # Fixed mesh coordinates at the boundaries
        g_eqs += DirichletBC(mesh_y=0) @ "gas_surface"
        g_eqs += DirichletBC(mesh_x=True, mesh_y=True) @ "gas_infinity"  # Keep initial value for curved interface
        g_eqs += PoissonFarFieldMonopoleCondition(name="c", far_value=self.c_infty, coefficient=self.vap_diffusivity / self.c_sat * scale_factor('temporal') / scale_factor('spatial')) @ "gas_infinity"  # Mimic infinite mesh
        g_eqs += InitialCondition(c=self.c_infty)  # Initial vapor concentration in the gas phase initially
        g_eqs += RemeshWhen(RemeshingOptions(max_expansion=1.2, min_expansion=0.8)) # Control remeshing

        '''Droplet equations'''

        d_eqs = MeshFileOutput()  # Insert mesh
        d_eqs += PseudoElasticMesh()  # Mesh motion
        d_eqs += NavierStokesEquations(mass_density=self.droplet_density, dynamic_viscosity=self.droplet_viscosity, gravity=self.g, boussinesq=True) # Navier Stokes equations
        d_eqs += AdvectionDiffusionEquations(diffusivity=self.thermal_diffusivity, fieldnames="temperature", space="C2")  # Advection-Diffusion equation
        d_eqs += DirichletBC(mesh_y=0, velocity_y=0) @ "droplet_surface"  # Allow slip but fix mesh at substrate
        d_eqs += NavierStokesSlipLength(sliplength=self.sliplength) @ "droplet_surface"  # Navier-Slip condition
        d_eqs += DirichletBC(temperature=self.substrate_temperature) @ "droplet_surface"  # Fixed temperature at contact line
        d_eqs += DirichletBC(mesh_x=0, velocity_x=0) @ "droplet_axis"  # Fix mesh at axisymmetry and no normal flow
        d_eqs += InitialCondition(temperature=self.substrate_temperature)  # Droplet at substrate temperature
        d_eqs += RemeshWhen(RemeshingOptions(max_expansion=1.2, min_expansion=0.8))

        '''Interface'''

        d_eqs += SurfactantTransportEquationMassTransfer(surfactant_diffusivity=self.surfactants_diffusivity) @ "interface"  # Transport equation for surfactants concentration
        d_eqs += InitialCondition(Gamma=self.initial_G) @ "interface"
        d_eqs += DirichletBC(u_p_x=0, u_p_y=0) @ "interface/droplet_surface"
        d_eqs += DirichletBC(u_p_x=0, u_p_y=0) @ "interface/droplet_axis"
        surface_tension = self.surface_tension - self.get_Ma_G(symbolic=True) * var("Gamma") / self.initial_G * 72 * milli * newton / meter
        d_eqs += NavierStokesFreeSurface(surface_tension=surface_tension, mass_transfer_rate=self.evap_rate) @ "interface"  # Kinematic boundary condition, Laplace pressure and Marangoni stress
        g_eqs += DirichletBC(c=self.c_sat) @ "interface"
        d_eqs += ConnectMeshAtInterface() @ "interface"  # Connect moving liquid and gas meshes at interface
        d_eqs += NeumannBC(temperature=self.evap_rate * self.latent_heat * self.thermal_diffusivity / self.conductivity) @ "interface"  # Temperature flux through boundary
        # Different contact line dynamics
        if self.pinned_contact_line:  # if pinned
            # Pinned contact line means mesh_x is fixed
            d_eqs += EnforcedBC(velocity_x=partial_t(var("mesh_x")) - 0) @ "interface/droplet_surface"
        else:
            d_eqs += NavierStokesContactAngle(contact_angle=self.contact_angle) @ "interface/droplet_surface"  # Constant contact angle

        '''Plotting'''

        d_eqs += LocalExpressions(evap_rate=self.evap_rate) @ "interface"  # Output file
        d_eqs += StreamFunctionFromVelocity()  # Stream functions
        d_eqs += StreamFunctionFromVelocityInterface() @ 'interface'
        for boundary in ["droplet_surface", "droplet_axis", "interface"]:
            d_eqs += DirichletBC(streamfunc=0) @ boundary  # No outwards flow at boundaries

        # Integrate to find area and positive values of stream function
        d_eqs += IntegralObservables(_cross_area=1, _pos_stream_area=heaviside(var("streamfunc")), _coordinante_system=cartesian)
        # Find positive stream function areas
        d_eqs += IntegralObservables(positive_stream_area_fraction=lambda _pos_stream_area, _cross_area: _pos_stream_area / _cross_area)
        # Integrate the square of the velocity
        d_eqs += IntegralObservables(_v_sqr_int=dot(var("velocity"), var("velocity")), _volume=1)
        # Then divide it by the volume and take the square root => average velocity magnitude
        d_eqs += IntegralObservables(avg_velo=lambda _v_sqr_int, _volume: square_root(_v_sqr_int / _volume))
        d_eqs += IntegralObservables(_coordinante_system=cartesian, temperature_int=var('temperature'), len = 1) @ "interface/droplet_axis"
        d_eqs += IntegralObservables(temperature = lambda  temperature_int, len: temperature_int/len) @ "interface/droplet_axis"
        # Volume as integral observable
        d_eqs += IntegralObservables(volume=1)
        d_eqs += IntegralObservableOutput(filename="observables", first_column=["time"])  # output to write it to file
        d_eqs += IntegralObservables(amount_surf=var("Gamma")) @ "interface"
        d_eqs += IntegralObservableOutput(filename="observables_interface", first_column=["time"]) @ "interface"

        '''Add equations'''
        self.add_equations(d_eqs @ "droplet" + g_eqs @ "gas")




