from pyoomph import *
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.advection_diffusion import *
from pyoomph.equations.stokes_stream_func import *
from pyoomph.equations.ALE import *
from pyoomph.equations.poisson import *
from pyoomph.expressions.units import *  # units
from pyoomph.utils.dropgeom import DropletGeometry
from pyoomph.meshes.remesher import *  # to remesh at large distortions
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
        self.thermal_diffusivity = 1.4558e-7 * meter ** 2 / second  # Droplet diffusivity
        self.conductivity = 0.598 * kilogram * meter / (second ** 3 * kelvin)  # Droplet thermal conductivity
        self.g = 9.81 * vector(0, -1) * meter / second**2  # Gravity acceleration
        self.droplet_density = 998 * kilogram / meter ** 3  # Density of droplet
        self.Ra_T = 0.001 * kilogram / (meter**3 * kelvin)
        self.droplet_viscosity = 1 * milli * pascal * second  # Viscosity of droplet
        self.base_radius = 1 * milli * meter  # Base radius of droplet
        self.substrate_temperature = (20 + 273.15) * kelvin  # Isothermal substrate temperature

        # Gas properties
        self.vap_diffusivity = 25.55e-6 * meter ** 2 / second  # Gas diffusivity
        self.c_sat = 0.0173 * kilogram / meter ** 3  # Saturated partial density of vapor
        self.c_infty = 0.5 * self.c_sat  # Partial vapor density far away

        # Interface properties
        self.contact_angle = 120 * degree  # Initial contact angle
        self.surface_tension = 72 * milli * newton / meter  # Surface tension
        self.Ma_T = 0.001 * milli * newton / meter / kelvin  # Factor for change in surface tension with temperature
        self.Ma_G = 0.001 * milli * newton / meter  # Factor for change in surface tension with surfactants action
        self.average_G = 0.5  # Fix average amount of surfactant concentration
        self.surfactants_diffusivity = 1e-6 * meter ** 2 / second  # Surfactants diffusivity
        self.latent_heat = 2260 * kilo * joule / kilogram  # Latent heat for evaporation
        self.evap_rate = - self.vap_diffusivity * dot(grad(var("c", domain="gas")), var("normal"))  # Evaporation rate
        self.pinned_contact_line = True  # Pinned contact line or constant contact angle

        # Mesh properties
        self.resolution = 0.05 * milli * meter  # Resolution of mesh
        self.sliplength = 1 * micro * meter  # Slip length

        # Add the plotter
        self.plotter = Plotter(self)

    def define_problem(self):
        # Settings: Axisymmetric and typical scales
        self.set_coordinate_system(axisymmetric)
        self.set_scaling(temporal=1 * second, spatial=self.base_radius)
        self.set_scaling(velocity=scale_factor("spatial") / scale_factor("temporal"))
        self.set_scaling(pressure=self.surface_tension / scale_factor("spatial"))
        self.set_scaling(T=self.substrate_temperature)
        self.set_scaling(c=self.c_sat)

        # Insert mesh
        geom = DropletGeometry(base_radius=self.base_radius, contact_angle=self.contact_angle)
        mesh = Mesh(radius=geom.curv_radius, contact_angle=geom.contact_angle)
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
        mass_density = self.droplet_density - var("T") * self.Ra_T
        d_eqs += NavierStokesEquations(mass_density=mass_density, dynamic_viscosity=self.droplet_viscosity, gravity=self.g, boussinesq=True) # Navier Stokes equations
        d_eqs += AdvectionDiffusionEquations(diffusivity=self.thermal_diffusivity, fieldnames="T", space="C2")  # Advection-Diffusion equation
        d_eqs += DirichletBC(mesh_y=0, velocity_y=0) @ "droplet_surface"  # Allow slip but fix mesh at substrate
        d_eqs += NavierStokesSlipLength(sliplength=self.sliplength) @ "droplet_surface"  # Navier-Slip condition
        d_eqs += DirichletBC(T=self.substrate_temperature) @ "droplet_surface"  # Fixed temperature at contact line
        d_eqs += DirichletBC(mesh_x=0, velocity_x=0) @ "droplet_axis"  # Fix mesh at axisymmetry and no normal flow
        d_eqs += InitialCondition(T=self.substrate_temperature)  # Droplet at substrate temperature
        d_eqs += RemeshWhen(RemeshingOptions(max_expansion=1.2, min_expansion=0.8))

        '''Interface'''

        lagr_mult_eqs = GlobalLagrangeMultiplier(lagrange_multiplier=0)  # Add lagrange multiplier for surfactants
        lagr_mult_eqs += Scaling(lagrange_multiplier=scale_factor('Gamma')/scale_factor('temporal'))
        self.add_equations(lagr_mult_eqs @ "globals")  # Add to global equations
        surface_tension = self.surface_tension - self.Ma_T * var("T") - self.Ma_G * var("Gamma")
        d_eqs += SurfactantTransportEquation(var("lagrange_multiplier", domain="globals"), self.average_G,
                                             surfactant_diffusivity=self.surfactants_diffusivity) @ "interface"  # Transport equation for surfactants concentration
        d_eqs += NavierStokesFreeSurface(surface_tension=surface_tension, mass_transfer_rate=self.evap_rate) @ "interface"  # Kinematic boundary condition, Laplace pressure and Marangoni stress
        g_eqs += DirichletBC(c=self.c_sat) @ "interface"
        g_eqs += ConnectMeshAtInterface() @ "interface"  # Connect moving liquid and gas meshes at interface
        d_eqs += ConnectMeshAtInterface() @ "interface"  # Connect moving liquid and gas meshes at interface
        d_eqs += NeumannBC(T=self.evap_rate * self.latent_heat * self.thermal_diffusivity / self.conductivity) @ "interface"  # Temperature flux through boundary

        # Different contact line dynamics
        if self.pinned_contact_line:  # if pinned
            # Pinned contact line means mesh_x is fixed
            d_eqs += EnforcedBC(velocity_x=partial_t(var("mesh_x")) - 0) @ "interface/droplet_surface"
        else:
            d_eqs += NavierStokesContactAngle(contact_angle=self.contact_angle) @ "interface/droplet_surface"  # Constant contact angle

        '''Plotting'''

        d_eqs += LocalExpressions(evap_rate=self.evap_rate) @ "interface"  # Output file
        d_eqs += StreamFunctionFromVelocity()  # Stream functions
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
        # Volume as integral observable
        d_eqs += IntegralObservables(volume=1)
        d_eqs += IntegralObservableOutput(filename="observables", first_column=["time"])  # output to write it to file

        '''Add equations'''
        self.add_equations(d_eqs @ "droplet" + g_eqs @ "gas")




