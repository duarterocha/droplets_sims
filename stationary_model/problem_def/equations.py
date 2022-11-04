from pyoomph.generic import InterfaceEquations
from pyoomph.generic import Equations
from pyoomph.expressions.units import *
from pyoomph.expressions import *


class SurfactantTransportEquation(InterfaceEquations):
    def __init__(self, lagr_mult, average_value, surfactant_diffusivity=1):
        super(SurfactantTransportEquation, self).__init__()
        self.D = surfactant_diffusivity  # diffusivity
        self.lagr_mult = lagr_mult
        self.average_value = average_value

    def define_fields(self):
        self.define_scalar_field("Gamma", "C2")

    def define_residuals(self):
        u = var("velocity")  # velocity at the interface
        G, Gtest = var_and_test("Gamma")
        l, ltest = self.lagr_mult, testfunction(self.lagr_mult)
        self.add_residual(weak(partial_t(G) + div(u * G), Gtest))
        self.add_residual(weak(self.D * grad(G), grad(Gtest)))
        self.add_residual(weak(G - self.average_value, ltest) + weak(l, Gtest))


class StaticDropletInterface(Equations):
    def __init__(self, theta_dot, *, sigma = 1, contact_line_radius=1, contact_angle=120 * degree, evap_rate=1, evap_factor=1):
        super(StaticDropletInterface, self).__init__()
        self.sigma = sigma
        self.evap_rate = evap_rate
        self.evap_factor = evap_factor
        self.contact_angle = contact_angle
        self.contact_line_radius = contact_line_radius
        self.theta_dot = theta_dot

    def define_fields(self):
        self.define_scalar_field("kinbc", "C2")  # Kinematic BC lagrange multiplier
        self.define_scalar_field("un_motion", "C2")  # Motion of the droplet due to evaporation

    def define_residuals(self):

        # Define variables
        l, ltest = var_and_test("kinbc")
        u, utest = var_and_test("velocity")
        theta_dot, theta_dot_test = self.theta_dot, testfunction(self.theta_dot)
        un_motion, un_motion_test = var_and_test("un_motion")
        n = var("normal")

        # Interface motion via toroidal coordinates
        r = var("coordinate_x")
        z = var("coordinate_y")
        d1_sqr = (r + self.contact_line_radius) ** 2 + z ** 2
        d2_sqr = (r - self.contact_line_radius) ** 2 + z ** 2
        sigma_toro = pi - self.contact_angle
        tau_toro = subexpression(log((d1_sqr / d2_sqr) ** (1 / 2)))
        inter_geom_factor = self.contact_line_radius / (cosh(tau_toro) - cos(sigma_toro))

        # Add the residual for normal velocity factor
        self.add_residual(weak((inter_geom_factor - un_motion), un_motion_test))  # Normal velocity factor

        # Coupled Lagrange multipliers: theta_dot and lambda, L = lambda*(dot(u,n) - self.evap_rate*self.evap_factor
        # - self.theta_dot*un_motion))
        self.add_residual(weak(dot(u, n) - (self.evap_factor * self.evap_rate + theta_dot * un_motion), ltest))  # dL/dl
        self.add_residual(weak(l, dot(n, utest)))  # dL/du
        self.add_residual(weak(-l * un_motion, theta_dot_test))  # dL/dtheta_dot

        # Marangoni term
        self.add_residual(weak(self.sigma,  div(utest)))


