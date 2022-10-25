from pyoomph.generic import InterfaceEquations
from pyoomph.expressions import *


class SurfactantTransportEquation(InterfaceEquations):
    def __init__(self, lagr_mult, average_value, surfactant_diffusivity=1):
        super(SurfactantTransportEquation , self).__init__()
        self.D=surfactant_diffusivity # diffusivity
        self.lagr_mult = lagr_mult
        self.average_value = average_value

    def define_fields(self):
        self.define_scalar_field("Gamma","C2")

    def define_residuals(self):
        u=var("velocity") # velocity at the interface
        G,Gtest=var_and_test("Gamma")
        l, ltest = self.lagr_mult, testfunction(self.lagr_mult)
        self.add_residual(weak(partial_t(G)+div(u*G),Gtest))
        self.add_residual(weak(self.D*grad(G), grad(Gtest)))
        self.add_residual(weak(G - self.average_value, ltest) + weak(l, Gtest))