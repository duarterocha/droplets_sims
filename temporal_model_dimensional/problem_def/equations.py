from pyoomph.generic import InterfaceEquations
from pyoomph.expressions import *


class SurfactantTransportEquation(InterfaceEquations):
    def __init__(self, lagr_mult, average_value, surfactant_diffusivity=1):
        super(SurfactantTransportEquation , self).__init__()
        self.D=surfactant_diffusivity # diffusivity
        self.lagr_mult = lagr_mult
        self.average_value = average_value

    def define_fields(self):
        self.define_scalar_field("Gamma","C2", testscale=scale_factor('temporal')/scale_factor('Gamma'))

    def define_residuals(self):
        u=var("velocity") # velocity at the interface
        G,Gtest=var_and_test("Gamma")
        l, ltest = self.lagr_mult, testfunction(self.lagr_mult)
        self.add_residual(weak(partial_t(G, ALE="auto")+div(u*G),Gtest))
        self.add_residual(weak(self.D*grad(G), grad(Gtest)))
        self.add_residual(weak(G - self.average_value, ltest) + weak(l, Gtest))


class SurfactantTransportEquationMassTransfer(InterfaceEquations):
    def __init__(self, surfactant_diffusivity=1):
        super(SurfactantTransportEquationMassTransfer , self).__init__()
        self.D=surfactant_diffusivity # diffusivity

    def define_fields(self):
        self.define_scalar_field("Gamma","C2", testscale=scale_factor('temporal')/scale_factor('Gamma'))
        self.define_vector_field("u_p","C2", scale=scale_factor("velocity"), testscale=1/scale_factor("velocity"))

    def define_residuals(self):
        u = var("velocity")  # velocity
        u_p, u_p_test = var_and_test("u_p")
        R = var("mesh")  # mesh
        n = var("normal")
        u_p_value = u - dot(u, n) * n + dot(partial_t(R), n) * n
        G,Gtest=var_and_test("Gamma")
        self.add_residual(weak(u_p-u_p_value, u_p_test))
        self.add_residual(weak(partial_t(G, ALE="auto")+div(u_p*G), Gtest))
        self.add_residual(weak(self.D*grad(G), grad(Gtest)))