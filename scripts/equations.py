from pyoomph.generic import InterfaceEquations
from pyoomph.expressions import *


class DiffusionInfinity(InterfaceEquations):

    def __init__(self, **kwargs):
        super(DiffusionInfinity, self).__init__()
        self.inftyvals = {**kwargs}
        self.origin = vector([0])

    def define_residuals(self):
        n = self.get_normal()
        d = var("coordinate") - self.origin
        for fn, val in self.inftyvals.items():
            y, y_test = var_and_test(fn)
            self.add_residual(weak((y - val) * dot(n, d) / dot(d, d), y_test))