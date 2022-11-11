from pyoomph import *
from pyoomph.equations.advection_diffusion import *
from pyoomph.expressions.units import *


'''
Define geometry
'''


class Mesh(GmshTemplate):

    def __init__(self, size=1, radius=1, contact_angle=90 * degree):
        super(Mesh, self).__init__()
        self.size = size
        self.radius = radius
        self.contact_angle = contact_angle

    def define_geometry(self):
        S = self.size

        # points
        p_origin = (0, 0)
        p_center = (0, - self.radius * cos(self.contact_angle) * S)
        p_surface = self.point(self.radius * sin(self.contact_angle) * S, 0)
        p_apex = self.point(0, self.radius * (1 - cos(self.contact_angle)) * S)
        gas_radius = self.radius * S * 3.5
        p_surface_infinity = self.point(gas_radius, 0, size = 1.25)
        p_apex_infinity = self.point(0, gas_radius, size = 1.25)

        # lines
        self.create_lines(p_surface_infinity, "gas_surface", p_surface, "droplet_surface", p_origin, "droplet_axis",
                          p_apex, "gas_axis", p_apex_infinity)
        self.circle_arc(p_surface, p_apex, center=p_center, name="interface")
        self.circle_arc(p_surface_infinity, p_apex_infinity, center=p_origin, name="gas_infinity")

        # droplet
        self.plane_surface("droplet_surface", "droplet_axis", "interface", name="droplet")

        # gas phase
        self.plane_surface("gas_surface", "gas_infinity", "gas_axis", "interface", name="gas")
