from pyoomph import *
from pyoomph.output.plotting import *
from pyoomph.expressions.units import *
from pyoomph.expressions import *


class Plotter(MatplotlibPlotter):

    def __init__(self, problem, filetrunk="plot_{:05d}", fileext="png", eigenvector=None, eigenmode="abs"):
        super(Plotter, self).__init__(problem, filetrunk, fileext, eigenvector, eigenmode)
        self.xrange = 2 * problem.droplet_radius

    def define_plot(self):
        p = self.get_problem()
        self.background_color = "darkgrey"

        # view
        ymin = -0.2 * self.xrange
        ymax = 1.2 * self.xrange
        self.set_view(-self.xrange, ymin, self.xrange, ymax)

        # colorbars
        cb_v = self.add_colorbar("velocity [$\mu$m/s]", cmap="seismic", position="bottom right", factor=1000000)
        cb_vap = self.add_colorbar("vapor conc. [kg/m$^3$]", cmap="Blues", position="top right", factor=1)
        cb_T = self.add_colorbar("temperature [$^o$C]", cmap="coolwarm", position="bottom left", offset=-273.15)
        cb_T.useOffset = False
        cb_G = self.add_colorbar("surfactants", cmap="autumn", position="top left", factor=1)
        cb_streams = self.add_colorbar("stream_sign", cmap="Greys_r", position="center")  # black to white
        cb_streams.Ndisc = 2  # but only two discrete colors (black and white)
        cb_streams._vmin = -10000  # fixed symmetric extrema: <0 will be black, >0 will be white
        cb_streams._vmax = 10000
        cb_streams.invisible = True  # Do not show the colorbar

        # plots
        self.add_plot("droplet/velocity", colorbar=cb_v)
        self.add_plot("droplet/temperature", colorbar=cb_T, transform="mirror_x")
        self.add_plot("droplet/interface/Gamma", colorbar=cb_G, transform=[None, "mirror_x"])
        self.add_plot("gas/c", colorbar=cb_vap, transform=[None, "mirror_x"])
        self.add_plot("droplet/interface", linecolor="yellow", transform=[None, "mirror_x"])

        # streamlines
        streams = self.add_plot("droplet/velocity", mode="streamlines", linewidths=2, colorbar=cb_streams,
                                transform=[None, "mirror_x"])
        '''for stream in streams:  # Iterate over non-mirrored and mirrored streamlines
            stream.colorfield = "streamfunc"  # and set the colorfield to the stream function'''

        # evaporation arrows
        arrkey = self.add_arrow_key(position="top right", title="evap. rate [g/m$^2$s]", factor=1000)
        arrkey.ymargin += 0.15  # Shift it a bit down to not overlay with the colorbar
        plt = self.add_plot("droplet/interface/evap_rate", arrowkey=arrkey, transform=[None, "mirror_x"])

        # Info text
        droplet_volume = p.get_mesh("droplet").evaluate_observable("volume")
        #Ra = round(float(droplet_volume * p.gravity * p.droplet_density_ref / (p.thermal_diffusivity * p.droplet_viscosity)), 2)
        #Ma = round(float(droplet_volume ** rational_num(1,3) * p.surface_tension_ref / (p.thermal_diffusivity * p.droplet_viscosity)),2)
        #txt = "Ma={marangoni} \nRa={rayleigh}".format(marangoni=Ma, rayleigh=Ra)
        #parameters_txt = self.add_text(txt, position="top left", textsize=20,
        #                bbox=dict(boxstyle='round', facecolor='wheat', alpha=1))
        #parameters_txt.ymargin += 0.1
        time_txt = self.add_text("t={time} s".format(time=round(float(p.get_current_time() / second),2)),
                                 position="top center", textsize=20,
                                 bbox=dict(boxstyle='round', facecolor='white', alpha=1))
        time_txt.ymargin += 0.1
