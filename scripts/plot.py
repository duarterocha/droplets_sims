from pyoomph import *
from pyoomph.output.plotting import *


class Plotter(MatplotlibPlotter):

	def __init__(self, problem, filetrunk="plot_{:05d}", fileext="png", eigenvector=None, eigenmode="abs"):
		super(Plotter, self).__init__(problem, filetrunk, fileext, eigenvector, eigenmode)
		self.xrange = 1.25

	def define_plot(self):

		p = self.get_problem()
		self.background_color = "darkgrey"

		# view
		ymin = -0.2 * self.xrange
		ymax = 1.2 * self.xrange
		self.set_view(-self.xrange, ymin, self.xrange, ymax)

		# colorbars
		cb_v = self.add_colorbar("velocity", cmap="seismic", position="bottom right", factor=10)
		cb_vap = self.add_colorbar("vapor", cmap="Blues", position="top right", factor=1)
		cb_T = self.add_colorbar("temperature", cmap="coolwarm", position="bottom left")
		cb_streams = self.add_colorbar("stream_sign", cmap="Greys_r", position="center")  # black to white
		cb_streams.Ndisc = 2  # but only two discrete colors (black and white)
		cb_streams._vmin = -10000  # fixed symmetric extrema: <0 will be black, >0 will be white
		cb_streams._vmax = 10000
		cb_streams.invisible = True  # Do not show the colorbar

		# plots
		self.add_plot("droplet/velocity", colorbar=cb_v)
		self.add_plot("droplet/T", colorbar=cb_T, transform="mirror_x")
		self.add_plot("gas/c", colorbar=cb_vap, transform=[None, "mirror_x"])
		self.add_plot("droplet/interface", linecolor="yellow", transform=[None, "mirror_x"])

		# stream lines
		streams = self.add_plot("droplet/velocity", mode="streamlines", linewidths=2, colorbar=cb_streams,
								transform=[None, "mirror_x"])
		for stream in streams:  # Iterate over non-mirrored and mirrored stream lines
			stream.colorfield = "streamfunc"  # and set the colorfield to the stream function

		# evaporation arrows
		arrkey = self.add_arrow_key(position="top right", title="evaporation")
		arrkey.ymargin += 0.15  # Shift it a bit down to not overlay with the colorbar
		plt = self.add_plot("droplet/interface/evap_rate", arrowkey=arrkey, transform=[None, "mirror_x"])


		# Info text
		txt = "Ma={marangoni}, Ra={rayleigh}".format(marangoni=str(round(p.Ma.value,2)), rayleigh=str(round(p.Ra.value),2))
		self.add_text(txt, position="top left", textsize=20,
					  bbox=dict(boxstyle='round', facecolor='wheat', alpha=1))