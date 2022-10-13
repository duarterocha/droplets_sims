import matplotlib.pyplot as plt
from plotting import *

file = '/Users/duarterocha17/Desktop/droplet_sims/stability_analysis/param_scan/stability_analysis.txt'

title = "Competing Marangoni vs Rayleigh"
parameter_contour_plot(file, 'positive_stream_area_fraction', title)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)

title = "Real eigenvalues"
parameter_contour_plot(file, 'real_eigen', title, eigen_evalutation=True)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)

title = "Imaginary eigenvalues"
parameter_contour_plot(file, 'imag_eigen', title, eigen_evalutation=True)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)