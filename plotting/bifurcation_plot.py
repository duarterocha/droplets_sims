import matplotlib.pyplot as plt
from plotting import *

file = '/Users/duarterocha17/Desktop/droplet_sims/stability_analysis/param_scan/stability_analysis.txt'
title = "Competing Marangoni vs Rayleigh"
parameter_contour_plot(file, 'real_eigen', title, eigen_evalutation=True)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)
plt.show()