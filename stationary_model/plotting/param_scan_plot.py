from stationary_model.plotting import *

file = '/Users/duarterocha17/Desktop/droplet_sims/stability_analysis/param_scan/stability_analysis.txt'
title = "Competing Ma vs Ra"
parameter_contour_plot(file, 'positive_stream_area_fraction', title, eigen_evalutation=False)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)
plt.show()