from stationary_model.plotting import *

file_param_scan = '/Users/duarterocha17/Desktop/droplet_sims/stability_analysis/param_scan/stability_analysis.txt'
file_bifurcation = '/Users/duarterocha17/Desktop/droplet_sims/stability_analysis/bifurcation_tracking/bifurcation.txt'
title = "Bifurcations"
bifurcation_plot(file_bifurcation, file_param_scan)
plt.legend(loc = 'lower right')
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)
plt.show()