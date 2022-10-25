from stationary_model.plotting import *

'''file_param_scan = '/Users/duarterocha17/Desktop/droplet_sims/stability_analysis/param_scan/observables.txt'
title = "Average Velocity Magnitude"
parameter_avg_velo_contour_plot(file_param_scan, 'avg_velo[]', title)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)
plt.show()'''

title = "Average Velocity Magnitude"
avg_velo_plot('/Users/duarterocha17/Desktop/droplet_sims/stability_analysis/temporal_solve/observables.txt')
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)
plt.show()