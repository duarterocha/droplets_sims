from stationary_model.plotting import *

file = '/Users/duarterocha17/Desktop/droplet_sims/stability_analysis/stability_analysis/stability_analysis.txt'

title = "Competing Marangoni vs Rayleigh"
parameter_contour_plot(file, 'positive_stream_area_fraction', title)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)

title = "Real eigenvalues (0)"
parameter_contour_plot(file, 'real_eigen[0]', title, eigen_evalutation=True)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)

title = "Imaginary eigenvalues (0)"
parameter_contour_plot(file, 'imag_eigen[0]', title, eigen_evalutation=True)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)

'''title = "Real eigenvalues (1)"
parameter_contour_plot(file, 'real_eigen1', title, eigen_evalutation=True)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)

title = "Imaginary eigenvalues (1)"
parameter_contour_plot(file, 'imag_eigen1', title, eigen_evalutation=True)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)

title = "Real eigenvalues (2)"
parameter_contour_plot(file, 'real_eigen2', title, eigen_evalutation=True)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)

title = "Imaginary eigenvalues (2)"
parameter_contour_plot(file, 'imag_eigen2', title, eigen_evalutation=True)
plt.savefig(os.path.join(os.getcwd(), 'output_figures/'+title+".png"), dpi = 100)'''