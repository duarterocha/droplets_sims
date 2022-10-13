import pyoomph.solvers.petsc
from runs.param_scan import *
from problem_def import *
from utils import *
from plotting import *

if __name__ == "__main__":
    with DropletTempProblem() as problem:

        problem.set_c_compiler("tcc")
        problem.max_refinement_level = 0

        problem.set_output_directory(os.path.join(os.getcwd(), "param_scan"))
        output_dir = problem.get_output_directory()

        parameter_scan(problem, calculate_eigenvalues=True)

        '''parameter_contour_plot(os.path.join(output_dir, "stability_analysis.txt"), 'positive_stream_area_fraction')
        parameter_contour_plot(os.path.join(output_dir, "stability_analysis.txt"), 'real_eigen', title = 'Real Eigenvalues', eigen_evalutation=True)
        parameter_contour_plot(os.path.join(output_dir, "stability_analysis.txt"), 'imag_eigen',
                               title='Imaginary Eigenvalues', eigen_evalutation=True)'''
