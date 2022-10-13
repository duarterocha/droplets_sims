import pyoomph.solvers.petsc
from problem_def import *
from utils import *

if __name__ == "__main__":
    with DropletTempProblem() as problem:

        plot_eigen_values(problem, Plotter)
        problem.max_refinement_level = 0

        bifurcation_tracking(problem, init_guess_param_eigenvector=[10, 100])