from problem_def import *
import numpy
from pyoomph.expressions import *
from utils import *
import math

if __name__ == "__main__":
    with DropletTempProblem() as problem:
        problem.set_c_compiler("tcc")

        plot_eigen_values(problem, Plotter)
        problem.default_timestepping_scheme = "Newmark2"

        find_unstable_from_init_stable(problem)

        # Perturb solution
        problem.perturb_dofs(0.01*numpy.real(problem.get_last_eigenvectors()[0]))

        # Set time to 0
        problem.set_current_time(0)

        # Run temporal solve
        omega = numpy.imag(problem.get_last_eigenvalues()[0])
        problem.run(100 * math.pi / omega, startstep=0.02 * math.pi / omega)

        problem.output_at_increased_time()
