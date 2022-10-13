from problem_def import *
import numpy
from pyoomph.expressions import *
from os import listdir
import math

if __name__ == "__main__":
    with DropletTempProblem() as problem:
        problem.set_c_compiler("tcc")

        # output
        problem.plotter = [Plotter(problem)]  # Plot the normal solution
        problem.plotter.append(Plotter(problem, eigenvector=0, eigenmode="real",
                                       filetrunk="eigenreal_{:05d}"))  # real part of the eigenfunction
        problem.plotter.append(
            Plotter(problem, eigenvector=0, eigenmode="imag", filetrunk="eigenimag_{:05d}"))  # imag. part

        for p in problem.plotter:
            p.file_ext = ["png"]

        hysteresisfile = open(os.path.join(problem.get_output_directory(), "hysteresis.txt"), "w")

        problem.load_state("temporal_solve/_states/state_003560.dump")
        with problem.select_dofs() as dofs:
            dofs.unselect("gas") # Do not solve gas phase again

            problem.solve()
            problem.solve_eigenproblem(6) # Get last eigenvalues

            ds = -0.001
            while numpy.real(problem.get_last_eigenvalues()[0])<0:
                ds = problem.arclength_continuation("Ra", ds, max_ds=1000)
                problem.solve()
                problem.solve_eigenproblem(6)
                line = [problem.get_Ma(), problem.get_Ra(), problem.get_last_eigenvalues()[0]]
                hysteresisfile.write("\t".join(map(str, line)) + "\n")
                problem.output()


