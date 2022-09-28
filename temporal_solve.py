from scripts import *
import numpy
from pyoomph.expressions import *

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

        problem.set_Ma(100)
        problem.set_Ra(1000)
        Ra_value = problem.get_Ra()
        problem.solve(globally_convergent_newton=True, max_newton_iterations=100)
        problem.solve_eigenproblem(6)
        while problem.get_last_eigenvalues()[0] < 0:
            Ra_value += 1000
            problem.set_Ra(Ra_value)
            problem.solve()
            problem.solve_eigenproblem(6)
        problem.set_Ra(problem.get_Ra()+5000)
        problem.solve()
        problem.perturb_dofs(0.001 * problem.get_last_eigenvectors()[0])
        problem.solve()
        problem.solve_eigenproblem(6)
        omega = numpy.imag(problem.get_last_eigenvalues()[0])
        problem.set_current_time(0)
        problem.run(2 * pi / omega, startstep=0.2 * pi / omega)

        problem.output_at_increased_time()
