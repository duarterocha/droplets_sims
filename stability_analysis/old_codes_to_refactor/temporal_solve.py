from problem_def import *
import numpy
from pyoomph.expressions import *
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

        problem.default_timestepping_scheme = "Newmark2"

        problem.set_Ma(100)
        problem.set_Ra(1000)
        Ra_value = problem.get_Ra()
        problem.solve(globally_convergent_newton=True, max_newton_iterations=100)
        problem.solve_eigenproblem(6)
        eigenfile = open(os.path.join(problem.get_output_directory(), "eigenvalues.txt"), "w")
        while numpy.real(problem.get_last_eigenvalues()[0]) < 0:
            Ra_value += 1000
            problem.set_Ra(Ra_value)
            problem.solve()
            problem.solve_eigenproblem(6)
            line = [problem.get_Ma(), problem.get_Ra(), numpy.real(problem.get_last_eigenvalues()[0]),
                    numpy.imag(problem.get_last_eigenvalues()[0])]  # line to write
            eigenfile.write("\t".join(map(str, line)) + "\n")  # write to file eigenfile.flush()
        problem.output()
        problem.perturb_dofs(0.01*numpy.real(problem.get_last_eigenvectors()[0]))
        problem.set_current_time(0)
        omega = numpy.imag(problem.get_last_eigenvalues()[0])
        problem.run(300 * math.pi / omega, startstep=0.05 * math.pi / omega)

        problem.output_at_increased_time()
