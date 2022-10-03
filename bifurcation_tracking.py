import pyoomph.solvers.petsc
from scripts import *

if __name__ == "__main__":
    with DropletTempProblem() as problem:

        problem.plotter = [Plotter(problem)]  # Plot the normal solution
        problem.plotter.append(Plotter(problem, eigenvector=0, eigenmode="real",
                                       filetrunk="eigenreal_{:05d}"))  # real part of the eigenfunction
        problem.plotter.append(
            Plotter(problem, eigenvector=0, eigenmode="imag", filetrunk="eigenimag_{:05d}"))  # imag. part
        problem.plotter.append(
            Plotter(problem, eigenvector=0, eigenmode="abs", filetrunk="eigenabs_{:05d}"))  # magnitude

        for p in problem.plotter:
            p.file_ext = ["png"]  # plot both png and pdf for all plotters

        # Initial config
        problem.set_Ma(10)
        problem.set_Ra(100)

        problem.contact_angle = 120 * degree
        problem.max_refinement_level = 0
        problem.solve()

        Neigen = 6  # Numbers of eigenvalues to calculate
        problem.solve_eigenproblem(Neigen)  # get an eigenvector as guess
        # Write eigenvalues to file

        bifurcationfile = open(os.path.join(problem.get_output_directory(), "bifurcation.txt"), "w")

        with problem.select_dofs() as dofs:
            dofs.unselect("gas")  # Gas is not solved for

            problem.set_Ma(100)

            for parameter in problem.find_bifurcation_via_eigenvalues("Ra", 1000, epsilon=1e-4):
                print(parameter)
            problem.activate_bifurcation_tracking("Ra", "hopf",
                        eigenvector=problem.get_last_eigenvectors()[0], omega=numpy.imag(problem.get_last_eigenvalues()[0]))
            problem.solve()
            ds = 0.001
            while problem.get_Ma() < 1000:
                problem.output_at_increased_time()
                ds = problem.arclength_continuation("Ma", ds)
                line = [problem.get_Ma(), problem.get_Ra()]
                bifurcationfile.write("\t".join(map(str, line)) + "\n")

            problem.deactivate_bifurcation_tracking()

            problem.solve()
            problem.solve_eigenproblem(Neigen)
            problem.activate_bifurcation_tracking("Ra", "hopf",
                                                  eigenvector=problem.get_last_eigenvectors()[0],
                                                  omega=numpy.imag(problem.get_last_eigenvalues()[0]))
            problem.solve()

            ds = -10
            while problem.get_Ma() > 1:
                problem.output_at_increased_time()
                ds = problem.arclength_continuation("Ma", ds, max_ds=50)
                line = [problem.get_Ma(), problem.get_Ra()]
                bifurcationfile.write("\t".join(map(str, line)) + "\n")

