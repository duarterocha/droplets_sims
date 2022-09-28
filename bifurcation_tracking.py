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

        # Scanned parameter ranges in terms of 10^...\
        NumPerOrder = 5  # Five sample points between each power of 10
        # Ma ranges
        Ma_MinRangePowerTen = 0
        Ma_MaxRangePowerTen = 5
        # Ra ranges
        Ra_MinRangePowerTen = 0
        Ra_MaxRangePowerTen = 5

        # Create the sample points
        # Ma points
        Ma_NumPts = NumPerOrder * (Ma_MaxRangePowerTen - Ma_MinRangePowerTen) + 1
        Ma_ParamRange = numpy.logspace(Ma_MinRangePowerTen, Ma_MaxRangePowerTen, Ma_NumPts, endpoint=True)
        # Ra points
        Ra_NumPts = NumPerOrder * (Ra_MaxRangePowerTen - Ra_MinRangePowerTen) + 1
        Ra_ParamRange = numpy.logspace(Ra_MinRangePowerTen, Ra_MaxRangePowerTen, Ra_NumPts, endpoint=True)

        # Initial config
        problem.set_Ma(Ma_ParamRange[0])
        problem.set_Ra(Ra_ParamRange[0])

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