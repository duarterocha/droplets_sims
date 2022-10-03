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

        problem.max_refinement_level = 0
        problem.initialise()

        Neigen = 6  # Numbers of eigenvalues to calculate
        #problem.solve_eigenproblem(Neigen)  # get an eigenvector as guess
        # Write eigenvalues to file

        bifurcationfile = open(os.path.join(problem.get_output_directory(), "bifurcation.txt"), "w")
        
        #raise RuntimeError("Fill in your state file below and remove this RuntimeError!")
        problem.load_state("hysteresis/_states/state_003608.dump")
        
        with problem.select_dofs() as dofs:
            dofs.unselect("gas")  # Gas is not solved for
            problem.solve() # Solve once more
            problem.solve_eigenproblem(Neigen) # solve eigenproblem. Since the gas is now deactivated, the eigenvector has different size!
            problem.activate_bifurcation_tracking("Ra","fold",eigenvector=problem.get_last_eigenvectors()[0])
            problem.solve()
            
            ds = -1
            while problem.get_Ma() >= 1:
                problem.output_at_increased_time()
                ds = problem.arclength_continuation("Ma", ds, max_ds=10)
                line = [problem.get_Ma(), problem.get_Ra()]
                bifurcationfile.write("\t".join(map(str, line)) + "\n")

            problem.deactivate_bifurcation_tracking()

            problem.solve()
            problem.solve_eigenproblem(Neigen)
            problem.activate_bifurcation_tracking("Ra", "hopf",
                                                  eigenvector=problem.get_last_eigenvectors()[0],
                                                  omega=numpy.imag(problem.get_last_eigenvalues()[0]))
            problem.solve()

            ds = 1
            while problem.get_Ma() <= 1000:
                problem.output_at_increased_time()
                ds = problem.arclength_continuation("Ma", ds, max_ds=50)
                line = [problem.get_Ma(), problem.get_Ra()]
                bifurcationfile.write("\t".join(map(str, line)) + "\n")


