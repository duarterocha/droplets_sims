import pyoomph.solvers.petsc
from problem_def import *

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

        

        bifurcationfile=None
        
        
        def out_bifurc_file():
            line = [problem.get_Ma(), problem.get_Ra()]
            evs=problem.get_last_eigenvalues()
            for ev in evs:
            	line+=[ numpy.real(ev),numpy.imag(ev)]
            bifurcationfile.write("\t".join(map(str, line)) + "\n")
            bifurcationfile.flush()                
            problem.output_at_increased_time()
        
        def do_scan(ds,MaBound,BoundSign):
            problem.load_state("hysteresis/_states/state_003608.dump")

            with problem.select_dofs() as dofs:
                dofs.unselect("gas")  # Gas is not solved for
                problem.solve() # Solve once more
                problem.solve_eigenproblem(Neigen) # solve eigenproblem. Since the gas is now deactivated, the eigenvector has different size!            
                bifurc_ev=problem.get_last_eigenvectors()[0]
                problem.activate_bifurcation_tracking("Ra","fold",eigenvector=bifurc_ev)
                while problem.get_Ma()*BoundSign >= MaBound*BoundSign:                                   
                   problem.deactivate_bifurcation_tracking() # Unfortunately, we must deactivate bifurcation tracking to see more eigenvalues
                   problem.solve_eigenproblem(Neigen) # solve eigenproblem. Since the gas is now deactivated, the eigenvector has different size!            
                   bifurc_ev=problem.get_last_eigenvectors()[0]
                   out_bifurc_file()
		             
		             # Go one step with bifurcation tracking
                   problem.activate_bifurcation_tracking("Ra","fold",eigenvector=bifurc_ev)
                   problem.solve()
                   ds = problem.arclength_continuation("Ma", ds)
                   problem.reset_arc_length_parameters()                
                   problem.deactivate_bifurcation_tracking()

        bifurcationfile = open(os.path.join(problem.get_output_directory(), "bifurcation_down.txt"), "w")            
        do_scan(-1,1,1)    

        

