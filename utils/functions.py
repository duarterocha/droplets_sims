from pyoomph import *
from pyoomph.utils.paramscan import *
from problem_def import *

'''
Parameter range
'''


def parameter_range(numperoder=5, **kwargs):
    if len(kwargs.items()) == 0:
        raise (RuntimeError("Please enter at least a parameter"))

    # Store the ranges in a dictionary
    param_dict = {}

    for key, value in kwargs.items():
        # Set limits for log scale for each parameter
        key_min_range_power_ten = 0
        key_max_range_power_ten = value

        # Get number of points to consider for each parameter
        key_num_pts = numperoder * (key_max_range_power_ten - key_min_range_power_ten) + 1
        # Get parameter range
        key_param_range = numpy.logspace(key_min_range_power_ten, key_max_range_power_ten, key_num_pts, endpoint=True)
        param_dict[key] = key_param_range

    return (param_dict)


'''
Output for stability analysis
'''


def set_title_output_with_eigen(output_file):
    title = ["Strength", "Ma", "Ra", "real_eigen", "imag_eigen", "positive_stream_area_fraction"]
    output_file.write("\t".join(map(str, title)) + "\n")


def output_with_eigen(problem, output_file, Neigen=6):
    eigvals, eigvects = problem.solve_eigenproblem(Neigen, shift=0)  # solve for 6 eigenvalues with zero shift
    line = [problem.get_Strength(), problem.get_Ma(), problem.get_Ra(), eigvals[0].real, eigvals[0].imag]
    print(problem.get_mesh("droplet").evaluate_observable("positive_stream_area_fraction"))
    line += [problem.get_mesh("droplet").evaluate_observable("positive_stream_area_fraction")]  # line to write
    output_file.write("\t".join(map(str, line)) + "\n")  # write to file eigenfile.flush()
    problem.output_at_increased_time()  # and write the output


'''
Parameter scan
'''


def parameter_scan(problem_name, calculate_eigenvalues=False):
    # Get parameter range for all parameters
    param_range = parameter_range(Ma=5, Ra=5)
    Ma_param_range = param_range["Ma"]
    Ra_param_range = param_range["Ra"]

    # Initial config
    problem_name.set_Strength(0)
    problem_name.set_Ma(Ma_param_range[0])
    problem_name.set_Ra(Ra_param_range[0])

    problem_name.solve()  # Solve the full problem

    # If also performing stability analysis, create a file to write output
    if calculate_eigenvalues:
        eigenfile = open(os.path.join(problem_name.get_output_directory(), "stability_analysis.txt"), "w")
        set_title_output_with_eigen(problem_name, eigenfile)

    # The gas phase does not change. So we remove it from solving
    with problem_name.select_dofs() as dofs:
        dofs.unselect("gas")

        # We increase the Marangoni number and do one scan in Ra up, then at the next Marangoni number, we scan Ra down
        for Ma_for_Ra_up, Ma_for_Ra_down in zip(Ma_param_range[::2], Ma_param_range[1::2]):
            problem_name.go_to_param(Ma=Ma_for_Ra_up)
            for Ra in Ra_param_range:
                problem_name.go_to_param(Ra=Ra)
                problem_name.solve(globally_convergent_newton=True, max_newton_iterations=100)
                if calculate_eigenvalues:
                    output_with_eigen(problem_name, eigenfile)
                else:
                    problem_name.output()

            problem_name.go_to_param(Ma=Ma_for_Ra_down)
            for Ra in reversed(Ra_param_range):
                problem_name.go_to_param(Ra=Ra)
                problem_name.solve(globally_convergent_newton=True, max_newton_iterations=100)
                if calculate_eigenvalues:
                    output_with_eigen(problem_name, eigenfile)
                else:
                    problem_name.output()


def plot_eigen_values(problem_name, plotter_object):
    # Plot the normal solution
    problem_name.plotter = [plotter_object(problem_name)]
    # Real part of the eigenfunction
    problem_name.plotter.append(
        plotter_object(problem_name, eigenvector=0, eigenmode="real", filetrunk="eigenreal_{:05d}"))
    # Imaginary part of the eigenfunction
    problem_name.plotter.append(
        plotter_object(problem_name, eigenvector=0, eigenmode="imag", filetrunk="eigenimag_{:05d}"))
    # Magnitude of the eigenfunction
    problem_name.plotter.append(
        plotter_object(problem_name, eigenvector=0, eigenmode="abs", filetrunk="eigenabs_{:05d}"))

    for p in problem_name.plotter:
        p.file_ext = ["png"]  # plot both png and pdf for all plotters


'''
Bifurcation tracking
'''


def set_title_output_with_bifurcation(output_file):
    title = ["Ma", "Ra"]
    output_file.write("\t".join(map(str, title)) + "\n")


def jump_on_bifurcation(problem_name, bifurcation_type, eigenvector, omega, bifurcationfile, ds=0.001, lower_Ma=1,
                        upper_Ma=10000):
    # Activate bifurcation tracking and solve
    problem_name.activate_bifurcation_tracking("Ra", bifurcation_type, eigenvector=eigenvector, omega=omega)
    problem_name.solve()

    # Arc length continuation do jump on bifurcations, going up
    if ds > 0:
        while problem_name.get_Ma() < upper_Ma:
            problem_name.output_at_increased_time()
            ds = problem_name.arclength_continuation("Ma", ds)
            line = [problem_name.get_Ma(), problem_name.get_Ra()]
            bifurcationfile.write("\t".join(map(str, line)) + "\n")
    # Arc length continuation do jump on bifurcations, going down
    else:
        while problem_name.get_Ma() > lower_Ma:
            problem_name.output_at_increased_time()
            ds = problem_name.arclength_continuation("Ma", ds)
            line = [problem_name.get_Ma(), problem_name.get_Ra()]
            bifurcationfile.write("\t".join(map(str, line)) + "\n")

    problem_name.deactivate_bifurcation_tracking()


def bifurcation_tracking(problem_name, init_guess_param_eigenvector=[10, 100], starting_Ma = 100, bifurcation_type="hopf"):
    # Note: If the starting Ma number is around 10, it might not find the bifurcation

    # Initial config
    problem_name.set_Ma(init_guess_param_eigenvector[0])
    problem_name.set_Ra(init_guess_param_eigenvector[1])

    # First solve
    problem_name.solve()
    # Get an eigenvector as guess
    problem_name.solve_eigenproblem(6)
    eigenvector = problem_name.get_last_eigenvectors()[0]
    omega = None
    if bifurcation_type == 'hopf':
        omega = numpy.imag(problem_name.get_last_eigenvalues()[0])

    # Write eigenvalues to file
    bifurcationfile = open(os.path.join(problem_name.get_output_directory(), "bifurcation.txt"), "w")
    set_title_output_with_bifurcation(bifurcationfile)

    with problem_name.select_dofs() as dofs:
        dofs.unselect("gas")  # Gas is not solved for

        # Find bifurcation for the starting Ma
        problem_name.set_Ma(starting_Ma)
        for parameter in problem_name.find_bifurcation_via_eigenvalues("Ra", 1000, epsilon=1e-4):
            print(parameter)

        # Jump on bifurcations going up
        ds = 0.001
        jump_on_bifurcation(problem_name, bifurcation_type, eigenvector, omega, bifurcationfile, ds)

        # Solve for last parameter set
        problem_name.solve()
        problem_name.solve_eigenproblem(6)

        ds = -10
        jump_on_bifurcation(problem_name, bifurcation_type, eigenvector, omega, bifurcationfile, ds)
