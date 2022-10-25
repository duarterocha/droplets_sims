from stationary_model.utils import *
from stationary_model.plotting import *

if __name__ == "__main__":
    with DropletTempProblem() as problem:

        problem.set_c_compiler("tcc")
        problem.max_refinement_level = 0

        problem.set_output_directory(os.path.join(os.getcwd(), "../../stability_analysis/stability_analysis"))
        output_dir = problem.get_output_directory()

        parameter_scan(problem, calculate_eigenvalues=True, Neigen=10, step_per_order=10)
