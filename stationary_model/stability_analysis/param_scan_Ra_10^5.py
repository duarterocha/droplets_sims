from stationary_model.utils import *
from stationary_model.plotting import *

if __name__ == "__main__":
    with DropletTempProblem() as problem:

        problem.set_c_compiler("tcc")
        problem.max_refinement_level = 0

        problem.set_output_directory(os.path.join(os.getcwd(), "param_scan_Ra_10^5"))
        output_dir = problem.get_output_directory()

        param_scan_eigenvalues_fixed_Ra(problem, initial_Ra = 10**5, Neigen = 10)
