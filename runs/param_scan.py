from problem_def import *
from utils import *

if __name__ == "__main__":
    with DropletTempProblem() as problem:

        problem.set_c_compiler("tcc")
        problem.max_refinement_level = 0

        parameter_scan(problem)