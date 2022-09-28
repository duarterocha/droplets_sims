from scripts import *

if __name__ == "__main__":
    with DropletTempProblem() as problem:

        problem.set_c_compiler("tcc")

        problem.contact_angle = 120 * degree

        problem.Ma.value = 3 * 10 ** 2
        problem.Ra.value = 5 * 10 ** 3

        problem.solve(globally_convergent_newton=True,max_newton_iterations=100)  # solve and output
        problem.output_at_increased_time()