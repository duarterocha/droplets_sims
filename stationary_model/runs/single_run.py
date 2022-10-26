from stationary_model.problem_def import *

if __name__ == "__main__":
    with DropletTempProblem() as problem:

        problem.set_c_compiler("tcc")

        problem.contact_angle = 150 * degree
        problem.surfactants_bool = True

        problem.set_Strength(0)
        problem.set_Ma(0)
        problem.set_Ra(1)

        problem.solve(globally_convergent_newton=True,max_newton_iterations=100)  # solve and output
        problem.output_at_increased_time()