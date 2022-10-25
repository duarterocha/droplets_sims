from temporal_model.problem_def import *

if __name__ == "__main__":
    with DropletTempProblem() as problem:

        problem.set_c_compiler("tcc")

        problem.contact_angle = 150 * degree
        problem.surfactants_bool = True

        problem.set_Strength(1)
        problem.set_Ma(10)
        problem.set_Ra(1)

        problem.run(10,startstep=0.01,outstep=True,temporal_error=1)
        problem.output_at_increased_time()