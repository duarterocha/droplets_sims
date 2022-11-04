from temporal_model_dimensional.problem_def import *

if __name__ == "__main__":
    with EvapWaterDropletProblem() as problem:

        problem.set_c_compiler("tcc")

        problem.contact_angle = 120 * degree
        problem.max_refinement_level=0

        problem.run(10*second,startstep=0.001*second,outstep=True,temporal_error=1)
        problem.output_at_increased_time()