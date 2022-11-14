from temporal_model_dimensional.problem_def import *

if __name__ == "__main__":
    with EvapWaterDropletProblem() as problem:

        problem.set_c_compiler("tcc")

        problem.surfactants_diffusivity = 1e-9 * meter ** 2 / second
        problem.set_Ma_G(0.002)  # Surfactants reduce the surface tension of pure water by 1%

        problem.max_refinement_level = 0

        problem.run(50*second,startstep=0.001*second,outstep=True,temporal_error=1, out_initially=False, maxstep=2*second)
        problem.output_at_increased_time()