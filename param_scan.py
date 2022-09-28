from scripts import *


if __name__ == "__main__":
    with DropletTempProblem() as problem:
        problem.set_c_compiler("tcc")

        problem.contact_angle=150*degree

        # Scanned parameter ranges in terms of 10^...\
        NumPerOrder = 5  # Five sample points between each power of 10
        # Ma ranges
        Ma_MinRangePowerTen=0
        Ma_MaxRangePowerTen=5
        # Ra ranges
        Ra_MinRangePowerTen = 0
        Ra_MaxRangePowerTen = 5

        # Create the sample points
        # Ma points
        Ma_NumPts=NumPerOrder*(Ma_MaxRangePowerTen-Ma_MinRangePowerTen)+1
        Ma_ParamRange=numpy.logspace(Ma_MinRangePowerTen,Ma_MaxRangePowerTen,Ma_NumPts,endpoint=True)
        # Ra points
        Ra_NumPts = NumPerOrder * (Ra_MaxRangePowerTen - Ra_MinRangePowerTen) + 1
        Ra_ParamRange = numpy.logspace(Ra_MinRangePowerTen, Ra_MaxRangePowerTen, Ra_NumPts, endpoint=True)

        # Initial config
        problem.set_Ma(Ma_ParamRange[0])
        problem.set_Ra(Ra_ParamRange[0])

        problem.max_refinement_level = 0

        problem.solve() # Solve the full problem

        # The gas phase does not change. So we remove it from solving
        with problem.select_dofs() as dofs:
            dofs.unselect("gas")

            # We increase the Marangoni number and do one scan in Ra up, then at the next Marangoni number, we scan Ra down
            # Thereby, there are no large jumps in the parameters and the solution can be found easier
            # Furthermore, if neighboring up and down scan results are continuous, there is no hysteresis in Ra
            # Little trick to get two neighboring sample points in a single loop
            for Ma_for_Ra_up, Ma_for_Ra_down in zip(Ma_ParamRange[::2], Ma_ParamRange[1::2]):
                problem.go_to_param(Ma=Ma_for_Ra_up)
                for Ra in Ra_ParamRange:
                    problem.go_to_param(Ra=Ra)
                    problem.solve(globally_convergent_newton=True, max_newton_iterations=100)
                    problem.output()

                problem.go_to_param(Ma=Ma_for_Ra_down)
                for Ra in reversed(Ra_ParamRange):
                    problem.go_to_param(Ra=Ra)
                    problem.solve(globally_convergent_newton=True, max_newton_iterations=100)
                    problem.output()
