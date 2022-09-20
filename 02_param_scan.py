from scripts import *

if __name__ == "__main__":
    with DropletTempProblem() as problem:
        problem.set_c_compiler("tcc")

        problem.contact_angle=120*degree

        # Scanned parameter ranges in terms of 10^...
        MinRangePowerTen=0
        MaxRangePowerTen=5
        NumPerOrder=5 # Five sample points between each power of 10

        # Create the sample points
        NumPts=NumPerOrder*(MaxRangePowerTen-MinRangePowerTen)+1
        ParamRange=numpy.logspace(MinRangePowerTen,MaxRangePowerTen,NumPts,endpoint=True)

        # Initial config
        problem.Ma.value = ParamRange[0]
        problem.Ra.value = ParamRange[0]

        problem.solve() # Solve the full problem

        # The gas phase does not change. So we remove it from solving
        with problem.select_dofs() as dofs:
            dofs.unselect("gas")

            # We increase the Marangoni number and do one scan in Ra up, then at the next Marangoni number, we scan Ra down
            # Thereby, there are no large jumps in the parameters and the solution can be found easier
            # Furthermore, if neighboring up and down scan results are continuous, there is no hysteresis in Ra
            # Little trick to get two neighboring sample points in a single loop
            for Ma_for_Ra_up, Ma_for_Ra_down in zip(ParamRange[::2], ParamRange[1::2]):
                problem.Ma.value = Ma_for_Ra_up
                for Ra in ParamRange:
                    problem.Ra.value = Ra
                    problem.solve(globally_convergent_newton=True, max_newton_iterations=100)
                    problem.output()

                problem.Ma.value = Ma_for_Ra_down
                for Ra in reversed(ParamRange):
                    problem.Ra.value = Ra
                    problem.solve(globally_convergent_newton=True, max_newton_iterations=100)
                    problem.output()