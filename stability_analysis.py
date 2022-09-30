import pyoomph.solvers.petsc
from scripts import *

if __name__ == "__main__":
    with DropletTempProblem() as problem:
        # Scanned parameter ranges in terms of 10^...\
        NumPerOrder = 5  # Five sample points between each power of 10
        # Ma ranges
        Ma_MinRangePowerTen = 0
        Ma_MaxRangePowerTen = 3
        # Ra ranges
        Ra_MinRangePowerTen = 0
        Ra_MaxRangePowerTen = 5

        # Create the sample points
        # Ma points
        Ma_NumPts = NumPerOrder * (Ma_MaxRangePowerTen - Ma_MinRangePowerTen) + 1
        Ma_ParamRange = numpy.logspace(Ma_MinRangePowerTen, Ma_MaxRangePowerTen, Ma_NumPts, endpoint=True)
        # Ra points
        Ra_NumPts = NumPerOrder * (Ra_MaxRangePowerTen - Ra_MinRangePowerTen) + 1
        Ra_ParamRange = numpy.logspace(Ra_MinRangePowerTen, Ra_MaxRangePowerTen, Ra_NumPts, endpoint=True)

        # Initial config
        problem.set_Ma(Ma_ParamRange[0])
        problem.set_Ra(Ra_ParamRange[0])

        problem.contact_angle = 120 * degree
        problem.max_refinement_level = 0
        problem.solve()

        Neigen = 6  # Numbers of eigenvalues to calculate
        # Write eigenvalues to file
        eigenfile=open(os.path.join(problem.get_output_directory(),"eigenvalues.txt"),"w")

        def output_with_eigen():
            eigvals, eigvects = problem.solve_eigenproblem(Neigen, shift=0)  # solve for 6 eigenvalues with zero shift
            positive_stream_area_fraction=problem.get_mesh("droplet").evaluate_observable("positive_stream_area_fraction")
            line = [problem.get_Ma(), problem.get_Ra(), positive_stream_area_fraction, eigvals[0].real, eigvals[0].imag]  # line to write
            eigenfile.write("\t".join(map(str, line)) + "\n")  # write to file eigenfile.flush()
            problem.output_at_increased_time()  # and write the output


        with problem.select_dofs() as dofs:
            dofs.unselect("gas")  # Gas is not solved for
            for Ma_for_Ra_up, Ma_for_Ra_down in zip(Ma_ParamRange[::2], Ma_ParamRange[1::2]):
                    problem.go_to_param(Ma=Ma_for_Ra_up)
                    for Ra in Ra_ParamRange:
                        problem.solve()
                        problem.go_to_param(Ra=Ra)
                        output_with_eigen()

                    problem.go_to_param(Ma=Ma_for_Ra_down)
                    for Ra in reversed(Ra_ParamRange):
                        problem.solve()
                        problem.go_to_param(Ra=Ra)
                        output_with_eigen()

