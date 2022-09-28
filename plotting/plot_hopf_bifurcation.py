from scripts import *
import matplotlib.pyplot as plt
import numpy

my_file= "../bifurcation_tracking/bifurcation.txt"
try:
    my_data=numpy.loadtxt(my_file)
except:
    raise RuntimeError("Apparently you haven't calculated the bifurcation line yet")

if __name__ == "__main__":
    with DropletTempProblem() as problem:
        problem.plotter = [Plotter(problem)]  # Plot the normal solution
        problem.plotter.append(Plotter(problem, eigenvector=0, eigenmode="real",
                                          filetrunk="eigenreal_{:05d}"))  # real part of the eigenfunction
        problem.plotter.append(
            Plotter(problem, eigenvector=0, eigenmode="imag", filetrunk="eigenimag_{:05d}"))  # imag. part
        problem.plotter.append(
            Plotter(problem, eigenvector=0, eigenmode="abs", filetrunk="eigenabs_{:05d}"))  # magnitude

        for p in problem.plotter:
            p.file_ext = ["png", "pdf"]  # plot both png and pdf for all plotters

