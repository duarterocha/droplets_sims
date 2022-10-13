import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy
import pandas as pd
import os


def get_data_for_contour(file, *args):

    # Get data from file
    df = pd.read_csv(file, sep='\t')

    # Sort unsorted data
    df = df.sort_values(by=['Ma', 'Ra'])

    # Get values for plotting
    Strength = numpy.unique(df['Strength'].values)
    Ma = numpy.unique(df['Ma'].values)
    NumMa = len(Ma)
    Ra = numpy.unique(df['Ra'].values)
    NumRa = len(Ra)

    # In case a contour is desired, it should be specified the desired name for the data (key)
    # and the column name on the file
    if len(args) != 0:
        for v in args:
            zdata = df[str(v)].values

    # Remove data from unfinished calculation
    if len(zdata) < NumMa * NumRa:
        NumMa = NumMa - 1
        zdata = zdata[0:NumRa * NumMa]
        Ma = Ma[0:NumMa]

    # Make a 2d array
    zdata = zdata.reshape(NumMa, NumRa).transpose()

    return (Strength, Ma, Ra, zdata)

def personalize_cbars(ax, cnt, cnt_lines, strength_value, title, eigen_evalutation = False):

    # Set title
    ax.set_title(title)

    # Axis
    ax.loglog()
    ax.set_aspect('equal')
    ax.set_xlabel("Ma", fontsize=15)
    ax.set_ylabel("Ra", fontsize=15)

    # Colorbar, legend
    if eigen_evalutation:
        cnt_ticks, cnt_lines_ticks = [-1, 1], [0]
        cnt_yticklabels, cnt_lines_yticklabels = ['Negative', 'Positive'], ['0']
    else:
        cnt_ticks, cnt_lines_ticks = [0, 1], [0.5]
        cnt_yticklabels, cnt_lines_yticklabels = ['Ma$_{Dom}$', 'Ra$_{Dom}$'], ['Ma vs Ra']

    cbar_lines = plt.colorbar(cnt_lines, ticks=cnt_lines_ticks, shrink=0.05, aspect=1, pad=0)
    cbar_lines.ax.set_yticklabels(cnt_lines_yticklabels)
    cbar_lines.ax.tick_params(size=0)
    cbar = plt.colorbar(cnt, ticks=cnt_ticks, shrink=0.7)
    cbar.ax.set_yticklabels(cnt_yticklabels)

def parameter_contour_plot(file, plot_function = "positive_stream_area_fraction", title = "Competing Ma vs Ra", eigen_evalutation = False):

    # LaTEX fonts
    mpl.rcParams['text.usetex'] = True

    Strength, Ma, Ra, plot_func = get_data_for_contour(file, plot_function)

    if eigen_evalutation:
        plot_func = numpy.where(plot_func > 0, 1, plot_func)
        plot_func = numpy.where(plot_func < 0, -1, plot_func)

    # Contour from data
    fig, ax = plt.subplots()
    cnt = ax.contourf(Ma, Ra, plot_func, cmap = "summer")

    levels = [0.01, 0.99] if not eigen_evalutation else [-0.99, 0.99]
    cnt_lines = ax.contourf(Ma, Ra, plot_func, levels=levels, colors="olive", linestyles='dashed')

    # Axis, labels, colorbar
    personalize_cbars(ax, cnt, cnt_lines, Strength, title, eigen_evalutation)

def bifurcation_plot(files_bifurcation, file_param_scan):

    parameter_contour_plot(file_param_scan)

    def plotting_method(file):
        # Get data from file
        df = pd.read_csv(file, sep='\t')

        # Sort unsorted data
        df = df.sort_values(by=['Ma', 'Ra'])

        # Get Ma and Ra
        Ma = df['Ma'].values
        Ra = df['Ra'].values

        plt.plot(Ma, Ra, linestyle='dashed')

    if type(files_bifurcation) != 'list':
        plotting_method(files_bifurcation)
    else:
        for file in files_bifurcation:
            plotting_method(file)








