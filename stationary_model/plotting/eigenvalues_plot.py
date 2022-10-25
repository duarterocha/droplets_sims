import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy
import pandas as pd
import os

# Output file
file = "/Users/duarterocha17/Desktop/droplet_sims/stability_analysis/param_scan/stability_analysis.txt"

# Get data from file
df = pd.read_csv(file, sep='\t')

# Sort unsorted data
df = df.sort_values(by=['Ma', 'Ra'])
df = df.loc[df['Ra'] == 100000]
df = df.loc[df['Ma'] <= 1000]

# Get values for plotting
Strength = numpy.unique(df['Strength'].values)
Ma = numpy.unique(df['Ma'].values)
Ra = numpy.unique(df['Ra'].values)

# LaTEX fonts
mpl.rcParams['text.usetex'] = True

real_eigen, imag_eigen, legend_real, legend_imag = {},{},{},{}
for i in range(5):
    real_eigen[i] = df['real_eigen['+str(i)+']'].values
    imag_eigen[i] = df['imag_eigen['+str(i)+']'].values
    legend_real[i] = "$\lambda_{}$".format(i)
    legend_imag[i] = "$\omega_{}$".format(i)


def plot(function, legend, title, ylabel):
    fig, ax = plt.subplots(figsize=(9, 8))
    for i in range(5):
        ax.set_title(title)
        ax.plot(Ma, function[i], label = legend[i])
        ax.set_xlabel("Ma", fontsize=15)
        ax.set_ylabel(ylabel, fontsize=15)
        ax.legend(loc = 'upper right')
        plt.xscale('symlog')
        plt.yscale('symlog')

plot(real_eigen, legend_real, 'Real Eigenvalues', ylabel='$\lambda$')
plot(imag_eigen, legend_imag, 'Imaginary Eigenvalues', ylabel='$\omega$')
plt.show()



