from scripts import *
import matplotlib.pyplot as plt
import numpy

my_file= "../temporal_solve/observables.txt"
try:
    my_data=numpy.loadtxt(my_file)
except:
    raise RuntimeError("Apparently you haven't calculated the temporal analysis yet")

time = my_data[:,0]
avg_velo = my_data[:,3]

fig,ax = plt.subplots()
ax.plot(time, avg_velo)
ax.set_xlabel("Time (s)")
ax.set_ylabel("Average velocity")
plt.show()


