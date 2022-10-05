import matplotlib.pyplot as plt
import numpy

stream_func_file= "../stability_analysis/observables.txt"
first_bif_file = "../bifurcation_tracking/observables.txt"
second_bif_file = "../bifurcation_tracking_hysteresis/observables.txt"

try:
    stream_func_data=numpy.loadtxt(stream_func_file)
    first_bif_curve=numpy.loadtxt(first_bif_file)
    second_bif_curve=numpy.loadtxt(second_bif_file)
except:
    raise RuntimeError("You made some mistake")


# Stream Function


# Due to the up-down scan, the data is not ordered well
sort_stream_func_data=stream_func_data[numpy.lexsort((stream_func_data[:,2], stream_func_data[:,1]))]

# Find the number of samples, i.e. the unique values
UniqueMa=numpy.unique(sort_stream_func_data[:,1])
NumMa=len(UniqueMa)
UniqueRa=numpy.unique(sort_stream_func_data[:,2])
NumRa=len(UniqueRa)

positive_stream_func = sort_stream_func_data[:, -1]

if len(positive_stream_func)<NumRa*NumMa:
    NumMa=NumMa-1
    positive_stream_func = positive_stream_func[0:NumRa*NumMa]
    UniqueMa=UniqueMa[0:NumMa]

positive_stream_func = positive_stream_func.reshape(NumMa, NumRa).transpose()

fig, ax = plt.subplots(figsize=(12,6), dpi=150)
cnt_stream_func = ax.contourf(UniqueMa, UniqueRa, positive_stream_func)

# Title and colorbar
ax.set_title("Stream Function")
cbar = fig.colorbar(cnt_stream_func, ticks=[0, 1])
cbar.ax.set_yticklabels(['0', '1'])
ax.set_ylabel("Ra")
ax.set_xlabel("Ma")
ax.loglog() # log-log plot
ax.set_box_aspect(1)


# First Bifurcation

Ma_first_curve = first_bif_curve[:,1]
Ra_first_curve = first_bif_curve[:,2]
#ax.plot(Ma_first_curve,Ra_first_curve, color="red")

# Second Bifurcation

Ma_second_curve = second_bif_curve[:,1]
Ra_second_curve = second_bif_curve[:,2]
#ax.plot(Ma_second_curve,Ra_second_curve, color="red", linestyle="--")

plt.xlim([UniqueMa[0], UniqueMa[-1]])
plt.ylim([UniqueRa[0], UniqueRa[-1]])
fig.savefig('bifurcations.png')
plt.show()