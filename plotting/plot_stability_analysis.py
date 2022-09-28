import matplotlib.pyplot as plt
import numpy

my_file= "../stability_analysis/eigenvalues.txt"
try:
    my_data=numpy.loadtxt(my_file)
except:
    raise RuntimeError("Apparently you haven't calculated the eigenvalues yet")

# Due to the up-down scan, the data is not ordered well
sort_data=my_data[numpy.lexsort((my_data[:,1], my_data[:,0]))]
# Now it is
#       print(sort_data); exit()

# Find the number of samples, i.e. the unique values
UniqueMa=numpy.unique(sort_data[:,0])
NumMa=len(UniqueMa)
UniqueRa=numpy.unique(sort_data[:,1])
NumRa=len(UniqueRa)

# Plot the positive stream covered area
positive_stream_func=sort_data[:,2]
eigenvalues_real=sort_data[:,3]
eigenvalues_imag=sort_data[:,4]
# However, while the parameter simulation is running, it might not be complete
# Remove the incomplete Ra scan
if len(eigenvalues_real)<NumRa*NumMa:
    NumMa=NumMa-1
    positive_stream_func = positive_stream_func[0:NumRa*NumMa]
    eigenvalues_real=eigenvalues_real[0:NumRa*NumMa]
    eigenvalues_imag=eigenvalues_imag[0:NumRa*NumMa]
    UniqueMa=UniqueMa[0:NumMa]

for value in range(len(eigenvalues_real)):
    # get positive real values
    if eigenvalues_real[value] > 0:
        eigenvalues_real[value] = 1
    elif eigenvalues_real[value] < 0:
        eigenvalues_real[value] = -1
    if eigenvalues_imag[value] > 0:
        eigenvalues_imag[value] = 1

# Make a 2d array
positive_stream_func=positive_stream_func.reshape(NumMa,NumRa).transpose()
eigenvalues_real=eigenvalues_real.reshape(NumMa,NumRa).transpose()
eigenvalues_imag=eigenvalues_imag.reshape(NumMa,NumRa).transpose()

# Plot Ma,Ra,PosStreamFraction
fig, ax = plt.subplots(1, 3, figsize=(12,6), sharey=True)
cnt_stream_func = ax[0].contourf(UniqueMa, UniqueRa,positive_stream_func)
cnt_real = ax[1].contourf(UniqueMa,UniqueRa,eigenvalues_real)
cnt_imag = ax[2].contourf(UniqueMa,UniqueRa,eigenvalues_imag)
# Title and colorbar
ax[0].set_title("Positive Stream Function")
ax[1].set_title("Real Eigenvalues")
ax[2].set_title("Imaginary Eigenvalues")
plt.colorbar(cnt_stream_func, ax=ax[0], orientation="horizontal", ticks=[0, 1])
plt.colorbar(cnt_real, ax=ax[1], orientation="horizontal", ticks=[-1, 0, 1])
plt.colorbar(cnt_imag, ax=ax[2], orientation="horizontal", ticks=[0, 1])
ax[0].set_ylabel("Ra")
# Labels
for x in range(3):
    ax[x].loglog() # log-log plot
    ax[x].set_xlabel("Ma")
    ax[x].set_box_aspect(1)
fig.savefig(my_file[:-4]+'.png')
plt.show()