import matplotlib.pyplot as plt
import numpy


my_file= "../param_scan/observables.txt"
try:
    my_data=numpy.loadtxt(my_file)
except:
    raise RuntimeError("Apparently you haven't calculated the parameter scan yet")

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
ZData=sort_data[:,3]
# However, while the parameter simulation is running, it might not be complete
# Remove the incomplete Ra scan
if len(ZData)<NumMa*NumRa:
    NumMa=NumMa-1
    ZData=ZData[0:NumMa*NumRa]
    UniqueMa=UniqueMa[0:NumMa]

# Make a 2d array
ZData=ZData.reshape(NumMa,NumRa).transpose()


# Plot Ma,Ra,PosStreamFraction
plt.contourf(UniqueMa,UniqueRa,ZData,levels=30)
plt.loglog() # log-log plot
# Labels and colorbar
plt.xlabel("Ma")
plt.ylabel("Ra")
plt.colorbar(label="fraction of pos. stream func.")
plt.contour(UniqueMa,UniqueRa,ZData,levels=[0.01,0.99],colors="black")
plt.savefig(my_file[:-4]+'.png')
plt.show()