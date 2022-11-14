import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

main_path = "parallel_runs"
sub_dirs = next(os.walk('parallel_runs'))[1]
Ma_G_dirs = [i for i in sub_dirs if i.startswith('Ma_G_')]
diffusivity_dirs = [i for i in sub_dirs if i.startswith('diffusivity_')]

Ma_G_values = []
for file in Ma_G_dirs:
    path = os.path.join(os.path.join(main_path, file), 'observables.txt')
    df = pd.read_csv(path, sep="\t")
    Ma_G_value = float(file[file.rfind("_")+1:])
    avg_velo = df['avg_velo[m/s]'].values[-1]
    Ma_G_values.append([Ma_G_value, avg_velo])

diffusivity_values = []
default_diffusivity = 1e-10
for file in diffusivity_dirs:
    path = os.path.join(os.path.join(main_path, file), 'observables.txt')
    df = pd.read_csv(path, sep="\t")
    diffusivity_value = default_diffusivity + (1e-6-default_diffusivity)*float(file[file.rfind("_")+1:])/0.02
    avg_velo = df['avg_velo[m/s]'].values[-1]
    diffusivity_values.append([diffusivity_value, avg_velo])

Ma_G_values = np.array(Ma_G_values)
df_Ma_G = pd.DataFrame(data={'Ma_G': Ma_G_values[:,0], 'avg_velo': Ma_G_values[:,1]})
df_Ma_G=df_Ma_G.sort_values(by="Ma_G", ignore_index=True)
df_Ma_G=df_Ma_G.loc[(df_Ma_G['Ma_G'] != 0.002) & (df_Ma_G['Ma_G'] !=0.003)]
print(df_Ma_G)

diffusivity_values = np.array(diffusivity_values)
df_diffusivity = pd.DataFrame(data={'diffusivity': diffusivity_values[:,0], 'avg_velo': diffusivity_values[:,1]})
df_diffusivity=df_diffusivity.sort_values(by="diffusivity", ignore_index=True)


fig,ax = plt.subplots(ncols=2, dpi=100, figsize=(15,5))
ax[0].plot(df_Ma_G['Ma_G'], df_Ma_G['avg_velo'], color='red', label="Ma$_\Gamma$")
ax[1].plot(df_diffusivity['diffusivity'], df_diffusivity['avg_velo'], label="Diffusivity")
ax[0].set_xlabel('Ma$_\Gamma$')
ax[1].set_xlabel('Diffusivity [m$^2$/s]')
ax[0].set_ylabel('Average velocity [m/s]')
ax[1].set_xscale("log")
plt.legend()
plt.show()
