import openmc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import glob
from IPython.display import Image

# from IPython.display import image

listKeff = ["098", "096", "093","088"]
fig, normAx = plt.subplots()
# relDiffAx = normAx.twinx()
normAx.set_xlabel("Radius [cm]")
normAx.set_ylabel("Normalized Flux [-]")

# TODO radius is currently bins, not cm. May take some work to change x axis size for all the lines. Cylmesh is currently set to radii[3] so 227.5cm
width = [125, 2.5, 75, 25]
radii = [width[0], width[0]+width[1], width[0]+width[1]+width[2] , width[0]+width[1]+width[2]+width[3]]
height = 300

# x_ax_geom = range(0, 100)# *radii[3]
x_ax_geom = [float(x)*radii[3]/100 for x in range(0, 100)]

# Critical (eigenvalue) run
scores = ['flux', 'heating', 'fission', '(n,Xt)']

spFS = openmc.StatePoint("./statepoint.100.h5")

flux_heat_fission_mesh = spFS.tallies[1].mean

tally_for_df = spFS.get_tally(id=1)
df_total = tally_for_df.get_pandas_dataframe(nuclides=False)
# df_flux = df_total[df_total['score'] == 'flux']
df_flux_FS = df_total[df_total['score'] == scores[0]]
mean_FS_R = df_flux_FS['mean'].values


df_avg = 0
df_avg_FS = np.zeros((100,1))

    
for i in range(1,100):
    df_avg_FS = df_avg_FS + df_flux_FS[df_flux_FS['mesh 1']['z'] == i]['mean'].values.reshape((100,1))

df_norm_FS = (df_avg_FS-np.min(df_avg_FS))/(np.max(df_avg_FS)-np.min(df_avg_FS))

normFS = mean_FS_R/np.max(mean_FS_R)

# rel_difference = df_norm_FS-df_norm_KEFF

# rel_difference_25 = np.zeros([25,1])
# for i in range(0,25):
#     for j in range(1,4):
#         rel_difference_25[i,0] = rel_difference_25[i,0] + rel_difference[4*i+j,0]

# relDiffAx.set_ylabel("Relative Difference Between Flux Distributions")

normAx.plot(x_ax_geom, df_norm_FS, label="$k_{eff}=$")
# normAx.plot(x_ax_geom, df_norm_FS, label="k_eff="+j)
# relDiffAx.plot(rel_difference)

# plt.legend(loc='lower right')
plt.show()




