import csv
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt

# path='/content/drive/MyDrive/chtcRunsParaOpenmc5_30_23/referencePWR_SNF_comp/paraFuelMany2/output.txt'
  
path=f"keff.txt"

results = []
with open(path) as csvfile:
    multData = csv.reader(csvfile)
    for row in multData:
        # print(','.join(row))
        results.append(row)

npResults = np.asarray(results).astype(float)

keffs = npResults[:,1]
ann_heights = npResults[:,0]
ann_height_errs = npResults[:,2]

coreDiam = np.multiply(ann_heights, 1.5)

# temp_ax.set_title("Multiplier Radius vs Eta")
figure, ax = plt.subplots()
ax.set_xlabel("Core Outer Diameter")
ax.set_ylabel(r"Multiplication Factor ($k_{eff}$)")
eta_lines = plt.plot(coreDiam, keffs, 'c')
figure.savefig(fname=f"sweep_keff.png")

plt.axhline(y=0.95, color='grey', linestyle='--', label=r"$k_{eff}$=0.95")
plt.legend()
plt.show()
