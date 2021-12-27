import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import csv
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit

mpl.rcParams["text.usetex"] = True


fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

periph_276_phi=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[0],skiprows=1)
periph_276_data=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[1],skiprows=1)
periph_13_phi=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[2],skiprows=1)
periph_13_data=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[3],skiprows=1)


# plt.plot(periph_276_phi, periph_276_data, color = 'red', linewidth=7, linestyle = '-',label=fr'$peripheral, 2.76TeV$')
# plt.plot(periph_13_phi, periph_13_data, color = 'blue', linewidth=7, linestyle = '-',label=fr'$peripheral, 13TeV$')

plt.scatter(periph_276_phi, periph_276_data, color = 'red', s=100,label=fr'$peripheral, 2.76TeV$')
plt.scatter(periph_13_phi, periph_13_data, color = 'blue', s=100,label=fr'$peripheral, 13TeV$')

plt.ylabel(r'$Y(\Delta \phi)$', size = 70)
plt.xlabel(r'$\Delta \phi$', size=70)

plt.title(r'$ATLAS,\quad Y^{periph} \, (0 \leq N^{rec}_{ch} < 20)$', fontsize = 75)

plt.minorticks_on()

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)

plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('peripheral.png')

fig.clear()



