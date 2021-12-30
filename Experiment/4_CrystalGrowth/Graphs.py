import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import csv
import matplotlib.ticker as ticker

from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit


batio3_angle=np.loadtxt('BaTio3.csv',delimiter=',',usecols=[0],skiprows=26)
batio3_inten=np.loadtxt('BaTio3.csv',delimiter=',',usecols=[1],skiprows=26)

catio3_angle=np.loadtxt('CaTio3.csv',delimiter=',',usecols=[0],skiprows=26)
catio3_inten=np.loadtxt('CaTio3.csv',delimiter=',',usecols=[1],skiprows=26)

srtio3_angle=np.loadtxt('SrTiO3.csv',delimiter=',',usecols=[0],skiprows=26)
srtio3_inten=np.loadtxt('SrTiO3.csv',delimiter=',',usecols=[1],skiprows=26)

batio3_Tem=np.loadtxt('BaTiO3_phase.csv',delimiter=',',usecols=[0],skiprows=1)
batio3_Cap=np.loadtxt('BaTiO3_phase.csv',delimiter=',',usecols=[1],skiprows=1)

# print(catio3_inten)


mpl.rcParams["text.usetex"] = True


fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)


plt.plot(batio3_Tem, batio3_Cap, color = 'blue', linewidth=4, label=fr'$BaTiO_3$')
# plt.scatter(batio3_angle, batio3_inten, color = 'blue', s = 20, label=fr'$BaTiO_3$')
# plt.yscale('log')
ax.minorticks_on()
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
plt.xlim(30,140)
# ax.set_xlim(10,90)
# ax.set_ylim(0,20)
ax.set_xlabel(r'$Temperature \, (^\circ C)$',size=50)
ax.set_ylabel(r'$Capacity (pF)$',size=50)
# ax.legend(fontsize=45,framealpha=False, loc = 'upper left')
plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.tight_layout()

fig.savefig('BaTiO3_phase.png')

fig.clear()

# plt.plot(batio3_angle, batio3_inten, color = 'blue', linewidth=4, label=fr'$BaTiO_3$')
plt.scatter(batio3_angle, batio3_inten, color = 'blue', s = 20, label=fr'$BaTiO_3$')
plt.yscale('log')
ax.minorticks_on()
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
plt.xlim(10,90)
ax.set_xlim(10,90)
# ax.set_ylim(0,20)
ax.set_xlabel(r'$2Theta \, (^\circ)$',size=50)
ax.set_ylabel(r'$counts$',size=50)
# ax.legend(fontsize=45,framealpha=False, loc = 'upper left')
plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.tight_layout()

fig.savefig('BaTiO3.png')

fig.clear()


# plt.plot(catio3_angle, catio3_inten, color = 'blue', linewidth=4, label=fr'$CaTiO_3$')
plt.scatter(catio3_angle, catio3_inten, color = 'blue', s = 20, label=fr'$CaTiO_3$')
plt.yscale('log')
ax.minorticks_on()
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
plt.xlim(10,90)
ax.set_xlim(10,90)
# ax.set_ylim(0,20)
ax.set_xlabel(r'$2Theta \, (^\circ)$',size=50)
ax.set_ylabel(r'$counts$',size=50)
# ax.legend(fontsize=45,framealpha=False, loc = 'upper left')
plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.tight_layout()

fig.savefig('CaTiO3.png')

fig.clear()

# plt.plot(srtio3_angle, srtio3_inten, color = 'blue', linewidth=4, label=fr'$SrTiO_3$')
plt.scatter(srtio3_angle, srtio3_inten, color = 'blue', s = 20, label=fr'$SrTiO_3$')
plt.yscale('log')
ax.minorticks_on()
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
plt.xlim(10,90)
ax.set_xlim(10,90)
# ax.set_ylim(0,20)
ax.set_xlabel(r'$2Theta \, (^\circ)$',size=50)
ax.set_ylabel(r'$counts$',size=50)
# ax.legend(fontsize=45,framealpha=False, loc = 'upper left')
plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.tight_layout()

fig.savefig('SrTiO3.png')

fig.clear()

# plt.plot(batio3_angle, batio3_inten, color = 'red', linewidth=4, label=fr'$BaTiO_3$')
# plt.plot(catio3_angle, catio3_inten, color = 'green', linewidth=4, label=fr'$CaTiO_3$')
# plt.plot(srtio3_angle, srtio3_inten, color = 'blue', linewidth=4, label=fr'$SrTiO_3$')

plt.scatter(batio3_angle, batio3_inten, color = 'red', s = 20, label=fr'$BaTiO_3$')
plt.scatter(catio3_angle, catio3_inten, color = 'green', s = 20, label=fr'$CaTiO_3$')
plt.scatter(srtio3_angle, srtio3_inten, color = 'blue', s = 20, label=fr'$SrTiO_3$')

ax.minorticks_on()
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
plt.xlim(10,90)
ax.set_xlim(10,90)
plt.yscale('log')
# ax.set_ylim(0,20)
ax.set_xlabel(r'$2Theta \, (^\circ)$',size=50)
ax.set_ylabel(r'$counts$',size=50)
plt.legend(fontsize=60,framealpha=False)
plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.tight_layout()

fig.savefig('All.png')