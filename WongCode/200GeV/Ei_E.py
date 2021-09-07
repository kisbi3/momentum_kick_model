import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

mpl.rcParams["text.usetex"] = True

fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

pt=np.loadtxt('Ei_E.csv',delimiter=',',usecols=[0],skiprows=1)
Ei=np.loadtxt('Ei_E.csv',delimiter=',',usecols=[1],skiprows=1)
Ef=np.loadtxt('Ei_E.csv',delimiter=',',usecols=[2],skiprows=1)
Ei_E=np.loadtxt('Ei_E.csv',delimiter=',',usecols=[3],skiprows=1)
E_Ei=np.loadtxt('Ei_E.csv',delimiter=',',usecols=[4],skiprows=1)


plt.xlabel(r'$p_{t}\quad(GeV)$',size=50)
plt.ylabel(r'$(1/N_{trig})dN_{ch}/p_{t}dp_{t}{}\quad(GeV^{-2})$',size=50)

plt.plot(pt,Ei, color="green",linewidth=7,linestyle = '-',label=r'$E_{i}$')
plt.plot(pt,Ef, color="blue",linewidth=7,linestyle = '-',label=r'$E_{f}$')
plt.plot(pt,Ei_E, color="red",linewidth=7,linestyle = '-',label=r'$E_{i}/E_{f}$')
plt.plot(pt,E_Ei, color="red",linewidth=7,linestyle = '--',label=r'$E_{f}/E_{i}$')

plt.minorticks_on()
plt.yscale('log')
ax.axis([0,4,0.1,100])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
plt.tight_layout()

fig.savefig('Ei_E.png')
# fig.savefig('Initial_Parton.png')
# plt.show()

fig.clear()

pt=np.loadtxt('Ei_E_phi.csv',delimiter=',',usecols=[0],skiprows=1)
Ef=np.loadtxt('Ei_E_phi.csv',delimiter=',',usecols=[1],skiprows=1)
Ei0=np.loadtxt('Ei_E_phi.csv',delimiter=',',usecols=[2],skiprows=1)
Ei5=np.loadtxt('Ei_E_phi.csv',delimiter=',',usecols=[3],skiprows=1)
Ei1=np.loadtxt('Ei_E_phi.csv',delimiter=',',usecols=[4],skiprows=1)
Eipi=np.loadtxt('Ei_E_phi.csv',delimiter=',',usecols=[5],skiprows=1)
Eipi2=np.loadtxt('Ei_E_phi.csv',delimiter=',',usecols=[6],skiprows=1)
Eipi4=np.loadtxt('Ei_E_phi.csv',delimiter=',',usecols=[7],skiprows=1)


plt.plot(pt,Ef, color="black",linewidth=7,linestyle = '--',label=r'$E_{f}$')
plt.plot(pt,Ei0, color="red",linewidth=7,linestyle = '-',label=r'$E_{i}\quad(phi=0.0)$')
plt.plot(pt,Ei5, color="green",linewidth=7,linestyle = '-',label=r'$E_{i}\quad(phi=0.5)$')
plt.plot(pt,Ei1, color="blue",linewidth=7,linestyle = '-',label=r'$E_{i}\quad(phi=1.0)$')
plt.plot(pt,Eipi, color="black",linewidth=7,linestyle = '-',label=r'$E_{i}\quad(phi=\pi)$')
plt.plot(pt,Eipi2, color="aqua",linewidth=7,linestyle = '-',label=r'$E_{i}\quad(phi=\pi/2)$')
plt.plot(pt,Eipi4, color="gold",linewidth=7,linestyle = '-',label=r'$E_{i}\quad(phi=\pi/4)$')

plt.xlabel(r'$p_{tf}\quad(GeV)$',size=50)
plt.ylabel(r'$E$',size=50)

plt.minorticks_on()
plt.yscale('log')
ax.axis([0,4,0.1,100])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
plt.tight_layout()

fig.savefig('Ei_E_phi.png')