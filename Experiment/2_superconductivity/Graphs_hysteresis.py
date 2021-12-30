import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import csv
import matplotlib.ticker as ticker

from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
# from scipy.optimize import curve_fit


mpl.rcParams["text.usetex"] = True


fig = plt.figure()
ax = plt.axes()


fig.set_size_inches(35, 16.534, forward=True)

hys2V02Hz_I=np.loadtxt('Hysteresis2V0.2Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys2V02Hz_B=np.loadtxt('Hysteresis2V0.2Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys2V04Hz_I=np.loadtxt('Hysteresis2V0.4Hz.txt',delimiter='\t',usecols=[2],skiprows=5)
hys2V04Hz_B=np.loadtxt('Hysteresis2V0.4Hz.txt',delimiter='\t',usecols=[6],skiprows=5)
hys2V06Hz_I=np.loadtxt('Hysteresis2V0.6Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys2V06Hz_B=np.loadtxt('Hysteresis2V0.6Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys2V08Hz_I=np.loadtxt('Hysteresis2V0.8Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys2V08Hz_B=np.loadtxt('Hysteresis2V0.8Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys2V10Hz_I=np.loadtxt('Hysteresis2V1.0Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys2V10Hz_B=np.loadtxt('Hysteresis2V1.0Hz.txt',delimiter='\t',usecols=[5],skiprows=5)

hys4V02Hz_I=np.loadtxt('Hysteresis4V0.2Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys4V02Hz_B=np.loadtxt('Hysteresis4V0.2Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys4V04Hz_I=np.loadtxt('Hysteresis4V0.4Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys4V04Hz_B=np.loadtxt('Hysteresis4V0.4Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys4V06Hz_I=np.loadtxt('Hysteresis4V0.6Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys4V06Hz_B=np.loadtxt('Hysteresis4V0.6Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys4V08Hz_I=np.loadtxt('Hysteresis4V0.8Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys4V08Hz_B=np.loadtxt('Hysteresis4V0.8Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys4V10Hz_I=np.loadtxt('Hysteresis4V1.0Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys4V10Hz_B=np.loadtxt('Hysteresis4V1.0Hz.txt',delimiter='\t',usecols=[5],skiprows=5)

hys6V02Hz_I=np.loadtxt('Hysteresis6V0.2Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys6V02Hz_B=np.loadtxt('Hysteresis6V0.2Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys6V04Hz_I=np.loadtxt('Hysteresis6V0.4Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys6V04Hz_B=np.loadtxt('Hysteresis6V0.4Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys6V06Hz_I=np.loadtxt('Hysteresis6V0.6Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys6V06Hz_B=np.loadtxt('Hysteresis6V0.6Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys6V08Hz_I=np.loadtxt('Hysteresis6V0.8Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys6V08Hz_B=np.loadtxt('Hysteresis6V0.8Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys6V10Hz_I=np.loadtxt('Hysteresis6V1.0Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys6V10Hz_B=np.loadtxt('Hysteresis6V1.0Hz.txt',delimiter='\t',usecols=[5],skiprows=5)

hys8V02Hz_I=np.loadtxt('Hysteresis8V0.2Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys8V02Hz_B=np.loadtxt('Hysteresis8V0.2Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys8V04Hz_I=np.loadtxt('Hysteresis8V0.4Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys8V04Hz_B=np.loadtxt('Hysteresis8V0.4Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys8V06Hz_I=np.loadtxt('Hysteresis8V0.6Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys8V06Hz_B=np.loadtxt('Hysteresis8V0.6Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys8V08Hz_I=np.loadtxt('Hysteresis8V0.8Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys8V08Hz_B=np.loadtxt('Hysteresis8V0.8Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys8V10Hz_I=np.loadtxt('Hysteresis8V1.0Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys8V10Hz_B=np.loadtxt('Hysteresis8V1.0Hz.txt',delimiter='\t',usecols=[5],skiprows=5)

hys10V02Hz_I=np.loadtxt('Hysteresis10V0.2Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys10V02Hz_B=np.loadtxt('Hysteresis10V0.2Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys10V04Hz_I=np.loadtxt('Hysteresis10V0.4Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys10V04Hz_B=np.loadtxt('Hysteresis10V0.4Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys10V06Hz_I=np.loadtxt('Hysteresis10V0.6Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys10V06Hz_B=np.loadtxt('Hysteresis10V0.6Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys10V08Hz_I=np.loadtxt('Hysteresis10V0.8Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys10V08Hz_B=np.loadtxt('Hysteresis10V0.8Hz.txt',delimiter='\t',usecols=[5],skiprows=5)
hys10V10Hz_I=np.loadtxt('Hysteresis10V1.0Hz.txt',delimiter='\t',usecols=[1],skiprows=5)
hys10V10Hz_B=np.loadtxt('Hysteresis10V1.0Hz.txt',delimiter='\t',usecols=[5],skiprows=5)

plt.plot(hys2V02Hz_I, hys2V02Hz_B, color = 'blue' , linewidth=1, linestyle = '-',label=fr'$2V, \, 0.2Hz$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Flux(Vs)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle='-',linewidth=2)
# plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hysteresis_2V_02Hz.png')

fig.clear()

plt.plot(hys2V02Hz_I, hys2V02Hz_B, color = 'blue' , linewidth=1, linestyle = '-',label=fr'$2V, \, 0.2Hz$')
plt.plot(hys2V04Hz_I, hys2V04Hz_B, color = 'red'  , linewidth=1, linestyle = '-',label=fr'$2V, \, 0.4Hz$')
plt.plot(hys2V06Hz_I, hys2V06Hz_B, color = 'green', linewidth=1, linestyle = '-',label=fr'$2V, \, 0.6Hz$')
plt.plot(hys2V08Hz_I, hys2V08Hz_B, color = 'black', linewidth=1, linestyle = '-',label=fr'$2V, \, 0.8Hz$')
plt.plot(hys2V10Hz_I, hys2V10Hz_B, color = 'cyan' , linewidth=1, linestyle = '-',label=fr'$2V, \, 1.0Hz$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Flux(Vs)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle='-',linewidth=2)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hysteresis_2V.png')

fig.clear()

plt.plot(hys4V02Hz_I, hys4V02Hz_B, color = 'blue' , linewidth=1, linestyle = '-',label=fr'$4V, \, 0.2Hz$')
plt.plot(hys4V04Hz_I, hys4V04Hz_B, color = 'red'  , linewidth=1, linestyle = '-',label=fr'$4V, \, 0.4Hz$')
plt.plot(hys4V06Hz_I, hys4V06Hz_B, color = 'green', linewidth=1, linestyle = '-',label=fr'$4V, \, 0.6Hz$')
plt.plot(hys4V08Hz_I, hys4V08Hz_B, color = 'black', linewidth=1, linestyle = '-',label=fr'$4V, \, 0.8Hz$')
plt.plot(hys4V10Hz_I, hys4V10Hz_B, color = 'cyan' , linewidth=1, linestyle = '-',label=fr'$4V, \, 1.0Hz$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Flux(Vs)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle='-',linewidth=2)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hysteresis_4V.png')

fig.clear()

plt.plot(hys6V02Hz_I, hys6V02Hz_B, color = 'blue' , linewidth=1, linestyle = '-',label=fr'$6V, \, 0.2Hz$')
plt.plot(hys6V04Hz_I, hys6V04Hz_B, color = 'red'  , linewidth=1, linestyle = '-',label=fr'$6V, \, 0.4Hz$')
plt.plot(hys6V06Hz_I, hys6V06Hz_B, color = 'green', linewidth=1, linestyle = '-',label=fr'$6V, \, 0.6Hz$')
plt.plot(hys6V08Hz_I, hys6V08Hz_B, color = 'black', linewidth=1, linestyle = '-',label=fr'$6V, \, 0.8Hz$')
plt.plot(hys6V10Hz_I, hys6V10Hz_B, color = 'cyan' , linewidth=1, linestyle = '-',label=fr'$6V, \, 1.0Hz$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Flux(Vs)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle='-',linewidth=2)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hysteresis_6V.png')

fig.clear()

plt.plot(hys8V02Hz_I, hys8V02Hz_B, color = 'blue' , linewidth=1, linestyle = '-',label=fr'$8V, \, 0.2Hz$')
plt.plot(hys8V04Hz_I, hys8V04Hz_B, color = 'red'  , linewidth=1, linestyle = '-',label=fr'$8V, \, 0.4Hz$')
plt.plot(hys8V06Hz_I, hys8V06Hz_B, color = 'green', linewidth=1, linestyle = '-',label=fr'$8V, \, 0.6Hz$')
plt.plot(hys8V08Hz_I, hys8V08Hz_B, color = 'black', linewidth=1, linestyle = '-',label=fr'$8V, \, 0.8Hz$')
plt.plot(hys8V10Hz_I, hys8V10Hz_B, color = 'cyan' , linewidth=1, linestyle = '-',label=fr'$8V, \, 1.0Hz$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Flux(Vs)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle='-',linewidth=2)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hysteresis_8V.png')

fig.clear()

plt.plot(hys10V02Hz_I, hys10V02Hz_B, color = 'blue' , linewidth=1, linestyle = '-',label=fr'$10V, \, 0.2Hz$')
plt.plot(hys10V04Hz_I, hys10V04Hz_B, color = 'red'  , linewidth=1, linestyle = '-',label=fr'$10V, \, 0.4Hz$')
plt.plot(hys10V06Hz_I, hys10V06Hz_B, color = 'green', linewidth=1, linestyle = '-',label=fr'$10V, \, 0.6Hz$')
plt.plot(hys10V08Hz_I, hys10V08Hz_B, color = 'black', linewidth=1, linestyle = '-',label=fr'$10V, \, 0.8Hz$')
plt.plot(hys10V10Hz_I, hys10V10Hz_B, color = 'cyan' , linewidth=1, linestyle = '-',label=fr'$10V, \, 1.0Hz$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Flux(Vs)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle='-',linewidth=2)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hysteresis_10V.png')

fig.clear()

plt.plot(hys2V02Hz_I, hys2V02Hz_B , color = 'blue' , linewidth=1, linestyle = '-',label=fr'$2V, \, 0.2Hz$')
plt.plot(hys4V02Hz_I, hys4V02Hz_B , color = 'red'  , linewidth=1, linestyle = '-',label=fr'$4V, \, 0.2Hz$')
plt.plot(hys6V02Hz_I, hys6V02Hz_B , color = 'green', linewidth=1, linestyle = '-',label=fr'$6V, \, 0.2Hz$')
plt.plot(hys8V02Hz_I, hys8V02Hz_B , color = 'black', linewidth=1, linestyle = '-',label=fr'$8V, \, 0.2Hz$')
plt.plot(hys10V02Hz_I, hys10V02Hz_B, color = 'cyan' , linewidth=1, linestyle = '-',label=fr'$10V, \, 0.2Hz$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Flux(Vs)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle='-',linewidth=2)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hysteresis_02Hz.png')

fig.clear()

plt.plot(hys2V04Hz_I, hys2V04Hz_B , color = 'blue' , linewidth=1, linestyle = '-',label=fr'$2V, \, 0.4Hz$')
plt.plot(hys4V04Hz_I, hys4V04Hz_B , color = 'red'  , linewidth=1, linestyle = '-',label=fr'$4V, \, 0.4Hz$')
plt.plot(hys6V04Hz_I, hys6V04Hz_B , color = 'green', linewidth=1, linestyle = '-',label=fr'$6V, \, 0.4Hz$')
plt.plot(hys8V04Hz_I, hys8V04Hz_B , color = 'black', linewidth=1, linestyle = '-',label=fr'$8V, \, 0.4Hz$')
plt.plot(hys10V04Hz_I, hys10V04Hz_B, color = 'cyan' , linewidth=1, linestyle = '-',label=fr'$10V, \, 0.4Hz$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Flux(Vs)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle='-',linewidth=2)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hysteresis_04Hz.png')

fig.clear()

plt.plot(hys2V06Hz_I, hys2V06Hz_B , color = 'blue' , linewidth=1, linestyle = '-', label=fr'$2V, \, 0.6Hz$')
plt.plot(hys4V06Hz_I, hys4V06Hz_B , color = 'red'  , linewidth=1, linestyle = '-', label=fr'$4V, \, 0.6Hz$')
plt.plot(hys6V06Hz_I, hys6V06Hz_B , color = 'green', linewidth=1, linestyle = '-', label=fr'$6V, \, 0.6Hz$')
plt.plot(hys8V06Hz_I, hys8V06Hz_B , color = 'black', linewidth=1, linestyle = '-', label=fr'$8V, \, 0.6Hz$')
plt.plot(hys10V06Hz_I,hys10V06Hz_B, color = 'cyan' , linewidth=1, linestyle = '-',label=fr'$10V, \, 0.6Hz$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Flux(Vs)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle='-',linewidth=2)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hysteresis_06Hz.png')

fig.clear()

plt.plot(hys2V08Hz_I, hys2V08Hz_B , color = 'blue' , linewidth=1, linestyle = '-', label=fr'$2V, \, 0.8Hz$')
plt.plot(hys4V08Hz_I, hys4V08Hz_B , color = 'red'  , linewidth=1, linestyle = '-', label=fr'$4V, \, 0.8Hz$')
plt.plot(hys6V08Hz_I, hys6V08Hz_B , color = 'green', linewidth=1, linestyle = '-', label=fr'$6V, \, 0.8Hz$')
plt.plot(hys8V08Hz_I, hys8V08Hz_B , color = 'black', linewidth=1, linestyle = '-', label=fr'$8V, \, 0.8Hz$')
plt.plot(hys10V08Hz_I,hys10V08Hz_B, color = 'cyan' , linewidth=1, linestyle = '-',label=fr'$10V, \, 0.8Hz$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Flux(Vs)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle='-',linewidth=2)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hysteresis_08Hz.png')

fig.clear()

plt.plot(hys2V10Hz_I, hys2V10Hz_B , color = 'blue' , linewidth=1, linestyle = '-', label=fr'$2V, \, 1.0Hz$')
plt.plot(hys4V10Hz_I, hys4V10Hz_B , color = 'red'  , linewidth=1, linestyle = '-', label=fr'$4V, \, 1.0Hz$')
plt.plot(hys6V10Hz_I, hys6V10Hz_B , color = 'green', linewidth=1, linestyle = '-', label=fr'$6V, \, 1.0Hz$')
plt.plot(hys8V10Hz_I, hys8V10Hz_B , color = 'black', linewidth=1, linestyle = '-', label=fr'$8V, \, 1.0Hz$')
plt.plot(hys10V10Hz_I,hys10V10Hz_B, color = 'cyan' , linewidth=1, linestyle = '-',label=fr'$10V, \, 1.0Hz$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Flux(Vs)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle='-',linewidth=2)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hysteresis_10Hz.png')

fig.clear()