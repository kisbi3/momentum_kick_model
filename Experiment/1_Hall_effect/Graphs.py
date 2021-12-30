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

mag=np.loadtxt('Data.csv',delimiter=',',usecols=[0],skiprows=3,max_rows=67)
elect_mag=np.loadtxt('Data.csv',delimiter=',',usecols=[1],skiprows=3,max_rows=67)
elect_mag_double=np.loadtxt('Data.csv',delimiter=',',usecols=[2],skiprows=3,max_rows=67)

plt.plot(elect_mag, mag, color = 'blue', linewidth=7, linestyle = '-',label=fr'$Current_1$')
plt.plot(elect_mag_double, mag, color = 'red', linewidth=7, linestyle = '-',label=fr'$Current_2$')

plt.xlabel(r'$Current(A)$',size=50)
plt.ylabel(r'$Magnetic\,\,Field(G)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.xlim(0,4)
plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Magnetic_Field.png')

fig.clear()

plt.plot(elect_mag, mag*elect_mag/elect_mag_double, color = 'blue', linewidth=7, linestyle = '-',label=fr'$Average\,\,Current$',zorder=2)

# def func(x,a,b):
#     return a*x+b

# popt1, pcov1 = curve_fit(f=func,xdata=elect_mag_double,ydata=mag)

# print(popt1)

# fitx=np.arange(0,4,0.01)
# plt.plot(fitx, func(fitx, *popt1), color='red', linewidth=7,zorder=1, linestyle = '--', label=fr'$Fitting-double\,\,check\,\,Current$')

plt.ylabel(r'$Magnetic\,\,Field(G)$',size=50)
plt.xlabel(r'$Current(A)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('doublecheck.png')



fig.clear()

plt.plot(mag, elect_mag-elect_mag_double, color = 'blue', linewidth=7, linestyle = '-',label=fr'$Current_1$')

plt.xlabel(r'$Magnetic\,\,Field(G)$',size=50)
plt.ylabel(r'$Current_1-Current_2$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Magnetic_Deviation.png')

fig.clear()

plt.plot(mag, ((elect_mag-elect_mag_double)/elect_mag)*100, color = 'blue', linewidth=7, linestyle = '-',label=fr'$Current_1$')

plt.xlabel(r'$Magnetic\,\,Field(G)$',size=50)
plt.ylabel(r'$relative\,error$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(0,4)
# plt.ylim(0,7000)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Magnetic_Deviation_percent.png')

fig.clear()

mag=np.loadtxt('Data.csv',delimiter=',',usecols=[5],skiprows=4,max_rows=30)
ntype_2mA=np.loadtxt('Data.csv',delimiter=',',usecols=[7],skiprows=4,max_rows=30)
ntype_4mA=np.loadtxt('Data.csv',delimiter=',',usecols=[12],skiprows=4,max_rows=30)
mag_6mA=np.loadtxt('Data.csv',delimiter=',',usecols=[15],skiprows=4,max_rows=23)
ntype_6mA=np.loadtxt('Data.csv',delimiter=',',usecols=[17],skiprows=4,max_rows=23)

plt.plot(mag, ntype_2mA, color = 'red', linewidth=7,zorder=2, linestyle = '--',label=fr'$n-type,\,\, 2mA$')
plt.plot(mag, ntype_4mA, color = 'green', linewidth=7,zorder=2, linestyle = '--',label=fr'$n-type,\,\, 4mA$')
plt.plot(mag_6mA, ntype_6mA, color = 'blue', linewidth=7,zorder=2, linestyle = '--',label=fr'$n-type,\,\, 6mA$')

def func(x,a):
    return a*x

popt1, pcov1 = curve_fit(f=func,xdata=mag,ydata=ntype_2mA)
popt2, pcov2 = curve_fit(f=func,xdata=mag,ydata=ntype_4mA)
popt3, pcov3 = curve_fit(f=func,xdata=mag_6mA,ydata=ntype_6mA)

print('n-type 2mA density',(-10**20*2*10)/(8.01*popt1))
print('n-type 4mA density',(-10**20*4*10)/(8.01*popt2))
print('n-type 6mA density',(-10**20*6*10)/(8.01*popt3))

fitx=np.arange(0,6000,0.01)
# print(popt1)
plt.plot(fitx, func(fitx, *popt1), color='magenta', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,2mA,\,\,Error={float(pcov1):.3E},\,\,({float(popt1):.3E})x$')
plt.plot(fitx, func(fitx, *popt2), color='lightgreen', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,4mA,\,\,Error={float(pcov2):.3E},\,\,({float(popt2):.3E})x$')
plt.plot(fitx, func(fitx, *popt3), color='cyan', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,6mA,\,\,Error={float(pcov3):.3E},\,\,({float(popt3):.3E})x$')



plt.xlabel(r'$Magnetic\,\,Field(G)$',size=50)
plt.ylabel(r'$Hall\,\,Voltage(mV)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.xlim(0,6000)
plt.ylim(-250,0)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False,loc='lower left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('n-type.png')

fig.clear()

fitx=np.arange(0,6000,0.01)
# print(popt1)
plt.plot(mag, ((ntype_2mA-func(mag, *popt1))/ntype_2mA)*100, color='red', linewidth=7,zorder=1, linestyle = '-', label=fr'$2mA$')
plt.plot(mag, ((ntype_4mA-func(mag, *popt2))/ntype_4mA)*100, color='green', linewidth=7,zorder=1, linestyle = '-', label=fr'$4mA$')
plt.plot(mag_6mA, ((ntype_6mA-func(mag_6mA, *popt3))/ntype_6mA)*100, color='blue', linewidth=7,zorder=1, linestyle = '-', label=fr'$6mA$')



plt.xlabel(r'$Magnetic\,\,Field(G)$',size=50)
plt.ylabel(r'$Voltage\,\,Difference$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.xlim(0,6000)
# plt.ylim(-250,0)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False,loc='lower left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hall_Voltage_Error_ntype.png')

fig.clear()

mag_6mA_2 = np.zeros(len(mag_6mA))

for i in range(len(mag_6mA)):
    mag_6mA_2[i] = mag_6mA[i]*elect_mag[2*i+1]/elect_mag_double[2*i+1]
double_mag_1 = np.zeros(len(mag))
for i in range(len(mag)):
    double_mag_1[i] = mag[i]*elect_mag[2*i+1]/elect_mag_double[2*i+1]
    # print(mag[i],elect_mag[2*i+1],elect_mag_double[2*i+1])

plt.plot(double_mag_1, ntype_2mA, color = 'red', linewidth=7,zorder=2, linestyle = '--',label=fr'$n-type,\,\, 2mA$')
plt.plot(double_mag_1, ntype_4mA, color = 'green', linewidth=7,zorder=2, linestyle = '--',label=fr'$n-type,\,\, 4mA$')
plt.plot(mag_6mA_2, ntype_6mA, color = 'blue', linewidth=7,zorder=2, linestyle = '--',label=fr'$n-type,\,\, 6mA$')

popt1, pcov1 = curve_fit(f=func,xdata=double_mag_1,ydata=ntype_2mA)
popt2, pcov2 = curve_fit(f=func,xdata=double_mag_1,ydata=ntype_4mA)
popt3, pcov3 = curve_fit(f=func,xdata=mag_6mA_2,ydata=ntype_6mA)

print('n-type 2mA density 2',(-10**20*2*10)/(8.01*popt1))
print('n-type 4mA density 2',(-10**20*4*10)/(8.01*popt2))
print('n-type 6mA density 2',(-10**20*6*10)/(8.01*popt3))

fitx=np.arange(0,6000,0.01)
# print(popt1)
plt.plot(fitx, func(fitx, *popt1), color='magenta', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,2mA,\,\,Error={float(pcov1):.3E},\,\,({float(popt1):.3E})x$')
plt.plot(fitx, func(fitx, *popt2), color='lightgreen', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,4mA,\,\,Error={float(pcov2):.3E},\,\,({float(popt2):.3E})x$')
plt.plot(fitx, func(fitx, *popt3), color='cyan', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,6mA,\,\,Error={float(pcov3):.3E},\,\,({float(popt3):.3E})x$')



plt.xlabel(r'$Magnetic\,\,Field(G)$',size=50)
plt.ylabel(r'$Hall\,\,Voltage(mV)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.xlim(0,6000)
plt.ylim(-250,0)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False,loc='lower left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('n-type_2.png')

fig.clear()

mag=np.loadtxt('Data.csv',delimiter=',',usecols=[21],skiprows=4,max_rows=30)
ptype_2mA=np.loadtxt('Data.csv',delimiter=',',usecols=[23],skiprows=4,max_rows=30)
ptype_4mA=np.loadtxt('Data.csv',delimiter=',',usecols=[28],skiprows=4,max_rows=30)
ptype_6mA=np.loadtxt('Data.csv',delimiter=',',usecols=[33],skiprows=4,max_rows=30)

plt.plot(mag, ptype_2mA, color = 'red', linewidth=7,zorder=2, linestyle = '--',label=fr'$p-type,\,\, 2mA$')
plt.plot(mag, ptype_4mA, color = 'green', linewidth=7,zorder=2, linestyle = '--',label=fr'$p-type,\,\, 4mA$')
plt.plot(mag, ptype_6mA, color = 'blue', linewidth=7,zorder=2, linestyle = '--',label=fr'$p-type,\,\, 6mA$')

def func(x,a):
    return a*x

popt1, pcov1 = curve_fit(f=func,xdata=mag,ydata=ptype_2mA)
popt2, pcov2 = curve_fit(f=func,xdata=mag,ydata=ptype_4mA)
popt3, pcov3 = curve_fit(f=func,xdata=mag,ydata=ptype_6mA)

popt4, pcov4 = curve_fit(f=func,xdata=mag,ydata=ptype_2mA)
popt5, pcov5 = curve_fit(f=func,xdata=mag,ydata=ptype_4mA)
popt6, pcov6 = curve_fit(f=func,xdata=mag,ydata=ptype_6mA)

print('p-type 2mA density',(-10**20*2*10)/(8.01*popt1))
print('p-type 4mA density',(-10**20*4*10)/(8.01*popt2))
print('p-type 6mA density',(-10**20*6*10)/(8.01*popt3))

fitx=np.arange(0,6000,0.01)
# print(popt1)
plt.plot(fitx, func(fitx, *popt1), color='magenta', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,2mA,\,\,Error={float(pcov1):.3E},\,\,({float(popt1):.3E})x$')
plt.plot(fitx, func(fitx, *popt2), color='lightgreen', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,4mA,\,\,Error={float(pcov2):.3E},\,\,({float(popt2):.3E})x$')
plt.plot(fitx, func(fitx, *popt3), color='cyan', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,6mA,\,\,Error={float(pcov3):.3E},\,\,({float(popt3):.3E})x$')


plt.xlabel(r'$Magnetic\,\,Field(G)$',size=50)
plt.ylabel(r'$Hall\,\,Voltage(mV)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.xlim(0,6000)
plt.ylim(0,140)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('p-type.png')

fig.clear()

fitx=np.arange(0,6000,0.01)
# print(popt1)
plt.plot(mag, ((ptype_2mA-func(mag, *popt1))/ptype_2mA)*100, color='red', linewidth=7,zorder=1, linestyle = '-', label=fr'$2mA$')
plt.plot(mag, ((ptype_4mA-func(mag, *popt2))/ptype_4mA)*100, color='green', linewidth=7,zorder=1, linestyle = '-', label=fr'$4mA$')
plt.plot(mag, ((ptype_6mA-func(mag, *popt3))/ptype_6mA)*100, color='blue', linewidth=7,zorder=1, linestyle = '-', label=fr'$6mA$')



plt.xlabel(r'$Magnetic\,\,Field(G)$',size=50)
plt.ylabel(r'$Voltage\,\,Difference$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.xlim(0,6000)
# plt.ylim(-250,0)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False,loc='lower left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Hall_Voltage_Error_ptype.png')

fig.clear()

double_mag_1 = np.zeros(len(mag))
for i in range(len(mag)):
    double_mag_1[i] = mag[i]*elect_mag[2*i+1]/elect_mag_double[2*i+1]
    # print(mag[i],elect_mag[2*i+1],elect_mag_double[2*i+1])

plt.plot(double_mag_1, ptype_2mA, color = 'red', linewidth=7,zorder=2, linestyle = '--',label=fr'$p-type,\,\, 2mA$')
plt.plot(double_mag_1, ptype_4mA, color = 'green', linewidth=7,zorder=2, linestyle = '--',label=fr'$p-type,\,\, 4mA$')
plt.plot(double_mag_1, ptype_6mA, color = 'blue', linewidth=7,zorder=2, linestyle = '--',label=fr'$p-type,\,\, 6mA$')


print('p-type 2mA density 2',(-10**20*2*10)/(8.01*popt1))
print('p-type 4mA density 2',(-10**20*4*10)/(8.01*popt2))
print('p-type 6mA density 2',(-10**20*6*10)/(8.01*popt3))

fitx=np.arange(0,6000,0.01)
# print(popt1)
plt.plot(fitx, func(fitx, *popt1), color='magenta', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,2mA,\,\,Error={float(pcov1):.3E},\,\,({float(popt1):.3E})x$')
plt.plot(fitx, func(fitx, *popt2), color='lightgreen', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,4mA,\,\,Error={float(pcov2):.3E},\,\,({float(popt2):.3E})x$')
plt.plot(fitx, func(fitx, *popt3), color='cyan', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,6mA,\,\,Error={float(pcov3):.3E},\,\,({float(popt3):.3E})x$')



plt.xlabel(r'$Magnetic\,\,Field(G)$',size=50)
plt.ylabel(r'$Hall\,\,Voltage(mV)$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.xlim(0,6000)
plt.ylim(0,140)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('p-type_2.png')
