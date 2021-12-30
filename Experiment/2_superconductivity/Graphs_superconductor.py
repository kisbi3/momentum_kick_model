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

#초전도성, 완전 도체 확인 그래프

# super_T_1=np.loadtxt('SuperconductorCT2RT.txt',delimiter='\t',usecols=[3],skiprows=5, encoding = 'UTF-8')
# super_V_1=np.loadtxt('SuperconductorCT2RT.txt',delimiter='\t',usecols=[2],skiprows=5, encoding = 'UTF-8')

# super_T_2=np.loadtxt('SuperconductorRT2CT.txt',delimiter='\t',usecols=[3],skiprows=5, encoding = 'UTF-8')
# super_V_2=np.loadtxt('SuperconductorRT2CT.txt',delimiter='\t',usecols=[2],skiprows=5, encoding = 'UTF-8')

# plt.plot(super_T_1, super_V_1 , color = 'blue' , linewidth=2, linestyle = '-', label=fr'$low \: to \: high$')
# plt.plot(super_T_2, super_V_2, color = 'red' , linewidth=2, linestyle = '-', label=fr'$high \: to \: low$')

# plt.xlabel(r'$Temperature \,\, (^\circ C)$',size=50)
# plt.ylabel(r'$Voltage\,\,(V)$',size=50)

# plt.minorticks_on()
# # plt.yscale('log')
# # ax.axis([-0.1,3.1,0,0.1])
# # plt.xlim(0,4)
# # plt.ylim(0,7000)

# ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

# plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
# plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


# plt.grid(color='silver',linestyle='-',linewidth=2)
# plt.legend(fontsize=45,framealpha=False,loc='upper left')
# # plt.legend(fontsize=45)
# plt.tight_layout()

# fig.savefig('superconductor.png')

# fig.clear()


#thermocouple 확인 그래프

thermocouple_T_Fe=np.loadtxt('thermocouple.csv',delimiter=',',usecols=[0],skiprows=2, max_rows=72)
thermocouple_V_Fe=np.loadtxt('thermocouple.csv',delimiter=',',usecols=[1],skiprows=2, max_rows=72)

thermocouple_T_NiCr=np.loadtxt('thermocouple.csv',delimiter=',',usecols=[3],skiprows=2, max_rows=60)
thermocouple_V_NiCr=np.loadtxt('thermocouple.csv',delimiter=',',usecols=[4],skiprows=2, max_rows=60)

thermocouple_T_Cu=np.loadtxt('thermocouple.csv',delimiter=',',usecols=[6],skiprows=2, max_rows=77)
thermocouple_V_Cu=np.loadtxt('thermocouple.csv',delimiter=',',usecols=[7],skiprows=2, max_rows=77)

def func(x,a):
    return a*(x-14.2)
def func1(x,a):
    return a*(x-17)

popt1, pcov1 = curve_fit(f=func,xdata=thermocouple_T_Fe,ydata=thermocouple_V_Fe)
popt2, pcov2 = curve_fit(f=func,xdata=thermocouple_T_NiCr,ydata=thermocouple_V_NiCr)
popt3, pcov3 = curve_fit(f=func1,xdata=thermocouple_T_Cu,ydata=thermocouple_V_Cu)

fitx=np.arange(14.2,85,0.01)
fitx2=np.arange(17,85,0.01)

plt.plot(fitx, func(fitx, *popt1), color='magenta', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,Fe,\,\,Error={float(pcov1):.3E},\,\,({float(popt1):.3E})x$')
plt.plot(fitx, func(fitx, *popt2), color='lightgreen', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,NiCr,\,\,Error={float(pcov2):.3E},\,\,({float(popt2):.3E})x$')
plt.plot(fitx2, func1(fitx2, *popt3), color='cyan', linewidth=7,zorder=1, linestyle = '-', label=fr'$Fitting,\,\,Cu,\,\,Error={float(pcov3):.3E},\,\,({float(popt3):.3E})x$')

plt.plot(thermocouple_T_Fe, thermocouple_V_Fe , color = 'blue' , linewidth=5, linestyle = '-', label=fr'$Fe$')
plt.plot(thermocouple_T_NiCr, thermocouple_V_NiCr , color = 'red' , linewidth=5, linestyle = '-', label=fr'$NiCr$')
plt.plot(thermocouple_T_Cu, thermocouple_V_Cu , color = 'green' , linewidth=5, linestyle = '-', label=fr'$Cu$')


plt.xlabel(r'$Temperature \,\, (^\circ C)$',size=50)
plt.ylabel(r'$Voltage\,\,(mV)$',size=50)

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

fig.savefig('thermocouple.png')