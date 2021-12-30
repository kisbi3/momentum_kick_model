import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import csv
import matplotlib.ticker as ticker

from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit

Si_thick=0.00045
Al_thick=0.000018
W_thick=0.00006


Si_Cur_Tem=np.loadtxt('data.csv',delimiter=',',usecols=[0],skiprows=4,max_rows=20)
Si_Cur_Cur=np.loadtxt('data.csv',delimiter=',',usecols=[1],skiprows=4,max_rows=20)
Si_Cur_Vol=np.loadtxt('data.csv',delimiter=',',usecols=[2],skiprows=4,max_rows=20)

Si_Tem_Tem=np.loadtxt('data.csv',delimiter=',',usecols=[3],skiprows=4,max_rows=18)
Si_Tem_Cur=np.loadtxt('data.csv',delimiter=',',usecols=[4],skiprows=4,max_rows=18)
Si_Tem_Vol=np.loadtxt('data.csv',delimiter=',',usecols=[5],skiprows=4,max_rows=18)

Al_Cur_Tem=np.loadtxt('data.csv',delimiter=',',usecols=[7],skiprows=4,max_rows=20)
Al_Cur_Cur=np.loadtxt('data.csv',delimiter=',',usecols=[8],skiprows=4,max_rows=20)
Al_Cur_Vol=np.loadtxt('data.csv',delimiter=',',usecols=[9],skiprows=4,max_rows=20)

W_Cur_Tem=np.loadtxt('data.csv',delimiter=',',usecols=[11],skiprows=4,max_rows=10)
W_Cur_Cur=np.loadtxt('data.csv',delimiter=',',usecols=[12],skiprows=4,max_rows=10)
W_Cur_Vol=np.loadtxt('data.csv',delimiter=',',usecols=[13],skiprows=4,max_rows=10)

# print(Si_Cur_Tem)
# print(Si_Cur_Cur)
# print(Si_Cur_Vol)

# print(Si_Cur_Tem_2)
# print(Si_Cur_Cur_2)
# print(Si_Cur_Vol_2)

# print(Al_Cur_Tem)
# print(Al_Cur_Cur)
# print(Al_Cur_Vol)

# print(W_Cur_Tem)
# print(W_Cur_Cur)
# print(W_Cur_Vol)


mpl.rcParams["text.usetex"] = True


# fig = plt.figure()
ax = plt.axes()
# ax1 = plt.axes()
# ax2 = ax1.twinx()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()


fig.set_size_inches(35, 16.534, forward=True)


def func(x,a):
    return a*x

popt1, pcov1 = curve_fit(f=func,xdata=Si_Cur_Cur,ydata=Si_Cur_Vol)
popt2, pcov2 = curve_fit(f=func,xdata=Al_Cur_Cur,ydata=Al_Cur_Vol)
popt3, pcov3 = curve_fit(f=func,xdata=W_Cur_Cur,ydata=W_Cur_Vol)

def func1(x,a,b):
    return a*x+b
popt4, pcov4 = curve_fit(f=func1,xdata=Si_Tem_Tem,ydata=Si_Tem_Vol)

s = 0.002       # 2mm
D = 0.01        # 1cm


def factor(alpha, T, t):
    f1 = np.log(np.sinh(t/s)/np.sinh(t/(2*s)))
    f2 = np.log(2)+np.log((((D/s)*(D/s))+3)/(((D/s)*(D/s))-3))
    return (np.pi*np.log(2)*(1+alpha*(T-20)))/(f1*f2)

factor_Si_Cur=factor(0, Si_Cur_Tem, Si_thick)

factor_Al=factor(0.004308, Al_Cur_Tem, Al_thick)
factor_W=factor(0.004403, W_Cur_Tem, W_thick)

R_Si_Cur = (factor_Si_Cur[0])*popt1*Si_thick       #알 수 없으므로 alpha=1로 고정
R_Al_Cur = (factor_Al)*popt2*Al_thick
R_W_Cur = (factor_W)*popt3*W_thick

print(f'Si resistance : {R_Si_Cur}')
print(f'Al resistance : {np.mean(R_Al_Cur)}')
print(f'W resistance : {np.mean(R_W_Cur)}')

fitx=np.arange(0,210,0.01)


ax1.plot(Si_Cur_Cur, Si_Cur_Vol, color = 'blue', linewidth=7, linestyle = '-',label=fr'$P \, type \, Si \,\, (B \, high \, doped)$')
ax1.plot(fitx, func(fitx, *popt1), color='cyan', linewidth=7,zorder=1, linestyle = ':', label=fr'$Fitting\,Si,\,\,Error={float(pcov1):.3E},\,\,({float(popt1):.3E})x$')
ax1.minorticks_on()
ax1.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
ax1.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
ax1.set_xlim(0,210)
ax1.set_ylim(0,20)
ax1.set_xlabel(r'$Current \,\, [mA]$',size=50)
ax1.set_ylabel(r'$Voltage \,\, [mV] \quad (Si)$',size=50)
ax1.legend(fontsize=45,framealpha=False, loc = 'upper left')
ax1.grid(color='silver',linestyle=':',linewidth=3)
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])



ax2.plot(Al_Cur_Cur, Al_Cur_Vol, color = 'red', linewidth=7, linestyle = '-',label=fr'$Al \, Foil$')
ax2.plot(W_Cur_Cur, W_Cur_Vol, color = 'green', linewidth=7, linestyle = '-',label=fr'$W \, Foil$')
ax2.plot(fitx, func(fitx, *popt2), color='magenta', linewidth=7,zorder=1, linestyle = ':', label=fr'$Fitting\,Al,\,\,Error={float(pcov2):.3E},\,\,({float(popt2):.3E})x$')
ax2.plot(fitx, func(fitx, *popt3), color='lightgreen', linewidth=7,zorder=1, linestyle = ':', label=fr'$Fitting\,W,\,\,Error={float(pcov3):.3E},\,\,({float(popt3):.3E})x$')
ax2.set_ylabel(r'$Voltage \,\, [mV] \quad (Al, \, W)$',size=50)
ax2.minorticks_on()

ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))
ax2.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
ax2.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)
ax2.set_xlim(0,210)
ax2.set_ylim(0,0.06)

# ax2.grid(color='silver',linestyle=':',linewidth=3)

ax2.legend(fontsize=45,framealpha=False, loc = 'lower right')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Cur_vs_Vol.png')

fig.clear()

# print(pcov4)
# print(popt4)

err=pcov4[0]
a=popt4[0]
b=popt4[1]

fitting = fitx*a+b
# print(fitting)
print(f'error = {err}')
print(f'a = {a}')
print(f'b = {b}')
factor_Si_Tem=factor(0, Si_Tem_Tem, Si_thick)
alpha  = (1-(Si_Tem_Cur*0.000018)/(Si_Tem_Vol*factor_Si_Tem[0]))/(Si_Tem_Tem-20)
factor_Si_Cur=factor(np.mean(alpha), Si_Cur_Tem, Si_thick)
Si_Cur_Cur = np.delete(Si_Cur_Cur,0)
Si_Cur_Cur = np.delete(Si_Cur_Cur,0)
Si_Cur_Cur = np.delete(Si_Cur_Cur,0)
Si_Cur_Cur = np.delete(Si_Cur_Cur,0)
Si_Cur_Vol = np.delete(Si_Cur_Vol,0)
Si_Cur_Vol = np.delete(Si_Cur_Vol,0)
Si_Cur_Vol = np.delete(Si_Cur_Vol,0)
Si_Cur_Vol = np.delete(Si_Cur_Vol,0)
factor_Si_Cur = np.delete(factor_Si_Cur,0)
factor_Si_Cur = np.delete(factor_Si_Cur,0)
factor_Si_Cur = np.delete(factor_Si_Cur,0)
factor_Si_Cur = np.delete(factor_Si_Cur,0)
print(f'Si resistance : {np.mean(factor_Si_Cur*Si_Cur_Vol*Si_thick/Si_Cur_Cur)}')

print(f'alpha = {np.mean(alpha)}')


plt.plot(Si_Tem_Tem, Si_Tem_Vol, color = 'blue', linewidth=7, linestyle = '-',label=fr'$P \, type \, Si \,\, (B \, high \, doped)$')
plt.plot(fitx, func1(fitx, a, b), color='cyan', linewidth=7,zorder=1, linestyle = ':', label=fr'$Fitting\,Si$')

plt.xlabel(r'$Temperature \,\, (^\circ C)$',size=50)
plt.ylabel(r'$Voltage \,\, [mV]$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
plt.xlim(30,120)
# plt.ylim(18,20.5)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, loc = 'upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Tem_vs_Vol.png')

fig.clear()

ax = plt.axes()
Si_Cur_Cur = np.delete(Si_Cur_Cur,0)
Si_Cur_Cur = np.delete(Si_Cur_Cur,0)
Si_Cur_Cur = np.delete(Si_Cur_Cur,0)
Si_Cur_Cur = np.delete(Si_Cur_Cur,0)
Si_Cur_Vol = np.delete(Si_Cur_Vol,0)
Si_Cur_Vol = np.delete(Si_Cur_Vol,0)
Si_Cur_Vol = np.delete(Si_Cur_Vol,0)
Si_Cur_Vol = np.delete(Si_Cur_Vol,0)
factor_Si_Cur = np.delete(factor_Si_Cur,0)
factor_Si_Cur = np.delete(factor_Si_Cur,0)
factor_Si_Cur = np.delete(factor_Si_Cur,0)
factor_Si_Cur = np.delete(factor_Si_Cur,0)

Al_Cur_Cur = np.delete(Al_Cur_Cur,0)
Al_Cur_Cur = np.delete(Al_Cur_Cur,0)
Al_Cur_Cur = np.delete(Al_Cur_Cur,0)
Al_Cur_Cur = np.delete(Al_Cur_Cur,0)
Al_Cur_Vol = np.delete(Al_Cur_Vol,0)
Al_Cur_Vol = np.delete(Al_Cur_Vol,0)
Al_Cur_Vol = np.delete(Al_Cur_Vol,0)
Al_Cur_Vol = np.delete(Al_Cur_Vol,0)
factor_Al = np.delete(factor_Al,0)
factor_Al = np.delete(factor_Al,0)
factor_Al = np.delete(factor_Al,0)
factor_Al = np.delete(factor_Al,0)

W_Cur_Cur = np.delete(W_Cur_Cur,0)
W_Cur_Cur = np.delete(W_Cur_Cur,0)
W_Cur_Cur = np.delete(W_Cur_Cur,0)
W_Cur_Cur = np.delete(W_Cur_Cur,0)
W_Cur_Vol = np.delete(W_Cur_Vol,0)
W_Cur_Vol = np.delete(W_Cur_Vol,0)
W_Cur_Vol = np.delete(W_Cur_Vol,0)
W_Cur_Vol = np.delete(W_Cur_Vol,0)
factor_W = np.delete(factor_W,0)
factor_W = np.delete(factor_W,0)
factor_W = np.delete(factor_W,0)
factor_W = np.delete(factor_W,0)



# print(Si_Cur_Cur)

plt.plot(Si_Cur_Cur, factor_Si_Cur*Si_Cur_Vol*Si_thick/Si_Cur_Cur, color = 'blue', linewidth=7, linestyle = '-',label=fr'$P \, type \, Si \,\, (B \, high \, doped)$')
# plt.plot(Al_Cur_Cur, factor_Al*Al_Cur_Vol*Al_thick/Al_Cur_Cur, color = 'red', linewidth=7, linestyle = '-',label=fr'$Al \, Foil$')
# plt.plot(W_Cur_Cur, factor_W*W_Cur_Vol*W_thick/W_Cur_Cur, color = 'green', linewidth=7, linestyle = '-',label=fr'$W \, Foil$')

# plt.plot(fitx, func1(fitx, a, b), color='cyan', linewidth=7,zorder=1, linestyle = ':', label=fr'$Fitting\,Si$')

# print(factor_Si_Cur*Si_Cur_Vol*Si_thick/Si_Cur_Cur)
# print(factor_Al*Al_Cur_Vol*Al_thick/Al_Cur_Cur)
# print(factor_W*W_Cur_Vol*W_thick/W_Cur_Cur)

print('\n')
print(f'Si resistance : {np.mean(factor_Si_Cur*Si_Cur_Vol*Si_thick/Si_Cur_Cur)}')
print(f'Al resistance : {np.mean(factor_Al*Al_Cur_Vol*Al_thick/Al_Cur_Cur)}')
print(f'W resistance : {np.mean(factor_W*W_Cur_Vol*W_thick/W_Cur_Cur)}')

plt.xlabel(r'$Current \,\, (mA)$',size=50)
plt.ylabel(r'$Resistance \,\, [\Omega m]$',size=50)

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])
# plt.xlim(25,200)
# plt.ylim(18,20.5)

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

# plt.ytic

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, loc = 'upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Cur_Resistance.png')
