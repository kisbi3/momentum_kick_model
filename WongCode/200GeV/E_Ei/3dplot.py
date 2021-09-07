import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
# from mpl_toolkits import mplot3d
from matplotlib import cm

mpl.rcParams["text.usetex"] = True

fig = plt.figure()
fig.set_size_inches(15, 10, forward=True)
ax3d = plt.axes(projection='3d')


detaf=np.loadtxt('E_Ei_3d.csv',delimiter=',',usecols=[0],skiprows=1)
dphif=np.loadtxt('E_Ei_3d.csv',delimiter=',',usecols=[1],skiprows=1)
E_Ei=np.loadtxt('E_Ei_3d.csv',delimiter=',',usecols=[2],skiprows=1)

ax3d = plt.axes(projection="3d")
# ax3d.plot_trisurf(detaf,dphif,E_Ei,cmap='jet', linewidths=0.5)
ax3d.plot_trisurf(dphif,detaf,E_Ei,cmap='jet', linewidths=0.5)

# plt.title(r'$\frac{E}{Ei}=\frac{\sqrt{p_{Tf}^2cosh(\eta_f)^2+m^2}}{\sqrt{p_{Ti}^2cosh(\eta_f)^2+m^2}}$', size = 40)
plt.title(r'$\frac{E}{Ei}=\frac{\sqrt{m_{\pi}^2+p_{Tf}^2}cosh(y_f)}{\sqrt{m_{\pi}^2+p_{Ti}^2}cosh(y_i)}$', size = 40)
ax3d.set_ylabel(r"$\Delta\eta$", size=25, labelpad = 25)
ax3d.set_xlabel(r"$\Delta\phi$", size=25, labelpad = 25)
# ax3d.set_zlabel(r"$\frac{E}{Ei}$", size=25, labelpad = 25)
# axes.set_zlabel("$\Delta\eta$")
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

plt.tight_layout()
# plt.show()

fig.savefig('E_Ei_3dplot_origin.png')

fig.clear()

fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

ptf3=np.loadtxt('rapidity_eta=3.csv',delimiter=',',usecols=[0],skiprows=1)
yi=np.loadtxt('rapidity_eta=3.csv',delimiter=',',usecols=[1],skiprows=1)
yf=np.loadtxt('rapidity_eta=3.csv',delimiter=',',usecols=[2],skiprows=1)
cosh3=np.loadtxt('rapidity_eta=3.csv',delimiter=',',usecols=[3],skiprows=1)

ptf0=np.loadtxt('rapidity_eta=0.csv',delimiter=',',usecols=[0],skiprows=1)
cosh0=np.loadtxt('rapidity_eta=0.csv',delimiter=',',usecols=[3],skiprows=1)

ptf15=np.loadtxt('rapidity_eta=1.5.csv',delimiter=',',usecols=[0],skiprows=1)
cosh15=np.loadtxt('rapidity_eta=1.5.csv',delimiter=',',usecols=[3],skiprows=1)

plt.plot(ptf3, cosh3, color = 'red', linewidth=7, linestyle = '-',label=r'$\eta = 3$')
plt.plot(ptf0, cosh0, color = 'black', linewidth=7, linestyle = '-',label=r'$\eta = 0$')
plt.plot(ptf15, cosh15, color = 'blue', linewidth=7, linestyle = '-',label=r'$\eta = 1.5$')

plt.xlabel(r'$p_{Tf}$',size=50)
# plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

plt.title(r'$\frac{cosh(y_f)}{cosh(y_i)}, \phi=0$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('coshy.png')