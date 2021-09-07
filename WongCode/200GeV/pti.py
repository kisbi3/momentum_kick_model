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

fig.set_size_inches(20, 16, forward=True)

pt=np.loadtxt('pti,q=1.csv',delimiter=',',usecols=[0],skiprows=1)
pti=np.loadtxt('pti,q=1.csv',delimiter=',',usecols=[1],skiprows=1)
rapid=np.loadtxt('pti,q=1.csv',delimiter=',',usecols=[2],skiprows=1)


plt.xlabel(r'$p_{tf}$',size=50)
plt.ylabel(r'$p_{ti}$',size=50)

plt.plot(pt,pti, color="green",linewidth=7,linestyle = '-',label=r'$p_{ti}$')
plt.plot(pt,rapid, color="red",linewidth=7,linestyle = '-',label=r'$rapidity$')

plt.minorticks_on()
# plt.yscale('log')
ax.axis([0,4,-15,10])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
plt.tight_layout()

fig.savefig('pti,q=1.png')

fig.clear()

ptf3d_q1=np.loadtxt('pti_3d,q=1.csv',delimiter=',',usecols=[0],skiprows=1)
phif3d_q1=np.loadtxt('pti_3d,q=1.csv',delimiter=',',usecols=[1],skiprows=1)
pti3d_q1=np.loadtxt('pti_3d,q=1.csv',delimiter=',',usecols=[2],skiprows=1)
pti3d_q1_ratio=np.loadtxt('pti_3d,q=1.csv',delimiter=',',usecols=[3],skiprows=1)



ax3d = plt.axes(projection="3d")
surf = ax3d.plot_trisurf(phif3d_q1, ptf3d_q1, pti3d_q1, cmap='jet', linewidths=0.5)

plt.title(r'$p_{ti}, q=1$', size = 40)
cbar = fig.colorbar(surf, shrink=0.7, aspect=4)
cbar.ax.tick_params(width=2,length=20, pad = 10, labelsize=35)
ax3d.set_xlabel(r"$\Delta\phi$", size=40, labelpad = 25)
ax3d.set_ylabel(r"$p_{tf}$", size=40, labelpad = 25)
ax3d.set_zlabel(r"$p_{ti}$", size=40, labelpad = 25)
# axes.set_zlabel("$\Delta\eta$")
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25)


plt.tight_layout()
# plt.show()

fig.savefig('pti_3d,q=1.png')

fig.clear()

ax3d = plt.axes(projection="3d")
surf = ax3d.plot_trisurf(phif3d_q1, ptf3d_q1, pti3d_q1_ratio, cmap='jet', linewidths=0.5)

cbar = fig.colorbar(surf, shrink=0.7, aspect=4)
cbar.ax.tick_params(width=2,length=20, pad = 10, labelsize=35)
plt.title(r'$p_{ti}, q=1\quad(ratio)$', size = 40)
ax3d.set_xlabel(r"$\Delta\phi$", size=40, labelpad = 25)
ax3d.set_ylabel(r"$p_{tf}$", size=40, labelpad = 25)
ax3d.set_zlabel(r"$p_{ti}/p_{tf}$", size=40, labelpad = 25)
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25)

plt.tight_layout()
# plt.show()

fig.savefig('pti_3d,q=1_ratio.png')

fig.clear()


ptf3d_q2=np.loadtxt('pti_3d,q=2.csv',delimiter=',',usecols=[0],skiprows=1)
phif3d_q2=np.loadtxt('pti_3d,q=2.csv',delimiter=',',usecols=[1],skiprows=1)
pti3d_q2=np.loadtxt('pti_3d,q=2.csv',delimiter=',',usecols=[2],skiprows=1)
pti3d_q2_ratio=np.loadtxt('pti_3d,q=2.csv',delimiter=',',usecols=[3],skiprows=1)


ax3d = plt.axes(projection="3d")
surf = ax3d.plot_trisurf(phif3d_q2, ptf3d_q2, pti3d_q2, cmap='jet', linewidths=0.5)

cbar = fig.colorbar(surf, shrink=0.7, aspect=4)
cbar.ax.tick_params(width=2,length=20, pad = 10, labelsize=35)
plt.title(r'$p_{ti}, q=2$', size = 40)
ax3d.set_xlabel(r"$\Delta\phi$", size=40, labelpad = 25)
ax3d.set_ylabel(r"$p_{tf}$", size=40, labelpad = 25)
ax3d.set_zlabel(r"$p_{ti}$", size=40, labelpad = 25)
# axes.set_zlabel("$\Delta\eta$")
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25)

plt.tight_layout()
# plt.show()

fig.savefig('pti_3d,q=2.png')

fig.clear()

ax3d = plt.axes(projection="3d")
surf = ax3d.plot_trisurf(phif3d_q2, ptf3d_q2, pti3d_q2_ratio, cmap='jet', linewidths=0.5)

cbar = fig.colorbar(surf, shrink=0.7, aspect=4)
cbar.ax.tick_params(width=2,length=20, pad = 10, labelsize=35)
plt.title(r'$p_{ti}, q=2\quad(ratio)$', size = 40)
ax3d.set_xlabel(r"$\Delta\phi$", size=40, labelpad = 25)
ax3d.set_ylabel(r"$p_{tf}$", size=40, labelpad = 25)
ax3d.set_zlabel(r"$p_{ti}/p_{tf}$", size=40, labelpad = 25)
# axes.set_zlabel("$\Delta\eta$")
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25)

plt.tight_layout()
# plt.show()

fig.savefig('pti_3d,q=2_ratio.png')


fig.clear()

rapid_etaf=np.loadtxt('rapidity_3d.csv',delimiter=',',usecols=[0],skiprows=1)
rapid_ptf=np.loadtxt('rapidity_3d.csv',delimiter=',',usecols=[1],skiprows=1)
rapidity=np.loadtxt('rapidity_3d.csv',delimiter=',',usecols=[2],skiprows=1)


ax3d = plt.axes(projection="3d")
# surf = ax3d.plot_trisurf(rapid_ptf, rapid_etaf, rapidity, cmap='jet', linewidths=0.5)
surf = ax3d.scatter(rapid_ptf, rapid_etaf, rapidity, cmap='jet', linewidths=0.5)

plt.title(r'$rapidity$', size = 40)
# cbar = fig.colorbar(surf, shrink=0.7, aspect=4)
# cbar.ax.tick_params(width=2,length=20, pad = 10, labelsize=35)
ax3d.set_xlabel(r"$p_{tf}$", size=40, labelpad = 25)
ax3d.set_ylabel(r"$\Delta\eta$", size=40, labelpad = 25)
ax3d.set_zlabel(r"$rapidity$", size=40, labelpad = 25)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25)
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25)


plt.tight_layout()
# plt.show()

fig.savefig('rapidity_3d.png')