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


detaf=np.loadtxt('3d_etapt_Ridge.csv',delimiter=',',usecols=[0],skiprows=1)
dptf=np.loadtxt('3d_etapt_Ridge.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge=np.loadtxt('3d_etapt_Ridge.csv',delimiter=',',usecols=[2],skiprows=1)

detaf_Jet=np.loadtxt('3d_etapt_Ridge_Jet.csv',delimiter=',',usecols=[0],skiprows=1)
dptf_Jet=np.loadtxt('3d_etapt_Ridge_Jet.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge_Jet=np.loadtxt('3d_etapt_Ridge_Jet.csv',delimiter=',',usecols=[2],skiprows=1)

detaf_E=np.loadtxt('3d_etapt_E.csv',delimiter=',',usecols=[0],skiprows=1)
dptf_E=np.loadtxt('3d_etapt_E.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge_E=np.loadtxt('3d_etapt_E.csv',delimiter=',',usecols=[2],skiprows=1)

detaf_remove_E=np.loadtxt('3d_etapt_Ridge_remove_E.csv',delimiter=',',usecols=[0],skiprows=1)
dptf_remove_E=np.loadtxt('3d_etapt_Ridge_remove_E.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge_remove_E=np.loadtxt('3d_etapt_Ridge_remove_E.csv',delimiter=',',usecols=[2],skiprows=1)


ax3d = plt.axes(projection="3d")
ax3d.plot_trisurf(detaf,dptf,Ridge,cmap='jet', linewidths=0.5)

plt.title(r'$Ridge$', size = 40)
ax3d.set_xlabel(r"$\Delta\eta$", size=25, labelpad = 25)
ax3d.set_ylabel(r"$P_T$", size=25, labelpad = 25)
# axes.set_zlabel("$\Delta\eta$")
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

plt.tight_layout()
# plt.show()

fig.savefig('Ridge_etapt_3dplot.png')



fig.clear()

ax3d = plt.axes(projection="3d")
ax3d.plot_trisurf(detaf_Jet,dptf_Jet,Ridge_Jet,cmap='jet', linewidths=0.5)

plt.title(r'$Ridge+Jet$', size = 40)
ax3d.set_xlabel(r"$\Delta\eta$", size=25, labelpad = 25)
ax3d.set_ylabel(r"$P_T$", size=25, labelpad = 25)
# axes.set_zlabel("$\Delta\eta$")
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

plt.tight_layout()
# plt.show()

fig.savefig('Ridge_Jet_etapt_3dplot.png')

fig.clear()

ax3d = plt.axes(projection="3d")
ax3d.plot_trisurf(detaf_E,dptf_E,Ridge_E,cmap='jet', linewidths=0.5)

plt.title(r'$E/E_i$', size = 40)
ax3d.set_xlabel(r"$\Delta\eta$", size=25, labelpad = 25)
ax3d.set_ylabel(r"$P_T$", size=25, labelpad = 25)
# axes.set_zlabel("$\Delta\eta$")
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

plt.tight_layout()
# plt.show()

fig.savefig('E_etapt_3dplot.png')

fig.clear()

ax3d = plt.axes(projection="3d")
ax3d.plot_trisurf(detaf_remove_E,dptf_remove_E,Ridge_remove_E,cmap='jet', linewidths=0.5)

plt.title(r'$Ridge, Remove\quad(E/E_i)$', size = 40)
ax3d.set_xlabel(r"$\Delta\eta$", size=25, labelpad = 25)
ax3d.set_ylabel(r"$P_T$", size=25, labelpad = 25)
# axes.set_zlabel("$\Delta\eta$")
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

plt.tight_layout()
# plt.show()

fig.savefig('Ridge_remove_E_etapt_3dplot.png')

fig.clear()



fig = plt.figure(figsize = (25,18))
ax = plt.axes()

fig.set_size_inches(23.384, 18, forward=True)

detaf_E_2d=np.loadtxt('eta_E.csv',delimiter=',',usecols=[0],skiprows=1)
Ridge_E_2d=np.loadtxt('eta_E.csv',delimiter=',',usecols=[1],skiprows=1)


plt.plot(detaf_E_2d, Ridge_E_2d, linewidth=3)

plt.minorticks_on()
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')

plt.grid(color='silver',linestyle=':',linewidth=3)
plt.tight_layout()

plt.title(r'$E/E_i$', size = 40)

plt.xlabel(r'$\eta$, size = 50')

fig.savefig('2d_E.png')