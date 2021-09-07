import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
# from mpl_toolkits import mplot3d
from matplotlib import cm

# mpl.rcParams["text.usetex"] = True

fig = plt.figure()
fig.set_size_inches(15, 10, forward=True)
ax3d = plt.axes(projection='3d')


detaf=np.loadtxt('3dplot_Ridge.csv',delimiter=',',usecols=[0],skiprows=1)
dphif=np.loadtxt('3dplot_Ridge.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge=np.loadtxt('3dplot_Ridge.csv',delimiter=',',usecols=[2],skiprows=1)

detaf_Jet=np.loadtxt('3dplot_Ridge_Jet.csv',delimiter=',',usecols=[0],skiprows=1)
dphif_Jet=np.loadtxt('3dplot_Ridge_Jet.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge_Jet=np.loadtxt('3dplot_Ridge_Jet.csv',delimiter=',',usecols=[2],skiprows=1)

detaf_E=np.loadtxt('3dplot_E.csv',delimiter=',',usecols=[0],skiprows=1)
dphif_E=np.loadtxt('3dplot_E.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge_E=np.loadtxt('3dplot_E.csv',delimiter=',',usecols=[2],skiprows=1)

detaf_remove_E=np.loadtxt('3dplot_Ridge_remove_E.csv',delimiter=',',usecols=[0],skiprows=1)
dphif_remove_E=np.loadtxt('3dplot_Ridge_remove_E.csv',delimiter=',',usecols=[1],skiprows=1)
Ridge_remove_E=np.loadtxt('3dplot_Ridge_remove_E.csv',delimiter=',',usecols=[2],skiprows=1)



ax3d = plt.axes(projection="3d")
ax3d.plot_trisurf(detaf,dphif,Ridge,cmap='jet', linewidths=0.5)

plt.title(r'$Ridge$', size = 40)
ax3d.set_xlabel(r"$\Delta\eta$", size=25, labelpad = 25)
ax3d.set_ylabel(r"$\Delta\phi$", size=25, labelpad = 25)
ax3d.set_zlabel(r"$dN/d\Delta \eta d\Delta\phi$", size=25, labelpad = 25)
# axes.set_zlabel("$\Delta\eta$")
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

plt.tight_layout()
# plt.show()

fig.savefig('Ridge_3dplot.png')



fig.clear()

ax3d = plt.axes(projection="3d")
ax3d.plot_trisurf(detaf_Jet,dphif_Jet,Ridge_Jet,cmap='jet', linewidths=0.5)

plt.title(r'$Ridge+Jet$', size = 40)
ax3d.set_xlabel(r"$\Delta\eta$", size=25, labelpad = 25)
ax3d.set_ylabel(r"$\Delta\phi$", size=25, labelpad = 25)
ax3d.set_zlabel(r"$dN/d\Delta \eta d\Delta\phi$", size=25, labelpad = 25)
# axes.set_zlabel("$\Delta\eta$")
plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

plt.tight_layout()
# plt.show()

fig.savefig('Ridge_Jet_3dplot.png')

# fig.clear()

# ax3d = plt.axes(projection="3d")
# ax3d.plot_trisurf(detaf_E,dphif_E,Ridge_E,cmap='jet', linewidths=0.5)

# plt.title(r'$E/E_i$', size = 40)
# ax3d.set_xlabel(r"$\Delta\eta$", size=25, labelpad = 25)
# ax3d.set_ylabel(r"$\Delta\phi$", size=25, labelpad = 25)
# # axes.set_zlabel("$\Delta\eta$")
# plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
# plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

# plt.tight_layout()
# # plt.show()

# fig.savefig('E_3dplot.png')

# fig.clear()

# ax3d = plt.axes(projection="3d")
# ax3d.plot_trisurf(detaf_remove_E,dphif_remove_E,Ridge_remove_E,cmap='jet', linewidths=0.5)

# plt.title(r'$Ridge, Remove\quad(E/E_i)$', size = 40)
# ax3d.set_xlabel(r"$\Delta\eta$", size=25, labelpad = 25)
# ax3d.set_ylabel(r"$\Delta\phi$", size=25, labelpad = 25)
# # axes.set_zlabel("$\Delta\eta$")
# plt.tick_params(axis='both',which='major',direction='in',width=2,length=30, pad = 10, labelsize=25, top = 'true', right = 'true')
# plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15, pad = 10, labelsize=25, top = 'true', right = 'true')

# plt.tight_layout()
# # plt.show()

# fig.savefig('Ridge_remove_E_3dplot.png')