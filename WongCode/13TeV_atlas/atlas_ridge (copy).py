import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit
import csv

# fig = plt.figure()
# ax = plt.axes()

# fig.set_size_inches(25, 20, forward=True)

# fig, axes = plt.subplots(nrows=2, ncols=5, sharex=True, sharey='row', figsize=(100, 30))
fig, axes = plt.subplots(nrows=1, ncols=5, sharex=True, figsize=(100, 30))
mpl.rcParams["text.usetex"] = True

param_F=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[1],skiprows=2)
param_G=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[2],skiprows=2)
param_v=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[3],skiprows=2)
periph_130_phi=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[2],skiprows=1)
periph_130_data=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[3],skiprows=1)
# print(param_F)

def ridge(phi, G, v):
    return G*(1+2*v*np.cos(2*phi))

near_side_phi = np.arange(-1., 1., 0.04)
totl_side_phi = np.arange(-1.5, 4.7, 0.04)

for i in range(1):
    for j in range(5):
        if i==1 and j==2:
            axes[i,j].tick_params(axis='both',which='major',direction='in',width=2, length=20, labelsize=45, top = 'true', right='true')
            axes[i,j].tick_params(axis='both',which='minor',direction='in',width=2, length=10, labelsize=45, top = 'true', right='true')
            axes[i,j].set_xlabel(r'$\Delta \phi$', size = 70)
            break
        # print(4-j)
        graph_ridge = ridge(totl_side_phi, param_G[4-j], param_v[4-j])
        templ_ridge = ridge(periph_130_phi, param_G[4-j], param_v[4-j])
        graph_perip = param_F[4-j]*periph_130_data
        axes[i,j].plot(totl_side_phi, graph_ridge, linewidth=10, color = 'orange', linestyle='--', label=r'$Y^{ridge}$')
        axes[i,j].plot(periph_130_phi, graph_perip+param_G[4-j], drawstyle='steps-mid', linewidth=10, color = 'blue', linestyle = '-',label=r'$FY^{periph}$')
        # axes[i+1,j].plot(periph_130_phi, graph_perip+param_G[4-j], drawstyle='steps-mid', linewidth=10, color = 'blue', linestyle = '-',label=r'$FY^{periph}$')
        axes[i,j].plot(periph_130_phi, templ_ridge+graph_perip, drawstyle='steps-mid', color = 'red', linewidth=10, linestyle = '-',label=r'$Y^{templ}(\Delta\phi)$')
        axes[i,j].set_xlim(-1.,1.)
        axes[i+1,j].set_xlim(-1.,1.)

        if j==0:
            axes[i,j].set_ylabel(r'$Y(\Delta \phi)$', size = 70)
            axes[i+1,j].set_ylabel(r'$Y(\Delta \phi)$', size = 70)
        if j==4:
            axes[i,j].set_title(r'$130 \leq N^{rec}_{ch}$', size = 70, pad=30)
        else:
            st = 10*i+90
            en = 10*i+100
            axes[i,j].set_title(str(st)+r'$\leq N^{rec}_{ch}<$'+str(en), size = 70, pad=30)
        # if i==1:
        #     axes[i,j].set_xlabel(r'$\Delta \phi$', size = 70)
        

        # plt.title(r'$90 \leq N^{rec}_{ch}<100$')

        axes[i,j].minorticks_on()
        axes[i,j].tick_params(axis='both',which='major',direction='in',width=2, length=20, labelsize=45, top = 'true', right='true')
        axes[i,j].tick_params(axis='both',which='minor',direction='in',width=2, length=10, labelsize=45, top = 'true', right='true')
        axes[i,j].legend(framealpha=False, fontsize = 60)
        axes[i,j].grid(color='silver',linestyle=':', linewidth=5)
        
        axes[i+1,j].minorticks_on()
        axes[i+1,j].tick_params(axis='both',which='major',direction='in',width=2, length=20, labelsize=45, top = 'true', right='true')
        axes[i+1,j].tick_params(axis='both',which='minor',direction='in',width=2, length=10, labelsize=45, top = 'true', right='true')
        axes[i+1,j].legend(framealpha=False, fontsize = 60)
        axes[i+1,j].grid(color='silver',linestyle=':', linewidth=5)

# plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout(h_pad=-1, w_pad=-1)

# plt.plot(totl_side_phi, ridge(totl_side_phi, param_G[1], param_v[1]), color = 'orange', linestyle = '--',label=r'$Y^{ridge}$')
# plt.plot(periph_130_phi, param_F[1]*periph_130_data, color = 'blue', linestyle = '--',label=r'$FY^{periph}$')



# # plt.xlim(-1.55,4.7)

# plt.ylabel(r'$Y(\Delta \phi)$')
# plt.xlabel(r'$\Delta \phi$')

# plt.title(r'$130 \leq N^{rec}_{ch}$')

# plt.minorticks_on()

# plt.tick_params(axis='both',which='major',direction='in', top = 'true')
#2plt.tick_params(axis='both',which='minor',direction='in', top = 'true')
#5lt.legend(framealpha=False)


# plt.grid(color='silver',linestyle=':')



# plt.plot(totl_side_phi, ridge(totl_side_phi, param_G[2], param_v[2]), color = 'orange', linestyle = '--',label=r'$Y^{ridge}$')
# plt.plot(periph_130_phi, param_F[2]*periph_130_data, color = 'blue', linestyle = '--',label=r'$FY^{periph}$')

# plt.ylabel(r'$Y(\Delta \phi)$')
# plt.xlabel(r'$\Delta \phi$')

# plt.title(r'$120 \leq N^{rec}_{ch}<130$')

# plt.minorticks_on()

# plt.tick_params(axis='both',which='major',direction='in', top = 'true')
#2plt.tick_params(axis='both',which='minor',direction='in', top = 'true')
#5lt.legend(framealpha=False)


# plt.grid(color='silver',linestyle=':')


# plt.plot(totl_side_phi, ridge(totl_side_phi, param_G[3], param_v[3]), color = 'orange', linestyle = '--',label=r'$Y^{ridge}$')
# plt.plot(periph_130_phi, param_F[3]*periph_130_data, color = 'blue', linestyle = '--',label=r'$FY^{periph}$')

# plt.ylabel(r'$Y(\Delta \phi)$')
# plt.xlabel(r'$\Delta \phi$')

# plt.title(r'$110 \leq N^{rec}_{ch}<120$')

# plt.minorticks_on()

# plt.tick_params(axis='both',which='major',direction='in', top = 'true')
#2plt.tick_params(axis='both',which='minor',direction='in', top = 'true')
#5lt.legend(framealpha=False)


# plt.grid(color='silver',linestyle=':')


# plt.plot(totl_side_phi, ridge(totl_side_phi, param_G[4], param_v[4]), color = 'orange', linestyle = '--',label=r'$Y^{ridge}$')
# plt.plot(periph_130_phi, param_F[4]*periph_130_data, color = 'blue', linestyle = '--',label=r'$FY^{periph}$')

# plt.ylabel(r'$Y(\Delta \phi)$')
# plt.xlabel(r'$\Delta \phi$')

# plt.title(r'$100 \leq N^{rec}_{ch}<110$')

# plt.minorticks_on()

# plt.tick_params(axis='both',which='major',direction='in', top = 'true')
#2plt.tick_params(axis='both',which='minor',direction='in', top = 'true')
#5lt.legend(framealpha=False)


# plt.grid(color='silver',linestyle=':')


# plt.plot(totl_side_phi, ridge(totl_side_phi, param_G[5], param_v[5]), color = 'orange', linestyle = '--',label=r'$Y^{ridge}$')
# plt.plot(periph_130_phi, param_F[5]*periph_130_data, color = 'blue', linestyle = '--',label=r'$FY^{periph}$')

# plt.ylabel(r'$Y(\Delta \phi)$')
# plt.xlabel(r'$\Delta \phi$')

# plt.title(r'$90 \leq N^{rec}_{ch}<100$')

# plt.minorticks_on()

# plt.tick_params(axis='both',which='major',direction='in', top = 'true')
#2plt.tick_params(axis='both',which='minor',direction='in', top = 'true')
#5lt.legend(framealpha=False)


# plt.grid(color='silver',linestyle=':')

# plt.tight_layout()

fig.savefig('./result/13TeV_my.png')

fig.clear()