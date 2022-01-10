import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit
import csv

periph_276_phi=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[0],skiprows=1)
periph_276_data=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[1],skiprows=1)
periph_130_phi=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[2],skiprows=1)
periph_130_data=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[3],skiprows=1)
phi_13TeV_90_up=np.loadtxt('./atlasgraphs/13TeV_90~.csv',delimiter=',',usecols=[0])
data_13TeV_90_up=np.loadtxt('./atlasgraphs/13TeV_90~.csv',delimiter=',',usecols=[1])
phi_13TeV_130_up=np.loadtxt('./atlasgraphs/13TeV_130~.csv',delimiter=',',usecols=[0])
data_13TeV_130_up=np.loadtxt('./atlasgraphs/13TeV_130~.csv',delimiter=',',usecols=[1])
phi_13TeV_120_130=np.loadtxt('./atlasgraphs/13TeV_120~130.csv',delimiter=',',usecols=[0])
data_13TeV_120_130=np.loadtxt('./atlasgraphs/13TeV_120~130.csv',delimiter=',',usecols=[1])
phi_13TeV_110_120=np.loadtxt('./atlasgraphs/13TeV_110~120.csv',delimiter=',',usecols=[0])
data_13TeV_110_120=np.loadtxt('./atlasgraphs/13TeV_110~120.csv',delimiter=',',usecols=[1])
phi_13TeV_100_110=np.loadtxt('./atlasgraphs/13TeV_100~110.csv',delimiter=',',usecols=[0])
data_13TeV_100_110=np.loadtxt('./atlasgraphs/13TeV_100~110.csv',delimiter=',',usecols=[1])
phi_13TeV_90_100=np.loadtxt('./atlasgraphs/13TeV_90~100.csv',delimiter=',',usecols=[0])
data_13TeV_90_100=np.loadtxt('./atlasgraphs/13TeV_90~100.csv',delimiter=',',usecols=[1])



mpl.rcParams["text.usetex"] = True

fig = plt.figure()
ax = plt.axes()
fig.set_size_inches(25, 20, forward=True)

def periph(phi, a, b, c, d):
    return a*np.cos(phi)+b*np.cos(2*phi)+c*np.cos(3*phi)+d

def atlas_ridge(phi, G, v):          #atlas는 0.5<pt<5 이므로 지금까지 그려온 그래프와 맞지 않음.
        return G*(1+2*v*np.cos(2*phi))

popt_periph_130TeV, pcov_periph_130TeV = curve_fit(periph, periph_130_phi, periph_130_data)
popt_periph_276TeV, pcov_periph_276TeV = curve_fit(periph, periph_276_phi, periph_276_data)

G_result=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[2],skiprows=1)
atlas_phi = np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV_CMS_ALICE/phiCorrelation_pt0-5.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result = np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV_CMS_ALICE/phiCorrelation_pt0-5.csv',delimiter=',',usecols=[3],skiprows=1)


Deltaphi = np.arange(-1.28,1.28,0.04)
v22_result_multi=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[0],skiprows=1)
v22_result=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[3],skiprows=1)

# for i in range(4):
Multi_90_up = atlas_ridge(Deltaphi, G_result[0], v22_result[0])
Multi_130_up = atlas_ridge(Deltaphi, G_result[1], v22_result[1])
Multi_120_130 = atlas_ridge(Deltaphi, G_result[2], v22_result[2])
Multi_110_120 = atlas_ridge(Deltaphi, G_result[3], v22_result[3])
Multi_100_110 = atlas_ridge(Deltaphi, G_result[4], v22_result[4])
Multi_90_100 = atlas_ridge(Deltaphi, G_result[5], v22_result[5])


fig.set_size_inches(40, 20, forward=True)

#단순 theory 계산값 -Czyam
plt.plot(atlas_phi, atlas_result-min(atlas_result), color = 'blue', linewidth=5, linestyle = '-',label=r'$result$')

#atlas의 Yridge를 계산한 값 -Czyam
plt.plot(Deltaphi, Multi_90_up-min(Multi_90_up), color = 'red', linewidth=5, linestyle = '--',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(Deltaphi, Multi_130_up-min(Multi_130_up), color = 'black', linewidth=5, linestyle = '--',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(Deltaphi, Multi_120_130-min(Multi_120_130), color = 'green', linewidth=5, linestyle = '--',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(Deltaphi, Multi_110_120-min(Multi_110_120), color = 'orange', linewidth=5, linestyle = '--',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(Deltaphi, Multi_100_110-min(Multi_100_110), color = 'blue', linewidth=5, linestyle = '--',label=r'$100 \leq N_{ch}^{rec}<110$')
plt.plot(Deltaphi, Multi_90_100-min(Multi_90_100), color = 'grey', linewidth=5, linestyle = '--',label=r'$90 \leq N_{ch}^{rec}<100$')

plt.legend(fontsize=45,framealpha=False,loc='upper right')

#단순히 data-Czyam
plt.plot(phi_13TeV_90_up, data_13TeV_90_up-min(data_13TeV_90_up), marker = 'o', markersize = 15, color = 'red', linewidth=5, linestyle = ':',label=r'$90 \leq N_{ch}^{rec}$')
plt.plot(phi_13TeV_130_up, data_13TeV_130_up-min(data_13TeV_130_up), marker = 'o', markersize = 15, color = 'black', linewidth=5, linestyle = ':',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(phi_13TeV_120_130, data_13TeV_120_130-min(data_13TeV_120_130), marker = 'o', markersize = 15, color = 'green', linewidth=5, linestyle = ':',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(phi_13TeV_110_120, data_13TeV_110_120-min(data_13TeV_110_120), marker = 'o', markersize = 15, color = 'orange', linewidth=5, linestyle = ':',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(phi_13TeV_100_110, data_13TeV_100_110-min(data_13TeV_100_110), marker = 'o', markersize = 15, color = 'blue', linewidth=5, linestyle = ':',label=r'$100 \leq N_{ch}^{rec}<110$')
plt.plot(phi_13TeV_90_100, data_13TeV_90_100-min(data_13TeV_90_100), marker = 'o', markersize = 15, color = 'grey', linewidth=5, linestyle = ':',label=r'$90 \leq N_{ch}^{rec}<100$')



plt.ylabel(r'$Y^{ridge}-C_{ZYAM}$', size = 70)
plt.xlabel(r'$\Delta\phi$', size=70)

plt.minorticks_on()

plt.ylim(-0.01,0.1)
plt.xlim(-1.3,1.3)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.legend(fontsize=45,framealpha=False,bbox_to_anchor=(1.13,0.4), ncol=2)
# plt.legend(fontsize=45,framealpha=False,loc='upper right')


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/high_multi-czyam.png')

fig.clear()


plt.plot(atlas_phi, atlas_result, color = 'blue', linewidth=5, linestyle = '-',label=r'$result$')
plt.plot(Deltaphi, Multi_130_up, color = 'black', linewidth=5, linestyle = '--',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(Deltaphi, Multi_120_130, color = 'green', linewidth=5, linestyle = '--',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(Deltaphi, Multi_110_120, color = 'orange', linewidth=5, linestyle = '--',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(Deltaphi, Multi_100_110, color = 'blue', linewidth=5, linestyle = '--',label=r'$100 \leq N_{ch}^{rec}<110$')

plt.ylabel(r'$Y^{ridge}$', size = 70)
plt.xlabel(r'$\Delta\phi$', size=70)

plt.minorticks_on()

# plt.ylim(-0.01,0.1)
plt.xlim(-1.3,1.3)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.legend(fontsize=45,framealpha=False,bbox_to_anchor=(1.,0.5))
plt.legend(fontsize=45,framealpha=False,loc='upper right')


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/high_multi_Yridge.png')

fig.clear()

# theory - Czyam    //    Yridge - Czyam 비교

F_result_multi=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[0],skiprows=1)
F_result=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[1],skiprows=1)

#Theory Result
theory_90_up = atlas_result-F_result[0]*periph(atlas_phi, *popt_periph_130TeV)
theory_130_up = atlas_result-F_result[1]*periph(atlas_phi, *popt_periph_130TeV)
theory_120_130 = atlas_result-F_result[2]*periph(atlas_phi, *popt_periph_130TeV)
theory_110_120 = atlas_result-F_result[3]*periph(atlas_phi, *popt_periph_130TeV)
theory_100_110 = atlas_result-F_result[4]*periph(atlas_phi, *popt_periph_130TeV)
theory_90_100 = atlas_result-F_result[5]*periph(atlas_phi, *popt_periph_130TeV)

plt.plot(atlas_phi, theory_90_up-min(theory_90_up), color = 'red', linewidth=5, linestyle = '-',label=r'$90 \leq N_{ch}^{rec}$')
plt.plot(atlas_phi, theory_130_up-min(theory_130_up), color = 'black', linewidth=5, linestyle = '-',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(atlas_phi, theory_120_130-min(theory_120_130), color = 'green', linewidth=5, linestyle = '-',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(atlas_phi, theory_110_120-min(theory_110_120), color = 'orange', linewidth=5, linestyle = '-',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(atlas_phi, theory_100_110-min(theory_100_110), color = 'blue', linewidth=5, linestyle = '-',label=r'$100 \leq N_{ch}^{rec}<110$')
plt.plot(atlas_phi, theory_90_100-min(theory_90_100), color = 'grey', linewidth=5, linestyle = '-',label=r'$90 \leq N_{ch}^{rec}<100$')

plt.legend(fontsize=45,framealpha=False,loc='upper right')

#experiment result
plt.plot(Deltaphi, Multi_90_up-min(Multi_90_up), color = 'red', linewidth=5, linestyle = '--',label=r'$90 \leq N_{ch}^{rec}$')
plt.plot(Deltaphi, Multi_130_up-min(Multi_130_up), color = 'black', linewidth=5, linestyle = '--',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(Deltaphi, Multi_120_130-min(Multi_120_130), color = 'green', linewidth=5, linestyle = '--',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(Deltaphi, Multi_110_120-min(Multi_110_120), color = 'orange', linewidth=5, linestyle = '--',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(Deltaphi, Multi_100_110-min(Multi_100_110), color = 'blue', linewidth=5, linestyle = '--',label=r'$100 \leq N_{ch}^{rec}<110$')
plt.plot(Deltaphi, Multi_90_100-min(Multi_90_100), color = 'grey', linewidth=5, linestyle = '--',label=r'$90 \leq N_{ch}^{rec}<100$')


plt.ylabel(r'$Y^{ridge}-C_{ZYAM}$', size = 70)
plt.xlabel(r'$\Delta\phi$', size=70)

plt.minorticks_on()

# plt.ylim(-0.01,0.1)
plt.xlim(-1.3,1.3)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.legend(fontsize=45,framealpha=False,bbox_to_anchor=(1.,0.5))



plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/high_multi_Yridge_czyam.png')

fig.clear()

#C_zyam 제외

theory_90_up = atlas_result-F_result[0]*periph(atlas_phi, *popt_periph_130TeV)
theory_130_up = atlas_result-F_result[1]*periph(atlas_phi, *popt_periph_130TeV)
theory_120_130 = atlas_result-F_result[2]*periph(atlas_phi, *popt_periph_130TeV)
theory_110_120 = atlas_result-F_result[3]*periph(atlas_phi, *popt_periph_130TeV)
theory_100_110 = atlas_result-F_result[4]*periph(atlas_phi, *popt_periph_130TeV)
theory_90_100 = atlas_result-F_result[5]*periph(atlas_phi, *popt_periph_130TeV)

plt.plot(atlas_phi, theory_90_up, color = 'red', linewidth=5, linestyle = '-',label=r'$90 \leq N_{ch}^{rec}$')
plt.plot(atlas_phi, theory_130_up, color = 'black', linewidth=5, linestyle = '-',label=r'$130 \leq N_{ch}^{rec}$')
plt.plot(atlas_phi, theory_120_130, color = 'green', linewidth=5, linestyle = '-',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(atlas_phi, theory_110_120, color = 'orange', linewidth=5, linestyle = '-',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(atlas_phi, theory_100_110, color = 'blue', linewidth=5, linestyle = '-',label=r'$100 \leq N_{ch}^{rec}<110$')
plt.plot(atlas_phi, theory_90_100, color = 'grey', linewidth=5, linestyle = '-',label=r'$90 \leq N_{ch}^{rec}<100$')

plt.legend(fontsize=45,framealpha=False)

#experiment result
# plt.plot(Deltaphi, Multi_90_up color = 'red', linewidth=5, linestyle = '--',label=r'$90 \leq N_{ch}^{rec}$')
# plt.plot(Deltaphi, Multi_130_up, color = 'black', linewidth=5, linestyle = '--',label=r'$130 \leq N_{ch}^{rec}$')
# plt.plot(Deltaphi, Multi_120_130, color = 'green', linewidth=5, linestyle = '--',label=r'$120 \leq N_{ch}^{rec}<130$')
# plt.plot(Deltaphi, Multi_110_120, color = 'orange', linewidth=5, linestyle = '--',label=r'$110 \leq N_{ch}^{rec}<120$')
# plt.plot(Deltaphi, Multi_100_110, color = 'blue', linewidth=5, linestyle = '--',label=r'$100 \leq N_{ch}^{rec}<110$')
# plt.plot(Deltaphi, Multi_90_100, color = 'grey', linewidth=5, linestyle = '--',label=r'$90 \leq N_{ch}^{rec}<100$')


plt.ylabel(r'$Y^{ridge}$', size = 70)
plt.xlabel(r'$\Delta\phi$', size=70)

plt.minorticks_on()

# plt.ylim(-0.01,0.1)
plt.xlim(-1.3,1.3)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.legend(fontsize=45,framealpha=False,bbox_to_anchor=(1.,0.5))



plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/13TeV_Ytheory.png')

fig.clear()