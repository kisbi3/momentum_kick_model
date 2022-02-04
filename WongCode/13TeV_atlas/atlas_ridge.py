import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit
import csv

mpl.rcParams["text.usetex"] = True

# data, template, ridge, peirpheral
param_F=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[1],skiprows=2)
param_G=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[2],skiprows=2)
param_v=np.loadtxt('./result/ridge_parameters.csv',delimiter=',',usecols=[3],skiprows=2)
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
periph_130_phi=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[2],skiprows=1)
periph_130_data=np.loadtxt('./atlasgraphs/Low_multiplicity.csv',delimiter=',',usecols=[3],skiprows=1)
# print(param_F)

# data 계산값들
atlas_phi_095 = np.loadtxt('phiCorrelation_pt0-5_95.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result_095 = np.loadtxt('phiCorrelation_pt0-5_95.csv',delimiter=',',usecols=[3],skiprows=1)
atlas_phi_105 = np.loadtxt('phiCorrelation_pt0-5_105.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result_105 = np.loadtxt('phiCorrelation_pt0-5_105.csv',delimiter=',',usecols=[3],skiprows=1)
atlas_phi_115 = np.loadtxt('phiCorrelation_pt0-5_115.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result_115 = np.loadtxt('phiCorrelation_pt0-5_115.csv',delimiter=',',usecols=[3],skiprows=1)
atlas_phi_125 = np.loadtxt('phiCorrelation_pt0-5_125.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result_125 = np.loadtxt('phiCorrelation_pt0-5_125.csv',delimiter=',',usecols=[3],skiprows=1)
atlas_phi_135 = np.loadtxt('phiCorrelation_pt0-5_135.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_result_135 = np.loadtxt('phiCorrelation_pt0-5_135.csv',delimiter=',',usecols=[3],skiprows=1)


# fig, axes = plt.subplots(nrows=2, ncols=5, sharex=True, sharey='row', figsize=(100, 30))

# template, data, ridge, peripheral 한번에 그리기 위한 figure 생성
fig1, (axes1, axes2, axes3) = plt.subplots(nrows=3, ncols=5,figsize=(100,20))

# 이론값과 Yridge를 비교할 그림 생성
fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(25, 20, forward=True)
# 이론값과의 그래프 색을 일치 시키기 위한 색
color_ridge = ['grey', 'blue', 'orange', 'green', 'black']

# axes2 = axes.twinx()


def ridge(phi, G, v):
    return G*(1+2*v*np.cos(2*phi))

near_side_phi = np.arange(-1., 1.04, 0.04)
totl_side_phi = np.arange(-1.5, 4.7, 0.04)
integrate_phi = np.arange(-1., 1., 0.001)
# print(integrate_phi)
ridge_integrate = np.zeros(5)

for j in range(5):
    # print(periph_130_data)
    graph_ridge = ridge(near_side_phi, param_G[4-j], param_v[4-j])
    templ_ridge = ridge(periph_130_phi, param_G[4-j], param_v[4-j])
    graph_perip = param_F[4-j]*periph_130_data
    graph_templ = templ_ridge + graph_perip

    # -Czyam인 상태에서의 적분값을 비교해야 하므로 -Czyam 분량의 사각형을 빼주어야 한다.
    ridge_integrate[j] = np.trapz(ridge(integrate_phi, param_G[4-j], param_v[4-j])-min(ridge(integrate_phi, param_G[4-j], param_v[4-j])), x = integrate_phi)
    # ridge_integrate[j] -= 2.*min(ridge(integrate_phi, param_G[4-j], param_v[4-j]))

    # axes2 = axes[j].twinx()
    axes1[j].plot(near_side_phi, graph_ridge, linewidth=10, color = 'orange', linestyle='--', label=r'$Y^{ridge}$') 
    axes1[j].plot(periph_130_phi, graph_perip, linewidth=10, color = 'blue', linestyle = '-',label=r'$FY^{periph}$')
    axes1[j].plot(periph_130_phi, graph_templ, linewidth=10, color = 'red', linestyle = '-',label=r'$Y^{templ}$')
    axes2[j].plot(near_side_phi, graph_ridge, linewidth=10, color = 'orange', linestyle='--', label=r'$Y^{ridge}$') 
    ax.plot(near_side_phi, graph_ridge-min(graph_ridge), linewidth=5, color = color_ridge[j], linestyle='--') 
    # axes2[j].plot(periph_130_phi, graph_perip, linewidth=10, color = 'blue', linestyle = '-',label=r'$FY^{periph}$')
    # axes2[j].plot(periph_130_phi, graph_templ, linewidth=10, color = 'red', linestyle = '-',label=r'$Y^{templ}$')
    # axes3[j].plot(near_side_phi, graph_ridge, linewidth=10, color = 'orange', linestyle='--', label=r'$Y^{ridge}$') 
    axes3[j].plot(periph_130_phi, graph_perip, linewidth=10, color = 'blue', linestyle = '-',label=r'$FY^{periph}$')
    # axes3[j].plot(periph_130_phi, graph_templ, linewidth=10, color = 'red', linestyle = '-',label=r'$Y^{templ}$')


    # axes1[j].set_ylim(min(graph_templ)-0.01, min(graph_templ)+0.07)
    # print(graph_perip)
    axes2[j].set_ylim(min(graph_ridge)-0.01, min(graph_ridge)+0.07)
    axes3[j].set_ylim(min(graph_perip)-0.01, min(graph_perip)+0.07)
    axes1[j].set_xlim(-1.,1.)
    axes2[j].set_xlim(-1.,1.)
    axes3[j].set_xlim(-1.,1.)
    axes3[j].set_xlabel(r'$\Delta \phi$', size = 70)
    # axes[j].set_yscale('log')
    # if j==0:
    #     axes1[j].set_ylabelr'$Y^{data}$', size = 80)
        
    if j==4:
        axes1[j].set_title(r'$130 \leq N^{rec}_{ch}$', size = 70, pad=30)
        axes1[j].scatter(phi_13TeV_130_up, data_13TeV_130_up, color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
        axes1[j].set_ylim(min(data_13TeV_130_up)-0.01, min(data_13TeV_130_up)+0.07)
        axes1[j].legend(framealpha=False, fontsize = 70, bbox_to_anchor=(1.05,1.0), loc='upper left')
        
    else:
        st = 10*j+90
        en = 10*j+100
        axes1[j].set_title(str(st)+r'$\leq N^{rec}_{ch}<$'+str(en), size = 70, pad=30)
    if j == 3:
        axes1[j].scatter(phi_13TeV_120_130, data_13TeV_120_130, color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
        axes1[j].set_ylim(min(data_13TeV_120_130)-0.01, min(data_13TeV_120_130)+0.07)
    elif j == 2:
        axes1[j].scatter(phi_13TeV_110_120, data_13TeV_110_120, color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
        axes1[j].set_ylim(min(data_13TeV_110_120)-0.01, min(data_13TeV_110_120)+0.07)
    elif j == 1:
        axes1[j].scatter(phi_13TeV_100_110, data_13TeV_100_110, color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
        axes1[j].set_ylim(min(data_13TeV_100_110)-0.01, min(data_13TeV_100_110)+0.07)
    elif j == 0:
        axes1[j].scatter(phi_13TeV_90_100, data_13TeV_90_100, color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
        axes1[j].set_ylim(min(data_13TeV_90_100)-0.01, min(data_13TeV_90_100)+0.07)
        axes2[j].set_ylabel(r'$Y(\Delta\phi)$', size = 80)


    # if i==1:
    #     axes[j].set_xlabel(r'$\Delta \phi$', size = 70)
    

    # plt.title(r'$90 \leq N^{rec}_{ch}<100$')

    axes1[j].minorticks_on()
    axes2[j].minorticks_on()
    axes3[j].minorticks_on()
    
    axes1[j].tick_params(axis='both',which='major',direction='in',width=2, length=20, labelsize=45, right = 'true')
    axes1[j].tick_params(axis='both',which='minor',direction='in',width=2, length=10, labelsize=45, right = 'true')

    axes2[j].tick_params(axis='y',which='major',direction='in',width=2, length=20, labelsize=45, right = 'true')
    axes2[j].tick_params(axis='y',which='minor',direction='in',width=2, length=10, labelsize=45, right = 'true')
    axes2[j].tick_params(axis='x',which='major',direction='in',width=2, length=0, labelsize=45, right = 'true')
    axes2[j].tick_params(axis='x',which='minor',direction='in',width=2, length=0, labelsize=45, right = 'true')
    
    axes3[j].tick_params(axis='both',which='major',direction='in',width=2, length=20, labelsize=45, right = 'true')
    axes3[j].tick_params(axis='both',which='minor',direction='in',width=2, length=10, labelsize=45, right = 'true')
    
    axes1[j].grid(color='silver',linestyle=':', linewidth=5)
    axes2[j].grid(color='silver',linestyle=':', linewidth=5)
    axes3[j].grid(color='silver',linestyle=':', linewidth=5)

    axes1[j].spines.bottom.set_visible(False)
    axes2[j].spines.top.set_visible(False)
    axes2[j].spines.bottom.set_visible(False)
    axes3[j].spines.top.set_visible(False)

    axes1[j].xaxis.tick_top()
    axes1[j].tick_params(labeltop=False)  # don't put tick labels at the top
    # axes2[j].get_xaxis().set_visible(False)
    axes2[j].xaxis.tick_top()
    axes2[j].tick_params(labeltop = False)
    axes3[j].xaxis.tick_bottom()

    d = .5
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=40, linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    axes1[j].plot([0, 1], [0, 0], transform=axes1[j].transAxes, **kwargs)
    axes2[j].plot([0, 1], [0, 0], [0, 1], [1, 1], transform=axes2[j].transAxes, **kwargs)
    # axes2[j].plot([0, 0], [1, 0], transform=axes2[j].transAxes, **kwargs)
    axes3[j].plot([0, 1], [1, 1], transform=axes3[j].transAxes, **kwargs)
    
fig1.tight_layout(h_pad = -12)

fig1.savefig('./result/13TeV_my.png')
# , transparent = True

fig1.clear()


print(ridge_integrate)

theory_norm = np.loadtxt('phiCorr_norm.csv', delimiter=',', usecols=[1])

#그냥 numpy로 적분하자.
theor_integrate = np.zeros(5)
theor_integrate[0] = np.trapz(atlas_result_095-min(atlas_result_095), x = atlas_phi_095)
theor_integrate[1] = np.trapz(atlas_result_105-min(atlas_result_105), x = atlas_phi_105)
theor_integrate[2] = np.trapz(atlas_result_115-min(atlas_result_115), x = atlas_phi_115)
theor_integrate[3] = np.trapz(atlas_result_125-min(atlas_result_125), x = atlas_phi_125)
theor_integrate[4] = np.trapz(atlas_result_135-min(atlas_result_135), x = atlas_phi_135)

# 이론값에서도 -Czyam인 상태에서의 적분값이 필요하므로 이에 해당하는 직사각형을 빼준다.
# theory_norm[0] -= 2*min(atlas_result_135)
# theory_norm[1] -= 2*min(atlas_result_125)
# theory_norm[2] -= 2*min(atlas_result_115)
# theory_norm[3] -= 2*min(atlas_result_105)
# theory_norm[4] -= 2*min(atlas_result_095)

# ratio_095 = ridge_integrate[0]/theory_norm[4]
# ratio_105 = ridge_integrate[1]/theory_norm[3]
# ratio_115 = ridge_integrate[2]/theory_norm[2]
# ratio_125 = ridge_integrate[3]/theory_norm[1]
# ratio_135 = ridge_integrate[4]/theory_norm[0]

# theory_norm_095 = (atlas_result_095-min(atlas_result_095))*ratio_095
# theory_norm_105 = (atlas_result_105-min(atlas_result_105))*ratio_105
# theory_norm_115 = (atlas_result_115-min(atlas_result_115))*ratio_115
# theory_norm_125 = (atlas_result_125-min(atlas_result_125))*ratio_125
# theory_norm_135 = (atlas_result_135-min(atlas_result_135))*ratio_135


ratio_095 = ridge_integrate[0]/theor_integrate[0]
ratio_105 = ridge_integrate[1]/theor_integrate[1]
ratio_115 = ridge_integrate[2]/theor_integrate[2]
ratio_125 = ridge_integrate[3]/theor_integrate[3]
ratio_135 = ridge_integrate[4]/theor_integrate[4]

theory_norm_095 = (atlas_result_095-min(atlas_result_095))*ratio_095
theory_norm_105 = (atlas_result_105-min(atlas_result_105))*ratio_105
theory_norm_115 = (atlas_result_115-min(atlas_result_115))*ratio_115
theory_norm_125 = (atlas_result_125-min(atlas_result_125))*ratio_125
theory_norm_135 = (atlas_result_135-min(atlas_result_135))*ratio_135

print(ratio_095, ratio_105, ratio_115, ratio_125, ratio_135)

plt.plot(atlas_phi_095, theory_norm_095, color = color_ridge[0], linewidth=5, linestyle = '-',label=r'$90 \,\, \leq N_{ch}^{rec}<100$')
plt.plot(atlas_phi_105, theory_norm_105, color = color_ridge[1], linewidth=5, linestyle = '-',label=r'$100 \leq N_{ch}^{rec}<110$')
plt.plot(atlas_phi_115, theory_norm_115, color = color_ridge[2], linewidth=5, linestyle = '-',label=r'$110 \leq N_{ch}^{rec}<120$')
plt.plot(atlas_phi_125, theory_norm_125, color = color_ridge[3], linewidth=5, linestyle = '-',label=r'$120 \leq N_{ch}^{rec}<130$')
plt.plot(atlas_phi_135, theory_norm_135, color = color_ridge[4], linewidth=5, linestyle = '-',label=r'$130 \leq N_{ch}^{rec}$')

plt.ylabel(r'$Y^{ridge}-C_{ZYAM}$', size = 70)
plt.xlabel(r'$\Delta\phi$', size=70)

plt.minorticks_on()

# plt.ylim(-0.01,0.1)
# plt.xlim(-1.3,1.3)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.legend(fontsize=45,framealpha=False,bbox_to_anchor=(1.13,0.4), ncol=2)
plt.legend(fontsize=45,framealpha=False,loc='upper right')


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/only_ridges.png')

fig.clear()


phi_13TeV_90_up=np.loadtxt('./atlasgraphs/13TeV_90~.csv',delimiter=',',usecols=[0])
data_13TeV_90_up=np.loadtxt('./atlasgraphs/13TeV_90~.csv',delimiter=',',usecols=[1])
atlas_phi_rezero=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV_CMS_ALICE/phiCorrelation_pt0-5.csv',delimiter=',',usecols=[0],skiprows=1)
atlas_dat_rezero=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV_CMS_ALICE/phiCorrelation_pt0-5.csv',delimiter=',',usecols=[3],skiprows=1)

plt.scatter(phi_13TeV_90_up, data_13TeV_90_up-min(data_13TeV_90_up), color = 'black', s=1000, marker='o',label=r'$Y^{data}-C_{ZYAM}$')
plt.plot(atlas_phi_rezero, atlas_dat_rezero-min(atlas_dat_rezero), color = 'blue', linewidth=5, linestyle = '-',label=r'$Y^{theory}-C_{ZYAM}$')

plt.title(r'$90 \leq N_{ch}^{rec}$', size = 70)
plt.ylabel(r'$Y^{ridge}-C_{ZYAM}$', size = 70)
plt.xlabel(r'$\Delta\phi$', size=70)

plt.minorticks_on()

plt.ylim(-0.001,0.04)
plt.xlim(-1.,1.)

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.legend(fontsize=45,framealpha=False,bbox_to_anchor=(1.13,0.4), ncol=2)
# plt.legend(fontsize=45,framealpha=False,loc='upper right')


plt.grid(color='silver',linestyle=':',linewidth=3)

plt.tight_layout()

fig.savefig('./result/rezero_atlas.png')

fig.clear()

fig1, axes1 = plt.subplots(nrows=1, ncols=5,figsize=(100,20),sharey='row')

axes1[0].scatter(phi_13TeV_90_100, data_13TeV_90_100-min(data_13TeV_90_100), color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
axes1[1].scatter(phi_13TeV_100_110, data_13TeV_100_110-min(data_13TeV_100_110), color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
axes1[2].scatter(phi_13TeV_110_120, data_13TeV_110_120-min(data_13TeV_110_120), color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
axes1[3].scatter(phi_13TeV_120_130, data_13TeV_120_130-min(data_13TeV_120_130), color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
axes1[4].scatter(phi_13TeV_130_up, data_13TeV_130_up-min(data_13TeV_130_up), color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
axes1[0].plot(atlas_phi_095, atlas_result_095-min(atlas_result_095), color = color_ridge[0], linewidth=5, linestyle = '-',label=r'$90 \,\, \leq N_{ch}^{rec}<100$')
axes1[1].plot(atlas_phi_105, atlas_result_105-min(atlas_result_105), color = color_ridge[1], linewidth=5, linestyle = '-',label=r'$100 \leq N_{ch}^{rec}<110$')
axes1[2].plot(atlas_phi_115, atlas_result_115-min(atlas_result_115), color = color_ridge[2], linewidth=5, linestyle = '-',label=r'$110 \leq N_{ch}^{rec}<120$')
axes1[3].plot(atlas_phi_125, atlas_result_125-min(atlas_result_125), color = color_ridge[3], linewidth=5, linestyle = '-',label=r'$120 \leq N_{ch}^{rec}<130$')
axes1[4].plot(atlas_phi_135, atlas_result_135-min(atlas_result_135), color = color_ridge[4], linewidth=5, linestyle = '-',label=r'$130 \leq N_{ch}^{rec}$')

axes1[0].set_title(r'$90 \leq N^{rec}_{ch} < 100$', size = 70, pad=30)
axes1[1].set_title(r'$100 \leq N^{rec}_{ch} < 110$', size = 70, pad=30)
axes1[2].set_title(r'$110 \leq N^{rec}_{ch} < 120$', size = 70, pad=30)
axes1[3].set_title(r'$120 \leq N^{rec}_{ch} < 130$', size = 70, pad=30)
axes1[4].set_title(r'$130 \leq N^{rec}_{ch}$', size = 70, pad=30)

axes1[0].set_ylabel(r'$Y^{ridge}-C_{ZYAM}$', size = 70)
for i in range(5):
    axes1[i].set_xlabel(r'$\Delta\phi$', size=70)

    axes1[i].minorticks_on()

    axes1[i].set_ylim(-0.001,0.07)
    axes1[i].set_xlim(-1.,1.)

    axes1[i].tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
    axes1[i].tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
    axes1[i].grid(color='silver',linestyle=':',linewidth=3)

fig1.tight_layout(h_pad = -12)

fig1.savefig('./result/rezero_atlas_multi.png')
# , transparent = True

fig1.clear()

axes1[0].scatter(phi_13TeV_90_100, data_13TeV_90_100-min(data_13TeV_90_100), color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
axes1[1].scatter(phi_13TeV_100_110, data_13TeV_100_110-min(data_13TeV_100_110), color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
axes1[2].scatter(phi_13TeV_110_120, data_13TeV_110_120-min(data_13TeV_110_120), color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
axes1[3].scatter(phi_13TeV_120_130, data_13TeV_120_130-min(data_13TeV_120_130), color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
axes1[4].scatter(phi_13TeV_130_up, data_13TeV_130_up-min(data_13TeV_130_up), color = 'black', s=1000, marker='o',label=r'$Y^{data}$')
axes1[0].plot(atlas_phi_095, atlas_result_095-min(atlas_result_095), color = color_ridge[0], linewidth=5, linestyle = '-',label=r'$90 \,\, \leq N_{ch}^{rec}<100$')
axes1[1].plot(atlas_phi_105, atlas_result_105-min(atlas_result_105), color = color_ridge[1], linewidth=5, linestyle = '-',label=r'$100 \leq N_{ch}^{rec}<110$')
axes1[2].plot(atlas_phi_115, atlas_result_115-min(atlas_result_115), color = color_ridge[2], linewidth=5, linestyle = '-',label=r'$110 \leq N_{ch}^{rec}<120$')
axes1[3].plot(atlas_phi_125, atlas_result_125-min(atlas_result_125), color = color_ridge[3], linewidth=5, linestyle = '-',label=r'$120 \leq N_{ch}^{rec}<130$')
axes1[4].plot(atlas_phi_135, atlas_result_135-min(atlas_result_135), color = color_ridge[4], linewidth=5, linestyle = '-',label=r'$130 \leq N_{ch}^{rec}$')

axes1[0].set_title(r'$90 \leq N^{rec}_{ch} < 100$', size = 70, pad=30)
axes1[1].set_title(r'$100 \leq N^{rec}_{ch} < 110$', size = 70, pad=30)
axes1[2].set_title(r'$110 \leq N^{rec}_{ch} < 120$', size = 70, pad=30)
axes1[3].set_title(r'$120 \leq N^{rec}_{ch} < 130$', size = 70, pad=30)
axes1[4].set_title(r'$130 \leq N^{rec}_{ch}$', size = 70, pad=30)

axes1[0].set_ylabel(r'$Y^{ridge}-C_{ZYAM}$', size = 70)
for i in range(5):
    axes1[i].set_xlabel(r'$\Delta\phi$', size=70)

    axes1[i].minorticks_on()

    axes1[i].set_ylim(-0.001,0.06)
    axes1[i].set_xlim(-1.,1.)

    axes1[i].tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
    axes1[i].tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
    axes1[i].grid(color='silver',linestyle=':',linewidth=3)

fig1.tight_layout(h_pad = -12)

# fig1.savefig('./result/rezero_atlas_q.png')
# , transparent = True

fig1.clear()


