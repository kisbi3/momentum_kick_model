import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit
import csv

mpl.rcParams["text.usetex"] = True

#atlas data들
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

#계산값들
q_080_phi=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q0.800000.csv',delimiter=',',usecols=[0], skiprows=1)
q_080_atl=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q0.800000.csv',delimiter=',',usecols=[3], skiprows=1)
q_085_phi=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q0.850000.csv',delimiter=',',usecols=[0], skiprows=1)
q_085_atl=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q0.850000.csv',delimiter=',',usecols=[3], skiprows=1)
q_090_phi=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q0.900000.csv',delimiter=',',usecols=[0], skiprows=1)
q_090_atl=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q0.900000.csv',delimiter=',',usecols=[3], skiprows=1)
q_095_phi=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q0.950000.csv',delimiter=',',usecols=[0], skiprows=1)
q_095_atl=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q0.950000.csv',delimiter=',',usecols=[3], skiprows=1)
q_100_phi=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.000000.csv',delimiter=',',usecols=[0], skiprows=1)
q_100_atl=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.000000.csv',delimiter=',',usecols=[3], skiprows=1)
q_105_phi=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.050000.csv',delimiter=',',usecols=[0], skiprows=1)
q_105_atl=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.050000.csv',delimiter=',',usecols=[3], skiprows=1)
q_110_phi=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.100000.csv',delimiter=',',usecols=[0], skiprows=1)
q_110_atl=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.100000.csv',delimiter=',',usecols=[3], skiprows=1)
q_120_phi=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.200000.csv',delimiter=',',usecols=[0], skiprows=1)
q_120_atl=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.200000.csv',delimiter=',',usecols=[3], skiprows=1)
q_130_phi=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.300000.csv',delimiter=',',usecols=[0], skiprows=1)
q_130_atl=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.300000.csv',delimiter=',',usecols=[3], skiprows=1)
q_140_phi=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.400000.csv',delimiter=',',usecols=[0], skiprows=1)
q_140_atl=np.loadtxt('./phiCorr_q/phiCorrelation_pt0-5_q1.400000.csv',delimiter=',',usecols=[3], skiprows=1)



fig1, axes1 = plt.subplots(nrows=1, ncols=5,figsize=(100,20),sharey='row')

axes1[0].scatter(phi_13TeV_90_100, data_13TeV_90_100-min(data_13TeV_90_100), color = 'black', s=1000, marker='o')
axes1[1].scatter(phi_13TeV_100_110, data_13TeV_100_110-min(data_13TeV_100_110), color = 'black', s=1000, marker='o')
axes1[2].scatter(phi_13TeV_110_120, data_13TeV_110_120-min(data_13TeV_110_120), color = 'black', s=1000, marker='o')
axes1[3].scatter(phi_13TeV_120_130, data_13TeV_120_130-min(data_13TeV_120_130), color = 'black', s=1000, marker='o')
axes1[4].scatter(phi_13TeV_130_up, data_13TeV_130_up-min(data_13TeV_130_up), color = 'black', s=1000, marker='o')

#양쪽 끝 색은 파란색, 파란색으로 그래프 사이 채우기, 그리고 중간의 그래프는 빨강!
color_q = ['blue', 'red']

#90<N<100
axes1[0].plot(q_080_phi, q_080_atl-min(q_080_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$low\,:\,q\,=\,0.8$')
axes1[0].plot(q_090_phi, q_090_atl-min(q_090_atl), color = color_q[1], linewidth=7, linestyle = '-',label=r'$mid\,:\,q\,=\,0.9$')
axes1[0].plot(q_100_phi, q_100_atl-min(q_100_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$high\,:\,q\,=\,1.0$')
axes1[0].fill_between(q_080_phi, q_080_atl-min(q_080_atl), q_100_atl-min(q_100_atl), color=color_q[0], alpha=0.4)
#100<N<110
axes1[1].plot(q_085_phi, q_085_atl-min(q_085_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$low\,:\,q\,=\,0.85$')
axes1[1].plot(q_095_phi, q_095_atl-min(q_095_atl), color = color_q[1], linewidth=7, linestyle = '-',label=r'$mid\,:\,q\,=\,0.95$')
axes1[1].plot(q_105_phi, q_105_atl-min(q_105_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$high\,:\,q\,=\,1.05$')
axes1[1].fill_between(q_085_phi, q_085_atl-min(q_085_atl), q_105_atl-min(q_105_atl), color=color_q[0], alpha=0.4)
#110<N<120
axes1[2].plot(q_100_phi, q_100_atl-min(q_100_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$low\,:\,q\,=\,1.0$')
axes1[2].plot(q_110_phi, q_110_atl-min(q_110_atl), color = color_q[1], linewidth=7, linestyle = '-',label=r'$mid\,:\,q\,=\,1.1$')
axes1[2].plot(q_120_phi, q_120_atl-min(q_120_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$high\,:\,q\,=\,1.2$')
axes1[2].fill_between(q_100_phi, q_100_atl-min(q_100_atl), q_120_atl-min(q_120_atl), color=color_q[0], alpha=0.4)
#120<N<130
axes1[3].plot(q_110_phi, q_110_atl-min(q_110_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$low\,:\,q\,=\,1.1$')
axes1[3].plot(q_120_phi, q_120_atl-min(q_120_atl), color = color_q[1], linewidth=7, linestyle = '-',label=r'$mid\,:\,q\,=\,1.2$')
axes1[3].plot(q_130_phi, q_130_atl-min(q_130_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$high\,:\,q\,=\,1.3$')
axes1[3].fill_between(q_110_phi, q_110_atl-min(q_110_atl), q_130_atl-min(q_130_atl), color=color_q[0], alpha=0.4)
#130<N
axes1[4].plot(q_120_phi, q_120_atl-min(q_120_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$low\,:\,q\,=\,1.2$')
axes1[4].plot(q_130_phi, q_130_atl-min(q_130_atl), color = color_q[1], linewidth=7, linestyle = '-',label=r'$mid\,:\,q\,=\,1.3$')
axes1[4].plot(q_140_phi, q_140_atl-min(q_140_atl), color = color_q[0], linewidth=5, linestyle = '-',label=r'$high\,:\,q\,=\,1.4$')
axes1[4].fill_between(q_120_phi, q_120_atl-min(q_120_atl), q_140_atl-min(q_140_atl), color=color_q[0], alpha=0.4)

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

    axes1[i].legend(framealpha=False, fontsize = 70)

fig1.tight_layout(h_pad = -12)

fig1.savefig('./result/rezero_atlas_q.png')
# , transparent = True

fig1.clear()


