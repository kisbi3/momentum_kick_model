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

# alice 논문에서는 near-side 정의를 -1.28<phi<1.28으로 두었다.
# fig.set_size_inches(50, 30, forward=True)
ali_dat=[]
ali_phi=[]
ali_err=[]
ali_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=129, max_rows=13))
ali_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=167, max_rows=13))
ali_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=205, max_rows=13))
ali_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[3],skiprows=12,max_rows=5))
ali_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=129, max_rows=13))
ali_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=167, max_rows=13))
ali_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=205, max_rows=13))
ali_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[0],skiprows=12,max_rows=5))
err_sta1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=129, max_rows=13)
err_sta2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=129, max_rows=13)
err_sys1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=129, max_rows=13)
err_sys2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=129, max_rows=13)
ali_err.append((err_sta1**2+err_sys1**2)**0.5)
ali_err.append((err_sta2**2+err_sys2**2)**0.5)
err_sta1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=167, max_rows=13)
err_sta2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=167, max_rows=13)
err_sys1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=167, max_rows=13)
err_sys2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=167, max_rows=13)
ali_err.append((err_sta1**2+err_sys1**2)**0.5)
ali_err.append((err_sta2**2+err_sys2**2)**0.5)
err_sta1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=205, max_rows=13)
err_sta2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=205, max_rows=13)
err_sys1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=205, max_rows=13)
err_sys2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=205, max_rows=13)
ali_err.append((err_sta1**2+err_sys1**2)**0.5)
ali_err.append((err_sta2**2+err_sys2**2)**0.5)
err_sta1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[4],skiprows=12,max_rows=5)
err_sta2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[5],skiprows=12,max_rows=5)
err_sys1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[6],skiprows=12,max_rows=5)
err_sys2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[7],skiprows=12,max_rows=5)
ali_err.append((err_sta1**2+err_sys1**2)**0.5)
ali_err.append((err_sta2**2+err_sys2**2)**0.5)

#pt 1~4를 해야하나?
cms_dat=[]
cms_phi=[]
cms_err=[]
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13))
cms_dat.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[1],skiprows=14, max_rows= 4))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13))
cms_phi.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[0],skiprows=14, max_rows= 4))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13))
cms_err.append(np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13))
err_sta1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[2],skiprows=14,max_rows=4)
err_sta2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[3],skiprows=14,max_rows=4)
err_sys1=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[4],skiprows=14,max_rows=4)
err_sys2=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[5],skiprows=14,max_rows=4)
cms_err.append((err_sta1**2+err_sys1**2)**0.5)
cms_err.append((err_sta2**2+err_sys2**2)**0.5)

#ALICE CMS ALICE CMS ... 순서
mom_dat=[]
mom_phi=[]
mom_dat.append(np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[1],skiprows=1))
mom_dat.append(np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[1],skiprows=1))
mom_dat.append(np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[1],skiprows=1))
mom_dat.append(np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[2],skiprows=1))
mom_dat.append(np.loadtxt('pTdis.csv',delimiter=',',usecols=[1],skiprows=1))
mom_dat.append(np.loadtxt('pTdis.csv',delimiter=',',usecols=[2],skiprows=1))
mom_phi.append(np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[0],skiprows=1))
mom_phi.append(np.loadtxt('pTdis.csv',delimiter=',',usecols=[0],skiprows=1))
ali_integral = np.trapz(ali_dat[3], x = ali_phi[3])
cms_integral = np.trapz(cms_dat[3], x = cms_phi[3])
mom_ali_integral = np.trapz(mom_dat[6], x = mom_phi[3])
mom_cms_integral = np.trapz(mom_dat[7], x = mom_phi[3])

fig1, axes1 = plt.subplots(nrows=1, ncols=5,figsize=(100,20))


for i in range(3):
        # print(i)
        # print(len(mom_phi[2*i+1]))
        # print(len(mom_dat[2*i+1]))
        axes1[i].plot(mom_phi[i], mom_dat[2*i] - min(mom_dat[2*i]), color = "red", linewidth=7, linestyle='-')
        axes1[i].plot(mom_phi[i], mom_dat[2*i+1] - min(mom_dat[2*i+1]), color = "black", linewidth=7, linestyle='-')
        axes1[i].errorbar(ali_phi[i], ali_dat[i]-min(ali_dat[i]), yerr=(ali_err[2*i],abs(ali_err[2*i+1])), color="red", markersize=20, marker='o', linestyle=' ', fillstyle='none', linewidth=5, capsize=15)
        axes1[i].errorbar(cms_phi[i], cms_dat[i]-min(cms_dat[i]), yerr=(cms_err[2*i],abs(cms_err[2*i+1])), color="black", markersize=20, marker='o', linestyle=' ', fillstyle='none', linewidth=5, capsize=15)
        axes1[i].set_xlabel(r'$\Delta\phi$', size=70)

        axes1[i].minorticks_on()

        axes1[i].set_ylim(-0.001,0.02)
        axes1[i].set_xlim(-1.28,1.28)

        axes1[i].tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
        axes1[i].tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
        axes1[i].grid(color='silver',linestyle=':',linewidth=3)

        axes1[i].legend(framealpha=False, fontsize = 70)

fig1.savefig('./paper_graph/Paper.png')
# x = np.arange(0,4,0.01)

# f1 = 0.32*np.exp(0.93*x)        #7TeV frnk
# f2 = 0.66*np.exp(0.71*x)        #13TeV frnk
# f3 = np.zeros(len(x))                          #200GeV

# plt.plot(x,f1, color = 'red', linewidth=7, label=r'$pp, \, 7TeV$')
# plt.plot(x,f2, color = 'green', linewidth=7, label=r'$pp, \, 13TeV$')
# plt.plot(x,f3+4, color = 'blue', linewidth=7, label=r'$AuAu, \, 0.2TeV$')

# plt.xlabel(r'$p_T$',size=70)
# plt.ylabel(r'$f_{R} \langle N \rangle $',size=70)

# plt.xlim(0,4)

# ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

# plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top='true')
# plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top='true')

# plt.grid(color='silver',linestyle=':',linewidth=5)
# plt.legend(fontsize=45,framealpha=False, loc='upper left')

# plt.tight_layout()

# fig.savefig('./paper_graph/frnk.png')