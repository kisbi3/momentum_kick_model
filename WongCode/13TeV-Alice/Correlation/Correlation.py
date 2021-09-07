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

# alice 논문에서는 near-side 정의를 -1.28<phi<1.28으로 두었다. 따라서 데이터는 약 -1.18<phi<1.18를 이용할 것이다.

fig.set_size_inches(35, 16.534, forward=True)

deltaphitable27=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=129, max_rows=13)
dNdphitable27=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=129, max_rows=13)
datamintable27=min(dNdphitable27)
table27_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=129, max_rows=13)
table27_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=129, max_rows=13)

table27_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=129, max_rows=13)
table27_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=129, max_rows=13)

table27_error1 = (table27_error1_stat**2+table27_error1_sys**2)**0.5
table27_error2 = (table27_error2_stat**2+table27_error2_sys**2)**0.5

resultphi=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[1],skiprows=1)
dNdphimin=min(resultdNdphi)

print(dNdphimin)


plt.errorbar(deltaphitable27,dNdphitable27-datamintable27, yerr=(table27_error1,abs(table27_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,13TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$result,\quad C_{ZYAM}=0.00501599$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

plt.xlabel(r'$\Delta\phi$',size=50)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
plt.title(r'$1.0<p_T<2.0$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_Correlation_czyam_pt1-2.png')

fig.clear()

# fig = plt.figure()
# ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

deltaphitable29=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=167, max_rows=13)
dNdphitable29=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=167, max_rows=13)
datamintable29=min(dNdphitable29)
table29_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=167, max_rows=13)
table29_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=167, max_rows=13)

table29_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=167, max_rows=13)
table29_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=167, max_rows=13)

table29_error1 = (table29_error1_stat**2+table29_error1_sys**2)**0.5
table29_error2 = (table29_error2_stat**2+table29_error2_sys**2)**0.5

resultphi=np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[1],skiprows=1)
dNdphimin=min(resultdNdphi)

print(dNdphimin)


plt.errorbar(deltaphitable29,dNdphitable29-datamintable29, yerr=(table29_error1,abs(table29_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,13TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$result,\quad C_{ZYAM}=0.00501599$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

plt.xlabel(r'$\Delta\phi$',size=50)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
plt.title(r'$2.0<p_T<3.0$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_Correlation_czyam_pt2-3.png')

fig.clear()

# fig = plt.figure()
# ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

deltaphitable31=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=205, max_rows=13)
dNdphitable31=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=205, max_rows=13)
datamintable31=min(dNdphitable31)
table31_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=205, max_rows=13)
table31_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=205, max_rows=13)

table31_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=205, max_rows=13)
table31_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=205, max_rows=13)

table31_error1 = (table31_error1_stat**2+table31_error1_sys**2)**0.5
table31_error2 = (table31_error2_stat**2+table31_error2_sys**2)**0.5

resultphi=np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[1],skiprows=1)
dNdphimin=min(resultdNdphi)

print(dNdphimin)


plt.errorbar(deltaphitable31,dNdphitable31-datamintable31, yerr=(table31_error1,abs(table31_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,13TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$result,\quad C_{ZYAM}=0.00501599$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

plt.xlabel(r'$\Delta\phi$',size=50)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
plt.title(r'$3.0<p_T<4.0$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_Correlation_czyam_pt3-4.png')

fig.clear()




# resultphi=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[0],skiprows=1)
# resultdNdphi=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[1],skiprows=1)
# dNdphimin=min(resultdNdphi)


# deltaphi_pt14=deltaphitable27
# dNdphi_pt14=(dNdphitable27+1.27)+(dNdphitable29+0.28)+(dNdphitable31+0.09)
#         # C_zyam 처리해야함 errorbar는 어떻게 할까
# mindNdphi_pt14=min(dNdphi_pt14)

# plt.errorbar(deltaphi_pt14,dNdphi_pt14-mindNdphi_pt14, yerr=((table27_error1**2+table29_error1**2+table31_error1**2)**0.5,(table27_error2**2+table29_error2**2+table31_error2**2)**0.5), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,13TeV \, CMS$',capsize=10)
# # ,fillstyle='none'
# # plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# # plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$result,\quad C_{ZYAM}=0.00501599$')
# # plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

# plt.xlabel(r'$\Delta\phi$',size=50)
# plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# # plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
# plt.title(r'$1.0<p_T<4.0$', fontsize = 60)
# plt.minorticks_on()
# # plt.yscale('log')
# # ax.axis([-0.1,3.1,0,0.1])

# ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

# plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
# plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


# plt.grid(color='silver',linestyle=':',linewidth=3)
# # plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
# # plt.legend(fontsize=45)
# plt.tight_layout()

# fig.savefig('1D_Correlation_czyam_pt1-4.png')

# fig.clear()


