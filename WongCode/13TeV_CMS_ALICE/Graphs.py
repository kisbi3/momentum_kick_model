import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import csv
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

mpl.rcParams["text.usetex"] = True

errors = 0

f = open('deviation.csv','w',newline='')
wr=csv.writer(f)

fig = plt.figure()
ax = plt.axes()

# alice 논문에서는 near-side 정의를 -1.28<phi<1.28으로 두었다.

fig.set_size_inches(35, 16.534, forward=True)

alice_deltaphitable27=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=129, max_rows=13)
alice_dNdphitable27=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=129, max_rows=13)
alice_datamintable27=min(alice_dNdphitable27)
alice_table27_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=129, max_rows=13)
alice_table27_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=129, max_rows=13)

alice_table27_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=129, max_rows=13)
alice_table27_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=129, max_rows=13)

alice_table27_error1 = (alice_table27_error1_stat**2+alice_table27_error1_sys**2)**0.5
alice_table27_error2 = (alice_table27_error2_stat**2+alice_table27_error2_sys**2)**0.5

cms_deltaphitable27=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
cms_dNdphitable27=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
cms_datamintable27=min(cms_dNdphitable27)
cms_table27_error1=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
cms_table27_error2=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table27.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

resultphi=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[0],skiprows=1)
alice_resultdNdphi=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[1],skiprows=1)
cms_resultdNdphi=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[2],skiprows=1)

alice_dNdphimin=min(alice_resultdNdphi)
cms_dNdphimin=min(cms_resultdNdphi)

# print(dNdphimin)
alice_czyam = alice_dNdphitable27-alice_datamintable27
alice_result_czyam = alice_resultdNdphi-alice_dNdphimin
cms_czyam = cms_dNdphitable27-cms_datamintable27
cms_result_czyam = cms_resultdNdphi-cms_dNdphimin
j = 4
alicemse = 0
cmsmse = 0
for i in range(len(alice_dNdphitable27)):
        # print(cms_deltaphitable27[i],alice_deltaphitable27[i],resultphi[round(j)])
        alicemse += abs(alice_czyam[i]-alice_result_czyam[round(j)])**2
        cmsmse += (cms_czyam[i]-cms_result_czyam[round(j)])**2
        j = 4+(i+1)*7.6699
alicemse = (alicemse/len(alice_dNdphitable27))**0.5
cmsmse = (cmsmse/len(cms_dNdphitable27))**0.5
print('1<pT<2')
print('ALICE error = ', alicemse)
print('CMS error = ', cmsmse)
wr.writerow(['1<pT<2','ALICE : ',alicemse])
wr.writerow(['1<pT<2','CMS : ',cmsmse])

errors += alicemse + cmsmse

# print(alicemse/len(alice_dNdphitable27))
# print(cmsmse/len(cms_dNdphitable27))


plt.errorbar(alice_deltaphitable27,alice_dNdphitable27-alice_datamintable27, yerr=(alice_table27_error1,abs(alice_table27_error2)), color="blue",markersize=20,marker='o',linestyle=' ',fillstyle='none',linewidth=5,label=r'$pp,13TeV \, ALICE$',capsize=10)
plt.errorbar(cms_deltaphitable27,cms_dNdphitable27-cms_datamintable27, yerr=(cms_table27_error1,abs(cms_table27_error2)), color="magenta",markersize=25,marker='o',linestyle=' ',fillstyle='none',linewidth=5,label=r'$pp,13TeV \, CMS$',capsize=15)

# msg1 = 'alice error = ', alicemse}
# print(msg1)


# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
plt.plot(resultphi, alice_resultdNdphi-alice_dNdphimin, color = 'blue', linewidth=7, linestyle = '-',label=fr'$result,\, ALICE,\, Deviation : {alicemse:.3E}$')
# plt.text(-0.5, 0.0025,fr"$ALICE\quad Error : {alicemse:.3E}$Deviationze=50)
# plt.text(-0.5, 0.,fr"$CMS\quad Error : {cmsmse:.3E}$", size=50)
plt.plot(resultphi, cms_resultdNdphi-cms_dNdphimin, color = 'magenta', linewidth=7, linestyle = '--',label=fr'$result,\, CMS,\, Deviation : {cmsmse:.3E}$')
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
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_PhiCorrelation_pt1-2.png')

fig.clear()

# fig = plt.figure()
# ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

alice_deltaphitable29=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=167, max_rows=13)
alice_dNdphitable29=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=167, max_rows=13)
alice_datamintable29=min(alice_dNdphitable29)
alice_table29_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=167, max_rows=13)
alice_table29_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=167, max_rows=13)

alice_table29_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=167, max_rows=13)
alice_table29_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=167, max_rows=13)

alice_table29_error1 = (alice_table29_error1_stat**2+alice_table29_error1_sys**2)**0.5
alice_table29_error2 = (alice_table29_error2_stat**2+alice_table29_error2_sys**2)**0.5

cms_deltaphitable29=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
cms_dNdphitable29=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
cms_datamintable29=min(cms_dNdphitable29)
cms_table29_error1=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
cms_table29_error2=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table29.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

resultphi=np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[0],skiprows=1)
alice_resultdNdphi=np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[1],skiprows=1)
cms_resultdNdphi=np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[2],skiprows=1)

alice_dNdphimin=min(alice_resultdNdphi)
cms_dNdphimin=min(cms_resultdNdphi)


# print(dNdphimin)
alice_czyam = alice_dNdphitable29-alice_datamintable29
alice_result_czyam = alice_resultdNdphi-alice_dNdphimin
cms_czyam = cms_dNdphitable29-cms_datamintable29
cms_result_czyam = cms_resultdNdphi-cms_dNdphimin
j = 4
alicemse = 0
cmsmse = 0
for i in range(len(alice_dNdphitable29)):
        # print(cms_deltaphitable27[i],alice_deltaphitable27[i],resultphi[round(j)])
        alicemse += abs(alice_czyam[i]-alice_result_czyam[round(j)])**2
        cmsmse += (cms_czyam[i]-cms_result_czyam[round(j)])**2
        j = 4+(i+1)*7.6699
alicemse = (alicemse/len(alice_dNdphitable29))**0.5
cmsmse = (cmsmse/len(cms_dNdphitable29))**0.5
print('2<pT<3')
print('ALICE error = ', alicemse)
print('CMS error = ', cmsmse)
wr.writerow(['2<pT<3','ALICE : ',alicemse])
wr.writerow(['2<pT<3','CMS : ',cmsmse])

errors += alicemse + cmsmse

plt.errorbar(alice_deltaphitable29,alice_dNdphitable29-alice_datamintable29, yerr=(alice_table29_error1,abs(alice_table29_error2)), color="blue",markersize=20,marker='o',fillstyle='none',linestyle=' ',linewidth=5,label=r'$pp,13TeV \, ALICE$',capsize=15)
plt.errorbar(cms_deltaphitable29,cms_dNdphitable29-cms_datamintable29, yerr=(cms_table29_error1,abs(cms_table29_error2)), color="magenta",markersize=25,marker='o',linestyle=' ',linewidth=5,fillstyle='none',label=r'$pp,13TeV \, CMS$',capsize=15)

# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
plt.plot(resultphi, alice_resultdNdphi-alice_dNdphimin, color = 'blue', linewidth=7, linestyle = '-',label=fr'$result,\, ALICE,\, Deviation : {alicemse:.3E}$')
plt.plot(resultphi, cms_resultdNdphi-cms_dNdphimin, color = 'magenta', linewidth=7, linestyle = '--',label=fr'$result,\, CMS,\, Deviation : {cmsmse:.3E}$')
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
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_PhiCorrelation_pt2-3.png')

fig.clear()

# fig = plt.figure()
# ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

alice_deltaphitable31=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=205, max_rows=13)
alice_dNdphitable31=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=205, max_rows=13)
alice_datamintable31=min(alice_dNdphitable31)
alice_table31_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=205, max_rows=13)
alice_table31_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=205, max_rows=13)

alice_table31_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=205, max_rows=13)
alice_table31_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=205, max_rows=13)

alice_table31_error1 = (alice_table31_error1_stat**2+alice_table31_error1_sys**2)**0.5
alice_table31_error2 = (alice_table31_error2_stat**2+alice_table31_error2_sys**2)**0.5

cms_deltaphitable31=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
cms_dNdphitable31=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
cms_datamintable31=min(cms_dNdphitable31)
cms_table31_error1=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
cms_table31_error2=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table31.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

resultphi=np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[0],skiprows=1)
alice_resultdNdphi=np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[1],skiprows=1)
cms_resultdNdphi=np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[2],skiprows=1)

alice_dNdphimin=min(alice_resultdNdphi)
cms_dNdphimin=min(cms_resultdNdphi)

alice_czyam = alice_dNdphitable31-alice_datamintable31
alice_result_czyam = alice_resultdNdphi-alice_dNdphimin
cms_czyam = cms_dNdphitable31-cms_datamintable31
cms_result_czyam = cms_resultdNdphi-cms_dNdphimin
j = 4
alicemse = 0
cmsmse = 0
for i in range(len(alice_dNdphitable31)):
        # print(cms_deltaphitable27[i],alice_deltaphitable27[i],resultphi[round(j)])
        alicemse += abs(alice_czyam[i]-alice_result_czyam[round(j)])**2
        cmsmse += (cms_czyam[i]-cms_result_czyam[round(j)])**2
        j = 4+(i+1)*7.6699
alicemse = (alicemse/len(alice_dNdphitable31))**0.5
cmsmse = (cmsmse/len(cms_dNdphitable31))**0.5
print('3<pT<4')
print('ALICE error = ', alicemse)
print('CMS error = ', cmsmse)

errors += alicemse + cmsmse

wr.writerow(['3<pT<4','ALICE : ',alicemse])
wr.writerow(['3<pT<4','CMS : ',cmsmse])

# print(dNdphimin)


plt.errorbar(alice_deltaphitable31,alice_dNdphitable31-alice_datamintable31, yerr=(alice_table31_error1,abs(alice_table31_error2)), color="blue",markersize=20,marker='o',linestyle=' ',fillstyle='none',linewidth=5,label=r'$pp,13TeV \, ALICE$',capsize=15)
plt.errorbar(cms_deltaphitable31,cms_dNdphitable31-cms_datamintable31, yerr=(cms_table31_error1,abs(cms_table31_error2)), color="magenta",markersize=25,marker='o',linestyle=' ', fillstyle='none',linewidth=5,label=r'$pp,13TeV \, CMS$',capsize=15)

# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
plt.plot(resultphi, alice_resultdNdphi-alice_dNdphimin, color = 'blue', linewidth=7, linestyle = '-',label=fr'$result,\,ALICE,\, Deviation : {alicemse:.3E}$')
plt.plot(resultphi, cms_resultdNdphi-cms_dNdphimin, color = 'magenta', linewidth=7, linestyle = '--',label=fr'$result,\,CMS,\, Deviation : {cmsmse:.3E}$')
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
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_PhiCorrelation_pt3-4.png')

fig.clear()

resultphi=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[0],skiprows=1)
alice_resultdNdphi=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[1],skiprows=1)
cms_resultdNdphi=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[2],skiprows=1)
cms_dNdphimin=min(cms_resultdNdphi)
alice_dNdphimin=min(alice_resultdNdphi)


deltaphi_pt14=cms_deltaphitable27
dNdphi_pt14=(cms_dNdphitable27+1.27)+(cms_dNdphitable29+0.28)+(cms_dNdphitable31+0.09)
        # C_zyam 처리해야함 errorbar는 어떻게 할까
mindNdphi_pt14=min(dNdphi_pt14)

cms_czyam = dNdphi_pt14-mindNdphi_pt14
cms_result_czyam = cms_resultdNdphi-cms_dNdphimin
j = 4
cmsmse = 0
for i in range(len(alice_dNdphitable31)):
        # print(cms_deltaphitable27[i],alice_deltaphitable27[i],resultphi[round(j)])
        cmsmse += (cms_czyam[i]-cms_result_czyam[round(j)])**2
        j = 4+(i+1)*7.6699
cmsmse = (cmsmse/len(cms_dNdphitable31))**0.5
print('1<pT<4')
print('CMS error = ', cmsmse)
wr.writerow(['1<pT<4','CMS : ',cmsmse])

errors += cmsmse

plt.errorbar(deltaphi_pt14,dNdphi_pt14-mindNdphi_pt14, yerr=((cms_table27_error1**2+cms_table29_error1**2+cms_table31_error1**2)**0.5,(cms_table27_error2**2+cms_table29_error2**2+cms_table31_error2**2)**0.5), color="magenta",markersize=25,fillstyle='none',marker='o',linestyle=' ',linewidth=5,label=r'$pp,13TeV \, CMS$',capsize=15)

# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
plt.plot(resultphi, alice_resultdNdphi-alice_dNdphimin, color = 'blue', linewidth=7, linestyle = '-',label=fr'$result,\,ALICE$')
plt.plot(resultphi, cms_resultdNdphi-cms_dNdphimin, color = 'magenta', linewidth=7, linestyle = '--',label=fr'$result,\,CMS,\, Deviation : {cmsmse:.3E}$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

plt.xlabel(r'$\Delta\phi$',size=50)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
plt.title(r'$1.0<p_T<4.0$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False,loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_PhiCorrelation_pt1-4.png')

fig.clear()


pt=np.loadtxt('pTdis.csv',delimiter=',',usecols=[0],skiprows=1)
aliceRidgedis=np.loadtxt('pTdis.csv',delimiter=',',usecols=[1],skiprows=1)
cmsRidgedis=np.loadtxt('pTdis.csv',delimiter=',',usecols=[2],skiprows=1)
alicejetdis=np.loadtxt('pTdis.csv',delimiter=',',usecols=[3],skiprows=1)
# cmsjetdis=np.loadtxt('pTdis.csv',delimiter=',',usecols=[4],skiprows=1)


Alice13TeVpt=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[0],skiprows=12,max_rows=5)
Alice13TeVptdis=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[3],skiprows=12,max_rows=5)
Alice13TeVptdis_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[4],skiprows=12,max_rows=5)
Alice13TeVptdis_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[5],skiprows=12,max_rows=5)
Alice13TeVptdis_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[6],skiprows=12,max_rows=5)
Alice13TeVptdis_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[7],skiprows=12,max_rows=5)

Alice13TeVptdis_error1 = pow(Alice13TeVptdis_error1_stat*Alice13TeVptdis_error1_stat+Alice13TeVptdis_error1_sys*Alice13TeVptdis_error1_sys,0.5)
Alice13TeVptdis_error2 = pow(Alice13TeVptdis_error2_stat*Alice13TeVptdis_error2_stat+Alice13TeVptdis_error2_sys*Alice13TeVptdis_error2_sys,0.5)

# print(Alice13TeVptdis_error1)
# print(Alice13TeVptdis_error2)

CMS13TeVpt=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[0],skiprows=14,max_rows=9)
CMS13TeVptdis=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[1],skiprows=14,max_rows=9)
CMS13TeVptdis_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[2],skiprows=14,max_rows=9)
CMS13TeVptdis_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[3],skiprows=14,max_rows=9)
CMS13TeVptdis_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[4],skiprows=14,max_rows=9)
CMS13TeVptdis_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[5],skiprows=14,max_rows=9)

CMS13TeVptdis_error1 = pow(CMS13TeVptdis_error1_stat*CMS13TeVptdis_error1_stat+CMS13TeVptdis_error1_sys*CMS13TeVptdis_error1_sys,0.5)
CMS13TeVptdis_error2 = pow(CMS13TeVptdis_error2_stat*CMS13TeVptdis_error2_stat+CMS13TeVptdis_error2_sys*CMS13TeVptdis_error2_sys,0.5)

alicemse=0
cmsmse = 0
alicemse = (Alice13TeVptdis[0]-aliceRidgedis[2])**2+(Alice13TeVptdis[1]-aliceRidgedis[3])**2+(Alice13TeVptdis[2]-aliceRidgedis[4])**2+(Alice13TeVptdis[3]-aliceRidgedis[5])**2
alicemse = (alicemse/4)**0.5

print('Y^Ridge')
for i in range(len(CMS13TeVpt)-2):
        # print(CMS13TeVpt[i],pt[i])
        cmsmse += (CMS13TeVptdis[i]-cmsRidgedis[i])**2
cmsmse += (CMS13TeVptdis[7]-cmsRidgedis[9])**2+(CMS13TeVptdis[8]-cmsRidgedis[20])**2
cmsmse = (cmsmse/9)**0.5

errors += alicemse + cmsmse

print('ALICE error = ',alicemse)
print('CMS error = ',cmsmse)

wr.writerow(['Y^Ridge','ALICE : ',alicemse])
wr.writerow(['Y^Ridge','CMS : ',cmsmse])

# yerr=(Alice13TeVptdis_error1,abs(Alice13TeVptdis_error2))
plt.errorbar(Alice13TeVpt, Alice13TeVptdis, yerr=(abs(Alice13TeVptdis_error2),Alice13TeVptdis_error1), color="blue",zorder=1,markersize=35,marker='v',linestyle=' ',linewidth=3,label=r'$result,\, ALICE, \,1.6<\vert\Delta\eta\vert<1.8$',capsize=10)
plt.errorbar(CMS13TeVpt, CMS13TeVptdis, yerr=(abs(CMS13TeVptdis_error2),CMS13TeVptdis_error1), color="black",zorder=1,markersize=35,marker='^',linestyle=' ',linewidth=3,label=r'$result,\, CMS, \,2.0<\vert\Delta\eta\vert<4.0$',capsize=10)

# plt.plot(pt,aliceRidgedis, color="blue",linewidth=5,linestyle = '-',label=r'$pp,Ridge,1.6<\vert\Delta\eta\vert<1.8$')
# plt.scatter(pt,aliceRidgedis, color="skyblue",s=1000,linewidths=6,facecolors='none',marker='D',label=fr'$pp,Ridge,1.6<\vert\Delta\eta\vert<1.8,\, Deviation : {alicemse:.3E}$',zorder=2)
# plt.scatter(pt,cmsRidgedis, color="crimson",s=1000,linewidths=6,facecolors='none',marker='s',label=fr'$pp,Ridge,2.0<\vert\Delta\eta\vert<4.0,\, Deviation : {cmsmse:.3E}$',zorder=2)


x=np.arange(0.45,10,0.01)
y1=0.06782985130264102*np.exp(-(0.9341561431595556/(x**2))-0.984985854370021*x)
y2=0.07314538280113343*np.exp((-0.8078484337720946/(x**2))-0.8777989424828683*x)

plt.plot(x,y1,color = 'skyblue', linewidth=7, linestyle = '-')
plt.plot(x,y2,color = 'red', linewidth=7, linestyle = '-')
# plt.plot(alicept,alicejetdis, color="blue",linewidth=5,linestyle = '--',label=r'$pp,Jet,1.6<\vert\Delta\eta\vert<1.8$')
# plt.plot(cmspt,cmsRidgedis, color="black",linewidth=5,linestyle = '-',label=r'$pp,Ridge,2.0<\vert\Delta\eta\vert<4.0$')


plt.xlabel(r'$p_{t}\quad(GeV)$',size=50)
# plt.ylabel(r'$(1/N_{trig})dN_{ch}/p_{t}dp_{t}{}\quad(GeV^{-2})$',size=50)
plt.ylabel(r'$Y^{Ridge}$',size=50)
# plt.plot(pt,jetdis,color="blue", linestyle = '-.', linewidth=7,label=r'$pp,Jet$')
# plt.plot(pt,Ridgedis+0.632*jetdis, color="black",linewidth=7,label=r'$p+p,Jet+Ridge$')

plt.minorticks_on()
# plt.yscale('log')
# ax.axis([0,4,0.0004,20])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
plt.tight_layout()

fig.savefig('pTdis.png')

print(f'average : {errors/9}')


fig.clear()

Alice13TeV_Y_near_jet_pt=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{near,JET}.csv',delimiter=',',usecols=[0],skiprows=12,max_rows=7)
Alice13TeV_Y_near_jet=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{near,JET}.csv',delimiter=',',usecols=[1],skiprows=12,max_rows=7)

Alice13TeV_Y_near_jet_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{near,JET}.csv',delimiter=',',usecols=[2],skiprows=12,max_rows=7)
Alice13TeV_Y_near_jet_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{near,JET}.csv',delimiter=',',usecols=[3],skiprows=12,max_rows=7)
Alice13TeV_Y_near_jet_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{near,JET}.csv',delimiter=',',usecols=[4],skiprows=12,max_rows=7)
Alice13TeV_Y_near_jet_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{near,JET}.csv',delimiter=',',usecols=[5],skiprows=12,max_rows=7)
Alice13TeV_Y_near_jet_min=min(Alice13TeV_Y_near_jet)

Alice13TeV_Y_near_jet_error1 = (Alice13TeV_Y_near_jet_error1_stat**2+Alice13TeV_Y_near_jet_error1_sys**2)**0.5
Alice13TeV_Y_near_jet_error2 = (Alice13TeV_Y_near_jet_error2_stat**2+Alice13TeV_Y_near_jet_error2_sys**2)**0.5

resultjet_phi=np.loadtxt('pTdis.csv',delimiter=',',usecols=[0],skiprows=1)
resultjet_dNdphi=np.loadtxt('pTdis.csv',delimiter=',',usecols=[3],skiprows=1)
resultjet_dNdphi_min=min(resultjet_dNdphi)


plt.errorbar(Alice13TeV_Y_near_jet_pt,Alice13TeV_Y_near_jet-Alice13TeV_Y_near_jet_min, yerr=(Alice13TeV_Y_near_jet_error1,Alice13TeV_Y_near_jet_error2), color="blue",markersize=20,marker='o',linestyle=' ',linewidth=5,label=r'$pp,13TeV \, ALICE$',capsize=15)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
plt.plot(resultjet_phi, resultjet_dNdphi-resultjet_dNdphi_min, color = 'blue', linewidth=7, linestyle = '-',label=r'$result$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

plt.xlabel(r'$p_{T,min}^{Jet}$',size=50)
plt.ylabel(r'$Y^{near}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
plt.title(r'$Y^{near}$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('Y^near.png')

fig.clear()

Alice13TeV_jetphi_over10=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=53,max_rows=13)
Alice13TeV_jetdNdphi_over10=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=53,max_rows=13)

Alice13TeV_jetdNdphi_over10_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=53,max_rows=13)
Alice13TeV_jetdNdphi_over10_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=53,max_rows=13)
Alice13TeV_jetdNdphi_over10_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=53,max_rows=13)
Alice13TeV_jetdNdphi_over10_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=53,max_rows=13)
Alice13TeV_jetdNdphi_over10_min=min(Alice13TeV_jetdNdphi_over10)

Alice13TeV_jetdNdphi_over10_error1 = (Alice13TeV_jetdNdphi_over10_error1_stat**2+Alice13TeV_jetdNdphi_over10_error1_sys**2)**0.5
Alice13TeV_jetdNdphi_over10_error2 = (Alice13TeV_jetdNdphi_over10_error2_stat**2+Alice13TeV_jetdNdphi_over10_error2_sys**2)**0.5

resultjet_phi=np.loadtxt('Jet_PhiCorrelation_pT_over10.csv',delimiter=',',usecols=[0],skiprows=1)
resultjet_dNdphi=np.loadtxt('Jet_PhiCorrelation_pT_over10.csv',delimiter=',',usecols=[1],skiprows=1)
resultjet_dNdphi_min=min(resultjet_dNdphi)


plt.errorbar(Alice13TeV_jetphi_over10,Alice13TeV_jetdNdphi_over10-Alice13TeV_jetdNdphi_over10_min, yerr=(Alice13TeV_jetdNdphi_over10_error1,Alice13TeV_jetdNdphi_over10_error2), color="blue",markersize=20,marker='o',linestyle=' ',linewidth=5,label=r'$pp,13TeV \, ALICE$',capsize=15)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
plt.plot(resultjet_phi, resultjet_dNdphi-resultjet_dNdphi_min, color = 'blue', linewidth=7, linestyle = '-',label=r'$result$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

plt.xlabel(r'$\Delta\phi$',size=50)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
plt.title(r'$1.0<p_T<4.0$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true', right = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true', right = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_Jet_PhiCorrelation_pt_over10.png')