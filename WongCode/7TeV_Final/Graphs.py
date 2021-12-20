import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import csv
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

mpl.rcParams["text.usetex"] = True

f = open('deviation.csv','w',newline='')
wr=csv.writer(f)

fig = plt.figure()
ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

# # fRNk=np.loadtxt('FittingParameters.csv',delimiter=',',usecols=[1],skiprows=8)

# deltaphitable26=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table26.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
# dNdphitable26=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table26.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
# datamintable26=min(dNdphitable26)
# table26_error1=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table26.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
# table26_error2=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table26.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

# # resultphi=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[0],skiprows=119, max_rows=237)
# # resultdNdphi=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[1],skiprows=119, max_rows=237)
# resultphi=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[0],skiprows=1)
# resultdNdphi=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[1],skiprows=1)
# resultJet=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[2],skiprows=1)
# dNdphimin=min(resultdNdphi)
# Ridge_Jetmin = min(resultdNdphi+resultJet)

# # plt.errorbar(deltaphitable26,dNdphitable26-datamintable26, yerr=(table26_error1,abs(table26_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# # plt.errorbar(deltaphitable26,dNdphitable26+4.77, yerr=(table26_error1,abs(table26_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# plt.errorbar(deltaphitable26,dNdphitable26, yerr=(table26_error1,abs(table26_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# # ,fillstyle='none'
# # plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# # plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# # frnk=fRNk[0]
# # plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
# # plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# # plt.plot(resultphi, resultdNdphi+resultJet-Ridge_Jetmin, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
# # plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

# plt.xlabel(r'$\Delta\phi$',size=50)
# plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# # plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
# plt.title(r'$0.1<p_T<1.0$', fontsize = 60)
# plt.minorticks_on()
# # plt.yscale('log')
# # ax.axis([-0.1,3.1,0,0.1])

# ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

# plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
# plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


# plt.grid(color='silver',linestyle=':',linewidth=3)
# # plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
# # plt.legend(fontsize=45)
# plt.tight_layout()

# fig.savefig('1D_Correlation_czyam_pt01-1.png')

# fig.clear()

errors = 0

deltaphitable28=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table28.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
dNdphitable28=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table28.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
datamintable28=min(dNdphitable28)
table28_error1=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table28.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
table28_error2=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table28.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

resultphi=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[1],skiprows=1)
resultJet=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[2],skiprows=1)
dNdphimin=min(resultdNdphi)
Ridge_Jetmin = min(resultdNdphi+resultJet)

# print(dNdphimin)

cms_czyam = dNdphitable28-datamintable28
cms_result_czyam = resultdNdphi-dNdphimin
j = 4
cmsmse = 0
for i in range(len(dNdphitable28)):
        # print(deltaphitable28[i], resultphi[round(j)])
        cmsmse += (cms_czyam[i]-cms_result_czyam[round(j)])**2
        j = 4+(i+1)*7.6699
cmsmse = (cmsmse/len(dNdphitable28))**0.5
print('1<pT<2')
print('CMS error = ', cmsmse)
wr.writerow(['1<pT<2',cmsmse])

errors +=  cmsmse

# plt.errorbar(deltaphitable28,dNdphitable28-datamintable28, yerr=(table28_error1,abs(table28_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# plt.errorbar(deltaphitable28,dNdphitable28+1.27, yerr=(table28_error1,abs(table28_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
plt.errorbar(deltaphitable28,dNdphitable28-datamintable28, yerr=(table28_error1,abs(table28_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet-Ridge_Jetmin, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=fr'$result, \, CMS$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

plt.xlabel(r'$\Delta\phi$',size=50)
# plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
# plt.title(r'$1.0<p_T<2.0$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False, loc='upper left')
plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1.5, 1))
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_Correlation_czyam_pt1-2.png')

fig.clear()

# fig = plt.figure()
# ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

deltaphitable30=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table30.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
dNdphitable30=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table30.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
datamintable30=min(dNdphitable30)
table30_error1=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table30.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
table30_error2=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table30.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

resultphi=np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[1],skiprows=1)
resultJet=np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[2],skiprows=1)
dNdphimin=min(resultdNdphi)
Ridge_Jetmin = min(resultdNdphi+resultJet)

# print(resultphi)
# print(dNdphimin)

cms_czyam = dNdphitable30-datamintable30
cms_result_czyam = resultdNdphi-dNdphimin
j = 4
cmsmse = 0
for i in range(len(dNdphitable30)):
        cmsmse += (cms_czyam[i]-cms_result_czyam[round(j)])**2
        j = 4+(i+1)*7.6699
cmsmse = (cmsmse/len(dNdphitable30))**0.5
print('2<pT<3')
print('CMS error = ', cmsmse)
wr.writerow(['2<pT<3',cmsmse])

errors +=  cmsmse


# plt.errorbar(deltaphitable30,dNdphitable30-datamintable30, yerr=(table30_error1,abs(table30_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# plt.errorbar(deltaphitable30,dNdphitable30+0.28, yerr=(table30_error1,abs(table30_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
plt.errorbar(deltaphitable30,dNdphitable30-datamintable30, yerr=(table30_error1,abs(table30_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)


# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet-Ridge_Jetmin, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=fr'$CMS,pp,7TeV,\, Deviation : {cmsmse:.3E}$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

plt.xlabel(r'$\Delta\phi$',size=50)
# plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
# plt.title(r'$2.0<p_T<3.0$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False, loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_Correlation_czyam_pt2-3.png')

fig.clear()

# fig = plt.figure()
# ax = plt.axes()

fig.set_size_inches(35, 16.534, forward=True)

deltaphitable32=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table32.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
dNdphitable32=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table32.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
datamintable32=min(dNdphitable32)
table32_error1=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table32.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
table32_error2=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table32.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

resultphi=np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[1],skiprows=1)
resultJet=np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[2],skiprows=1)
dNdphimin=min(resultdNdphi)
Ridge_Jetmin = min(resultdNdphi+resultJet)

# print(dNdphimin)

cms_czyam = dNdphitable32-datamintable32
cms_result_czyam = resultdNdphi-dNdphimin
j = 4
cmsmse = 0
for i in range(len(dNdphitable32)):
        cmsmse += (cms_czyam[i]-cms_result_czyam[round(j)])**2
        j = 4+(i+1)*7.6699
cmsmse = (cmsmse/len(dNdphitable32))**0.5
print('3<pT<4')
print('CMS error = ', cmsmse)
wr.writerow(['3<pT<4',cmsmse])

errors +=  cmsmse


# plt.errorbar(deltaphitable32,dNdphitable32-datamintable32, yerr=(table32_error1,abs(table32_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# plt.errorbar(deltaphitable32,dNdphitable32+0.09, yerr=(table32_error1,abs(table32_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
plt.errorbar(deltaphitable32,dNdphitable32-datamintable32, yerr=(table32_error1,abs(table32_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet-Ridge_Jetmin, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=fr'$CMS,pp,7TeV,\, Deviation : {cmsmse:.3E}$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

plt.xlabel(r'$\Delta\phi$',size=50)
# plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
# plt.title(r'$3.0<p_T<4.0$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False, loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_Correlation_czyam_pt3-4.png')


fig.clear()


resultphi=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[1],skiprows=1)
resultJet=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[2],skiprows=1)
dNdphimin=min(resultdNdphi)
Ridge_Jetmin = min(resultdNdphi+resultJet)


deltaphi_pt14=deltaphitable28
dNdphi_pt14=(dNdphitable28+1.27)+(dNdphitable30+0.28)+(dNdphitable32+0.09)
        # C_zyam 처리해야함 errorbar는 어떻게 할까
mindNdphi_pt14=min(dNdphi_pt14)

cms_czyam = dNdphi_pt14-mindNdphi_pt14
cms_result_czyam = resultdNdphi-dNdphimin
j = 4
cmsmse = 0
for i in range(len(deltaphi_pt14)):
        cmsmse += (cms_czyam[i]-cms_result_czyam[round(j)])**2
        j = 4+(i+1)*7.6699
cmsmse = (cmsmse/len(deltaphi_pt14))**0.5
print('1<pT<4')
print('CMS error = ', cmsmse)
wr.writerow(['1<pT<4',cmsmse])

errors +=  cmsmse

# plt.errorbar(deltaphi_pt14,dNdphi_pt14-mindNdphi_pt14, yerr=((table28_error1**2+table30_error1**2+table32_error1**2)**0.5,(table28_error2**2+table30_error2**2+table32_error2**2)**0.5), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
plt.errorbar(deltaphi_pt14,dNdphi_pt14-mindNdphi_pt14, yerr=((table28_error1**2+table30_error1**2+table32_error1**2)**0.5,(table28_error2**2+table30_error2**2+table32_error2**2)**0.5), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet-Ridge_Jetmin, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=fr'$CMS,pp,7TeV,\, Deviation : {cmsmse:.3E}$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

plt.xlabel(r'$\Delta\phi$',size=50)
# plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
# plt.title(r'$1.0<p_T<4.0$', fontsize = 60)
plt.minorticks_on()
# plt.yscale('log')
# ax.axis([-0.1,3.1,0,0.1])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
# plt.legend(fontsize=45,framealpha=False, loc='upper left')
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_Correlation_czyam_pt1-4.png')

fig.clear()



pt=np.loadtxt('pTdis.csv',delimiter=',',usecols=[0],skiprows=1)
aliceRidgedis=np.loadtxt('pTdis.csv',delimiter=',',usecols=[1],skiprows=1)
cmsRidgedis=np.loadtxt('pTdis.csv',delimiter=',',usecols=[2],skiprows=1)
alicejetdis=np.loadtxt('pTdis.csv',delimiter=',',usecols=[3],skiprows=1)
# cmsjetdis=np.loadtxt('pTdis.csv',delimiter=',',usecols=[4],skiprows=1)



CMS7TeVpt=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[0],skiprows=14,max_rows=4)
CMS7TeVptdis=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[1],skiprows=14,max_rows=4)
CMS7TeVptdis_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[2],skiprows=14,max_rows=4)
CMS7TeVptdis_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[3],skiprows=14,max_rows=4)
CMS7TeVptdis_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[4],skiprows=14,max_rows=4)
CMS7TeVptdis_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table34.csv',delimiter=',',usecols=[5],skiprows=14,max_rows=4)

CMS7TeVptdis_error1 = pow(CMS7TeVptdis_error1_stat*CMS7TeVptdis_error1_stat+CMS7TeVptdis_error1_sys*CMS7TeVptdis_error1_sys,0.5)
CMS7TeVptdis_error2 = pow(CMS7TeVptdis_error2_stat*CMS7TeVptdis_error2_stat+CMS7TeVptdis_error2_sys*CMS7TeVptdis_error2_sys,0.5)

cmsmse = 0

print('Y^Ridge')
for i in range(len(CMS7TeVpt)):
        # print(CMS7TeVpt[i],pt[i])
        cmsmse += (CMS7TeVptdis[i]-cmsRidgedis[i])**2
cmsmse = (cmsmse/4)**0.5

print('CMS error = ', cmsmse)

wr.writerow(['Y^Ridge',cmsmse])

# errors += alicemse + cmsmse

# print('ALICE error = ',alicemse)
# print('CMS error = ',cmsmse)

# wr.writerow(['Y^Ridge','ALICE : ',alicemse])
# wr.writerow(['Y^Ridge','CMS : ',cmsmse])

# yerr=(Alice13TeVptdis_error1,abs(Alice13TeVptdis_error2))
# plt.errorbar(Alice13TeVpt, Alice13TeVptdis, yerr=(abs(Alice13TeVptdis_error2),Alice13TeVptdis_error1), color="blue",zorder=1,markersize=35,marker='v',linestyle=' ',linewidth=3,label=r'$result,\, ALICE, \,1.6<\vert\Delta\eta\vert<1.8$',capsize=10)
plt.errorbar(CMS7TeVpt, CMS7TeVptdis, yerr=(abs(CMS7TeVptdis_error2),CMS7TeVptdis_error1), color="black",zorder=1,markersize=35,marker='^',linestyle=' ',linewidth=3,label=r'$result,\, CMS, \,2.0<\vert\Delta\eta\vert<4.0$',capsize=10)

# plt.plot(pt,aliceRidgedis, color="blue",linewidth=5,linestyle = '-',label=r'$pp,Ridge,1.6<\vert\Delta\eta\vert<1.8$')
# plt.scatter(pt,aliceRidgedis, color="skyblue",s=1000,linewidths=6,facecolors='none',marker='D',label=fr'$pp,Ridge,1.6<\vert\Delta\eta\vert<1.8,\, Deviation : {alicemse:.3E}$',zorder=2)
plt.scatter(pt,cmsRidgedis, color="crimson",s=1000,linewidths=6,facecolors='none',marker='s',label=fr'$pp,Ridge,2.0<\vert\Delta\eta\vert<4.0,\, Deviation : {cmsmse:.3E}$',zorder=2)

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

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
plt.tight_layout()

fig.savefig('pTdis.png')

fig.clear()


pt=np.loadtxt('dNdetaptdpt.csv',delimiter=',',usecols=[0],skiprows=1)
Ridgedis=np.loadtxt('dNdetaptdpt.csv',delimiter=',',usecols=[1],skiprows=1)
jetdis=np.loadtxt('dNdetaptdpt.csv',delimiter=',',usecols=[2],skiprows=1)


CMS7TeVpt=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/7TeV/HEPData-ins855299-v1-csv/Table4.csv',delimiter=',',usecols=[0],skiprows=13,max_rows=34)
CMS7TeVptdis=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/7TeV/HEPData-ins855299-v1-csv/Table4.csv',delimiter=',',usecols=[3],skiprows=13,max_rows=34)
CMS7TeVptdis_error1=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/7TeV/HEPData-ins855299-v1-csv/Table4.csv',delimiter=',',usecols=[4],skiprows=13,max_rows=34)
CMS7TeVptdis_error2=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/7TeV/HEPData-ins855299-v1-csv/Table4.csv',delimiter=',',usecols=[5],skiprows=13,max_rows=34)

plt.errorbar(CMS7TeVpt,CMS7TeVptdis, yerr=(CMS7TeVptdis_error1,abs(CMS7TeVptdis_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
plt.xlabel(r'$p_{t}\quad(GeV)$',size=50)
plt.ylabel(r'$(1/N_{trig})dN_{ch}/p_{t}dp_{t}{}\quad(GeV^{-2})$',size=50)
# plt.plot(pt,jetdis,color="blue", linestyle = '-.', linewidth=7,label=r'$pp,Jet$')
# plt.plot(pt,Ridgedis+0.632*jetdis, color="black",linewidth=7,label=r'$p+p,Jet+Ridge$')
plt.plot(pt,Ridgedis, color="red",linewidth=4,linestyle = '-',label=r'$p+p,Ridge$')


plt.minorticks_on()
plt.yscale('log')
ax.axis([0,4,0.0004,20])

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
# plt.ticklabel_format(axis='both',style='plain',useOffset=False)


plt.grid(color='silver',linestyle=':',linewidth=3)
plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
plt.tight_layout()

fig.savefig('dNptdpt.png')