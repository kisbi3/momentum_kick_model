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

fig.set_size_inches(35, 16.534, forward=True)

# fRNk=np.loadtxt('FittingParameters.csv',delimiter=',',usecols=[1],skiprows=8)

deltaphitable26=np.loadtxt('HEPData-ins1397173-v1-Table_26.csv',delimiter=',',usecols=[0],skiprows=19, max_rows=11)
dNdphitable26=np.loadtxt('HEPData-ins1397173-v1-Table_26.csv',delimiter=',',usecols=[1],skiprows=19, max_rows=11)
datamintable26=min(dNdphitable26)
table26_error1=np.loadtxt('HEPData-ins1397173-v1-Table_26.csv',delimiter=',',usecols=[2],skiprows=19, max_rows=11)
table26_error2=np.loadtxt('HEPData-ins1397173-v1-Table_26.csv',delimiter=',',usecols=[3],skiprows=19, max_rows=11)

# resultphi=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[0],skiprows=119, max_rows=237)
# resultdNdphi=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[1],skiprows=119, max_rows=237)
resultphi=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[1],skiprows=1)
resultJet=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[2],skiprows=1)
dNdphimin=min(resultdNdphi)
Ridge_Jetmin = min(resultdNdphi+resultJet)

# plt.errorbar(deltaphitable26,dNdphitable26-datamintable26, yerr=(table26_error1,abs(table26_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# plt.errorbar(deltaphitable26,dNdphitable26+4.77, yerr=(table26_error1,abs(table26_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
plt.errorbar(deltaphitable26,dNdphitable26, yerr=(table26_error1,abs(table26_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# frnk=fRNk[0]
# plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet-Ridge_Jetmin, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
plt.plot(resultphi, resultdNdphi+resultJet, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
# plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

plt.xlabel(r'$\Delta\phi$',size=50)
plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=50)

# plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 60)
plt.title(r'$0.1<p_T<1.0$', fontsize = 60)
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

fig.savefig('1D_Correlation_czyam_pt01-1.png')

fig.clear()


deltaphitable28=np.loadtxt('HEPData-ins1397173-v1-Table_28.csv',delimiter=',',usecols=[0],skiprows=19, max_rows=11)
dNdphitable28=np.loadtxt('HEPData-ins1397173-v1-Table_28.csv',delimiter=',',usecols=[1],skiprows=19, max_rows=11)
datamintable28=min(dNdphitable28)
table28_error1=np.loadtxt('HEPData-ins1397173-v1-Table_28.csv',delimiter=',',usecols=[2],skiprows=19, max_rows=11)
table28_error2=np.loadtxt('HEPData-ins1397173-v1-Table_28.csv',delimiter=',',usecols=[3],skiprows=19, max_rows=11)

resultphi=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[1],skiprows=1)
resultJet=np.loadtxt('phiCorrelation_pt1-2.csv',delimiter=',',usecols=[2],skiprows=1)
dNdphimin=min(resultdNdphi)
Ridge_Jetmin = min(resultdNdphi+resultJet)

print(dNdphimin)


# plt.errorbar(deltaphitable28,dNdphitable28-datamintable28, yerr=(table28_error1,abs(table28_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# plt.errorbar(deltaphitable28,dNdphitable28+1.27, yerr=(table28_error1,abs(table28_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
plt.errorbar(deltaphitable28,dNdphitable28, yerr=(table28_error1,abs(table28_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet-Ridge_Jetmin, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
plt.plot(resultphi, resultdNdphi+resultJet, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
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

deltaphitable30=np.loadtxt('HEPData-ins1397173-v1-Table_30.csv',delimiter=',',usecols=[0],skiprows=19, max_rows=11)
dNdphitable30=np.loadtxt('HEPData-ins1397173-v1-Table_30.csv',delimiter=',',usecols=[1],skiprows=19, max_rows=11)
datamintable30=min(dNdphitable30)
table30_error1=np.loadtxt('HEPData-ins1397173-v1-Table_30.csv',delimiter=',',usecols=[2],skiprows=19, max_rows=11)
table30_error2=np.loadtxt('HEPData-ins1397173-v1-Table_30.csv',delimiter=',',usecols=[3],skiprows=19, max_rows=11)

resultphi=np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[1],skiprows=1)
resultJet=np.loadtxt('phiCorrelation_pt2-3.csv',delimiter=',',usecols=[2],skiprows=1)
dNdphimin=min(resultdNdphi)
Ridge_Jetmin = min(resultdNdphi+resultJet)

# print(resultphi)
print(dNdphimin)


# plt.errorbar(deltaphitable30,dNdphitable30-datamintable30, yerr=(table30_error1,abs(table30_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# plt.errorbar(deltaphitable30,dNdphitable30+0.28, yerr=(table30_error1,abs(table30_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
plt.errorbar(deltaphitable30,dNdphitable30, yerr=(table30_error1,abs(table30_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)


# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet-Ridge_Jetmin, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
plt.plot(resultphi, resultdNdphi+resultJet, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
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

deltaphitable32=np.loadtxt('HEPData-ins1397173-v1-Table_32.csv',delimiter=',',usecols=[0],skiprows=19, max_rows=11)
dNdphitable32=np.loadtxt('HEPData-ins1397173-v1-Table_32.csv',delimiter=',',usecols=[1],skiprows=19, max_rows=11)
datamintable32=min(dNdphitable32)
table32_error1=np.loadtxt('HEPData-ins1397173-v1-Table_32.csv',delimiter=',',usecols=[2],skiprows=19, max_rows=11)
table32_error2=np.loadtxt('HEPData-ins1397173-v1-Table_32.csv',delimiter=',',usecols=[3],skiprows=19, max_rows=11)

resultphi=np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[1],skiprows=1)
resultJet=np.loadtxt('phiCorrelation_pt3-4.csv',delimiter=',',usecols=[2],skiprows=1)
dNdphimin=min(resultdNdphi)
Ridge_Jetmin = min(resultdNdphi+resultJet)

print(dNdphimin)


# plt.errorbar(deltaphitable32,dNdphitable32-datamintable32, yerr=(table32_error1,abs(table32_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# plt.errorbar(deltaphitable32,dNdphitable32+0.09, yerr=(table32_error1,abs(table32_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
plt.errorbar(deltaphitable32,dNdphitable32, yerr=(table32_error1,abs(table32_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet-Ridge_Jetmin, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
plt.plot(resultphi, resultdNdphi+resultJet, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
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


resultphi=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[1],skiprows=1)
resultJet=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[2],skiprows=1)
dNdphimin=min(resultdNdphi)
Ridge_Jetmin = min(resultdNdphi+resultJet)


deltaphi_pt14=deltaphitable28
dNdphi_pt14=(dNdphitable28+1.27)+(dNdphitable30+0.28)+(dNdphitable32+0.09)
        # C_zyam 처리해야함 errorbar는 어떻게 할까
mindNdphi_pt14=min(dNdphi_pt14)

# plt.errorbar(deltaphi_pt14,dNdphi_pt14-mindNdphi_pt14, yerr=((table28_error1**2+table30_error1**2+table32_error1**2)**0.5,(table28_error2**2+table30_error2**2+table32_error2**2)**0.5), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
plt.errorbar(deltaphi_pt14,dNdphi_pt14, yerr=((table28_error1**2+table30_error1**2+table32_error1**2)**0.5,(table28_error2**2+table30_error2**2+table32_error2**2)**0.5), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,7TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
# plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
# plt.plot(resultphi, resultdNdphi+resultJet-Ridge_Jetmin, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$Ridge$')
plt.plot(resultphi, resultJet, color = 'blue', linewidth=7, linestyle = '-',label=r'$Jet$')
plt.plot(resultphi, resultdNdphi+resultJet, color = 'black', linewidth=7, linestyle = '-',label=r'$Ridge+Jet$')
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
# plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
# plt.legend(fontsize=45)
plt.tight_layout()

fig.savefig('1D_Correlation_czyam_pt1-4.png')

fig.clear()


