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

deltaphitable25=np.loadtxt('HEPData-ins1397173-v1-Table_25.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
dNdphitable25=np.loadtxt('HEPData-ins1397173-v1-Table_25.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
datamintable25=min(dNdphitable25)
table25_error1=np.loadtxt('HEPData-ins1397173-v1-Table_25.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
table25_error2=np.loadtxt('HEPData-ins1397173-v1-Table_25.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

resultphi=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt01-1.csv',delimiter=',',usecols=[1],skiprows=1)
dNdphimin=min(resultdNdphi)

plt.errorbar(deltaphitable25,dNdphitable25-datamintable25, yerr=(table25_error1,abs(table25_error2)), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,13TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
# frnk=fRNk[0]
plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$result,f_R<N_k>=$')
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


deltaphitable27=np.loadtxt('HEPData-ins1397173-v1-Table_27.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
dNdphitable27=np.loadtxt('HEPData-ins1397173-v1-Table_27.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
datamintable27=min(dNdphitable27)
table27_error1=np.loadtxt('HEPData-ins1397173-v1-Table_27.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
table27_error2=np.loadtxt('HEPData-ins1397173-v1-Table_27.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

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

deltaphitable29=np.loadtxt('HEPData-ins1397173-v1-Table_29.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
dNdphitable29=np.loadtxt('HEPData-ins1397173-v1-Table_29.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
datamintable29=min(dNdphitable29)
table29_error1=np.loadtxt('HEPData-ins1397173-v1-Table_29.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
table29_error2=np.loadtxt('HEPData-ins1397173-v1-Table_29.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

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

deltaphitable31=np.loadtxt('HEPData-ins1397173-v1-Table_31.csv',delimiter=',',usecols=[0],skiprows=18, max_rows=13)
dNdphitable31=np.loadtxt('HEPData-ins1397173-v1-Table_31.csv',delimiter=',',usecols=[1],skiprows=18, max_rows=13)
datamintable31=min(dNdphitable31)
table31_error1=np.loadtxt('HEPData-ins1397173-v1-Table_31.csv',delimiter=',',usecols=[2],skiprows=18, max_rows=13)
table31_error2=np.loadtxt('HEPData-ins1397173-v1-Table_31.csv',delimiter=',',usecols=[3],skiprows=18, max_rows=13)

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


resultphi=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[0],skiprows=1)
resultdNdphi=np.loadtxt('phiCorrelation_pt1-4.csv',delimiter=',',usecols=[1],skiprows=1)
dNdphimin=min(resultdNdphi)


deltaphi_pt14=deltaphitable27
dNdphi_pt14=(dNdphitable27+1.27)+(dNdphitable29+0.28)+(dNdphitable31+0.09)
        # C_zyam 처리해야함 errorbar는 어떻게 할까
mindNdphi_pt14=min(dNdphi_pt14)

plt.errorbar(deltaphi_pt14,dNdphi_pt14-mindNdphi_pt14, yerr=((table27_error1**2+table29_error1**2+table31_error1**2)**0.5,(table27_error2**2+table29_error2**2+table31_error2**2)**0.5), color="blue",markersize=15,marker='o',linestyle=' ',linewidth=3,label=r'$pp,13TeV \, CMS$',capsize=10)
# ,fillstyle='none'
# plt.scatter(deltaphi,dNdphi, color="black", s = 80)
# plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
plt.plot(resultphi, resultdNdphi-dNdphimin, color = 'red', linewidth=7, linestyle = '-',label=r'$result,\quad C_{ZYAM}=0.00501599$')
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


