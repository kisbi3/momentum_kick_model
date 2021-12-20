import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import csv
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit
from multiprocessing import Process
import os


mpl.rcParams["text.usetex"] = True

errors = 0

f = open('deviation.csv','w',newline='')
wr=csv.writer(f)

fig = plt.figure()
ax = plt.axes()
# fig2,ax2 = plt.figure()

# alice 논문에서는 near-side 정의를 -1.28<phi<1.28으로 두었다.

fig.set_size_inches(35, 16.534, forward=True)
# fig.set_size_inches(50, 30, forward=True)


def functions(numbers):
        proc = os.getpid()
        if numbers == 1:
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

                # print(dNdphimin)z
                alice_czyam = alice_dNdphitable27-alice_datamintable27
                alice_result_czyam = alice_resultdNdphi-alice_dNdphimin
                cms_czyam = cms_dNdphitable27-cms_datamintable27
                cms_result_czyam = cms_resultdNdphi-cms_dNdphimin
                dphi = resultphi[1]-resultphi[0]
                dphi = 0.198/dphi
                j0 = j = 12
                alicemse_1 = 0
                cmsmse_1 = 0
                # print(resultphi)
                for i in range(len(alice_dNdphitable27)):
                        # print(cms_deltaphitable27[i],alice_deltaphitable27[i],resultphi[round(j)])
                        # alicemse_1 += abs(alice_czyam[i]-alice_result_czyam[round(j)])**2
                        # print(alice_deltaphitable27[i], resultphi[round(j)], alice_czyam[i], alice_result_czyam[round(j)])
                        # print(cms_deltaphitable27[i], resultphi[round(j)], cms_czyam[i], cms_result_czyam[round(j)])
                        alicemse_1 += abs(alice_czyam[i]-alice_result_czyam[round(j)])/alice_result_czyam[round(j)]
                        cmsmse_1 += abs(cms_czyam[i]-cms_result_czyam[round(j)])/cms_result_czyam[round(j)]
                        j = j0+(i+1)*dphi
                alicemse_1 = (alicemse_1/len(alice_dNdphitable27))*100
                cmsmse_1 = (cmsmse_1/len(cms_dNdphitable27))*100
                print('1<pT<2')
                print('ALICE error = ', alicemse_1)
                print('CMS error = ', cmsmse_1)
                wr.writerow(['1<pT<2','ALICE : ',alicemse_1])
                wr.writerow(['1<pT<2','CMS : ',cmsmse_1])

                errors_1 = 0
                errors_1 += alicemse_1 + cmsmse_1

                # print(alicemse_1/len(alice_dNdphitable27))
                # print(cmsmse_1/len(cms_dNdphitable27))


                plt.errorbar(alice_deltaphitable27,alice_dNdphitable27-alice_datamintable27, yerr=(alice_table27_error1,abs(alice_table27_error2)), color="blue",markersize=20,marker='o',linestyle=' ',fillstyle='none',linewidth=5,label=r'$pp,13TeV \, ALICE$',capsize=10)
                plt.errorbar(cms_deltaphitable27,cms_dNdphitable27-cms_datamintable27, yerr=(cms_table27_error1,abs(cms_table27_error2)), color="magenta",markersize=25,marker='o',linestyle=' ',fillstyle='none',linewidth=5,label=r'$pp,13TeV \, CMS$',capsize=15)

                # msg1 = 'alice error = ', alicemse_1}
                # print(msg1)


                # ,fillstyle='none'
                # plt.scatter(deltaphi,dNdphi, color="black", s = 80)
                # plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
                plt.plot(resultphi, alice_resultdNdphi-alice_dNdphimin, color = 'blue', linewidth=7, linestyle = '-',label=fr'$result,\, ALICE$')

                # ,\, Error : {alicemse_1:.2f}

                # plt.text(-0.5, 0.0025,fr"$ALICE\quad Error : {alicemse_1:.3E}$Deviationze=50)
                # plt.text(-0.5, 0.,fr"$CMS\quad Error : {cmsmse_1:.3E}$", size=50)
                plt.plot(resultphi, cms_resultdNdphi-cms_dNdphimin, color = 'magenta', linewidth=7, linestyle = '--',label=fr'$result,\, CMS$')

                # ,\, Error : {cmsmse_1:.2f}$

                # plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

                # plt.xlabel(r'$\Delta\phi$',size=70)

                plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=70)

                # plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 75)
                # plt.title(r'$1.0<p_T<2.0$', fontsize = 75)
                # plt.text(-0.25, 0.0, r'$1.0<p_T<2.0$', fontsize = 75)

                plt.minorticks_on()
                # plt.yscale('log')
                # ax.axis([-0.1,3.1,0,0.1])

                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

                plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top = 'true')
                plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top = 'true')
                # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


                plt.grid(color='silver',linestyle=':',linewidth=3)
                # plt.legend(fontsize=45,framealpha=False,loc='upper left')
                # plt.legend(fontsize=45,framealpha=False,bbox_to_anchor=(1,0.5))
                # plt.legend(fontsize=45)
                plt.tight_layout()

                fig.savefig('1D_PhiCorrelation_pt1-2.png')
                # plt.close()

                fig.clear()

# fig = plt.figure()
# ax = plt.axes()

# fig.set_size_inches(35, 16.534, forward=True)



        if numbers == 2:
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
                dphi = resultphi[1]-resultphi[0]
                dphi = 0.198/dphi
                j0 = j = 12
                alicemse_2 = 0
                cmsmse_2 = 0
                for i in range(len(alice_dNdphitable29)):
                        # print(cms_deltaphitable27[i],alice_deltaphitable27[i],resultphi[round(j)])
                        # alicemse_2 += abs(alice_czyam[i]-alice_result_czyam[round(j)])**2
                        # print(alice_deltaphitable29[i], resultphi[round(j)], alice_czyam[i], alice_result_czyam[round(j)])
                        # print(cms_deltaphitable29[i], resultphi[round(j)], cms_czyam[i], cms_result_czyam[round(j)])
                        alicemse_2 += abs(alice_czyam[i]-alice_result_czyam[round(j)])/alice_result_czyam[round(j)]
                        cmsmse_2 += abs(cms_czyam[i]-cms_result_czyam[round(j)])/cms_result_czyam[round(j)]
                        j = j0+(i+1)*dphi
                alicemse_2 = (alicemse_2/len(alice_dNdphitable29))*100
                cmsmse_2 = (cmsmse_2/len(cms_dNdphitable29))*100
                print('2<pT<3')
                print('ALICE error = ', alicemse_2)
                print('CMS error = ', cmsmse_2)
                wr.writerow(['2<pT<3','ALICE : ',alicemse_2])
                wr.writerow(['2<pT<3','CMS : ',cmsmse_2])

                errors_2 = 0
                errors_2 += alicemse_2 + cmsmse_2

                plt.errorbar(alice_deltaphitable29,alice_dNdphitable29-alice_datamintable29, yerr=(alice_table29_error1,abs(alice_table29_error2)), color="blue",markersize=20,marker='o',fillstyle='none',linestyle=' ',linewidth=5,label=r'$pp,13TeV \, ALICE$',capsize=15)
                plt.errorbar(cms_deltaphitable29,cms_dNdphitable29-cms_datamintable29, yerr=(cms_table29_error1,abs(cms_table29_error2)), color="magenta",markersize=25,marker='o',linestyle=' ',linewidth=5,fillstyle='none',label=r'$pp,13TeV \, CMS$',capsize=15)

                # ,fillstyle='none'
                # plt.scatter(deltaphi,dNdphi, color="black", s = 80)
                # plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
                plt.plot(resultphi, alice_resultdNdphi-alice_dNdphimin, color = 'blue', linewidth=7, linestyle = '-',label=fr'$result,\, ALICE$')

                # ,\, Error : {alicemse_2:.2f}

                plt.plot(resultphi, cms_resultdNdphi-cms_dNdphimin, color = 'magenta', linewidth=7, linestyle = '--',label=fr'$result,\, CMS$')

                # ,\, Error : {cmsmse_2:.2f}$

                # plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

                plt.xlabel(r'$\Delta\phi$',size=70)
                plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=70)

                # plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 75)
                plt.title(r'$2.0<p_T<3.0$', fontsize = 75)
                plt.minorticks_on()
                # plt.yscale('log')
                # ax.axis([-0.1,3.1,0,0.1])

                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

                plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
                plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
                # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


                plt.grid(color='silver',linestyle=':',linewidth=3)
                # plt.legend(fontsize=45,framealpha=False,loc='upper left')
                # plt.legend(fontsize=45)
                plt.tight_layout()

                fig.savefig('1D_PhiCorrelation_pt2-3.png')

                # fig.clear()

                # fig.set_size_inches(35, 16.534, forward=True)
                # plt.legend(fontsize=45,framealpha=False,loc='upper left')
                # fig.savefig('Legends')

                fig.clear()
                

                # fig = plt.figure()
                # ax = plt.axes()

        if numbers == 3:
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
                dphi = resultphi[1]-resultphi[0]
                dphi = 0.198/dphi
                j0 = j = 12
                alicemse_3 = 0
                cmsmse_3 = 0
                for i in range(len(alice_dNdphitable31)):
                        # print(cms_deltaphitable27[i],alice_deltaphitable27[i],resultphi[round(j)])
                        # alicemse_3 += abs(alice_czyam[i]-alice_result_czyam[round(j)])**2
                        # print(alice_deltaphitable27[i], resultphi[round(j)], alice_czyam[i], alice_result_czyam[round(j)])
                        # print(cms_deltaphitable31[i], resultphi[round(j)], cms_czyam[i], cms_result_czyam[round(j)])
                        alicemse_3 += abs(alice_czyam[i]-alice_result_czyam[round(j)])/alice_result_czyam[round(j)]
                        cmsmse_3 += abs(cms_czyam[i]-cms_result_czyam[round(j)])/cms_result_czyam[round(j)]
                        j = j0+(i+1)*dphi
                alicemse_3 = (alicemse_3/len(alice_dNdphitable31))*100
                cmsmse_3 = (cmsmse_3/len(cms_dNdphitable31))*100
                print('3<pT<4')
                print('ALICE error = ', alicemse_3)
                print('CMS error = ', cmsmse_3)

                errors_3 = 0
                errors_3 += alicemse_3 + cmsmse_3

                wr.writerow(['3<pT<4','ALICE : ',alicemse_3])
                wr.writerow(['3<pT<4','CMS : ',cmsmse_3])

                # print(dNdphimin)


                plt.errorbar(alice_deltaphitable31,alice_dNdphitable31-alice_datamintable31, yerr=(alice_table31_error1,abs(alice_table31_error2)), color="blue",markersize=20,marker='o',linestyle=' ',fillstyle='none',linewidth=5,label=r'$pp,13TeV \, ALICE$',capsize=15)
                plt.errorbar(cms_deltaphitable31,cms_dNdphitable31-cms_datamintable31, yerr=(cms_table31_error1,abs(cms_table31_error2)), color="magenta",markersize=25,marker='o',linestyle=' ', fillstyle='none',linewidth=5,label=r'$pp,13TeV \, CMS$',capsize=15)

                # ,fillstyle='none'
                # plt.scatter(deltaphi,dNdphi, color="black", s = 80)
                # plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
                plt.plot(resultphi, alice_resultdNdphi-alice_dNdphimin, color = 'blue', linewidth=7, linestyle = '-',label=fr'$result,\,ALICE$')

                # ,\, Error : {alicemse_3:.2f}

                plt.plot(resultphi, cms_resultdNdphi-cms_dNdphimin, color = 'magenta', linewidth=7, linestyle = '--',label=fr'$result,\,CMS$')

                # ,\, Error : {cmsmse_3:.2f}$

                # plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

                plt.xlabel(r'$\Delta\phi$',size=70)
                plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=70)

                # plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 75)
                plt.title(r'$3.0<p_T<4.0$', fontsize = 75)
                plt.minorticks_on()
                # plt.yscale('log')
                # ax.axis([-0.1,3.1,0,0.1])

                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

                plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
                plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
                # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


                plt.grid(color='silver',linestyle=':',linewidth=3)
                # plt.legend(fontsize=45,framealpha=False,loc='upper left')
                # plt.legend(fontsize=45)
                plt.tight_layout()

                fig.savefig('1D_PhiCorrelation_pt3-4.png')

                fig.clear()

        if numbers == 4:

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
                dphi = resultphi[1]-resultphi[0]
                dphi = 0.198/dphi
                j0 = j = 12
                cmsmse_4 = 0
                for i in range(len(alice_dNdphitable31)):
                        # print(cms_deltaphitable27[i],alice_deltaphitable27[i],resultphi[round(j)])
                        cmsmse_4 += abs(cms_czyam[i]-cms_result_czyam[round(j)])/cms_result_czyam[round(j)]
                        j = j0+(i+1)*dphi
                cmsmse_4 = (cmsmse_4/len(cms_dNdphitable31))*100
                print('1<pT<4')
                print('CMS error = ', cmsmse_4)
                wr.writerow(['1<pT<4','CMS : ',cmsmse_4])

                errors_4 = 0
                errors_4 += cmsmse_4

                plt.errorbar(deltaphi_pt14,dNdphi_pt14-mindNdphi_pt14, yerr=((cms_table27_error1**2+cms_table29_error1**2+cms_table31_error1**2)**0.5,(cms_table27_error2**2+cms_table29_error2**2+cms_table31_error2**2)**0.5), color="magenta",markersize=25,fillstyle='none',marker='o',linestyle=' ',linewidth=5,label=r'$pp,13TeV \, CMS$',capsize=15)

                # ,fillstyle='none'
                # plt.scatter(deltaphi,dNdphi, color="black", s = 80)
                # plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
                plt.plot(resultphi, alice_resultdNdphi-alice_dNdphimin, color = 'blue', linewidth=7, linestyle = '-',label=fr'$result,\,ALICE$')
                plt.plot(resultphi, cms_resultdNdphi-cms_dNdphimin, color = 'magenta', linewidth=7, linestyle = '--',label=fr'$result,\,CMS$')

                # ,\, Error : {cmsmse_4:.2f}$

                # plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

                plt.xlabel(r'$\Delta\phi$',size=70)
                plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=70)

                # plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 75)
                plt.title(r'$1.0<p_T<4.0$', fontsize = 75)
                plt.minorticks_on()
                # plt.yscale('log')
                # ax.axis([-0.1,3.1,0,0.1])

                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

                plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
                plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
                # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


                plt.grid(color='silver',linestyle=':',linewidth=3)
                # plt.legend(fontsize=45,framealpha=False,loc='upper left')
                # plt.legend(fontsize=45)
                plt.tight_layout()

                fig.savefig('1D_PhiCorrelation_pt1-4.png')

                fig.clear()

        if numbers == 5:
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

                Alice13TeVpt_low = np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[1],skiprows=12,max_rows=5)
                Alice13TeVpt_high = np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge}.csv',delimiter=',',usecols=[2],skiprows=12,max_rows=5)
                Alice13TeVdpt = Alice13TeVpt_high-Alice13TeVpt_low

                norm_Alice = 0
                for i in range(len(Alice13TeVdpt)):
                        norm_Alice += Alice13TeVptdis[i]*Alice13TeVdpt[i]
                        # print(Alice13TeVdpt[i])
                norm_Alice = 1/norm_Alice

                Alice13TeVptdis = Alice13TeVptdis*norm_Alice

                # print(Alice13TeVptdis_error1)
                # print(Alice13TeVptdis_error2)

                CMS13TeVpt               =np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[0],skiprows=14,max_rows=9)
                CMS13TeVptdis            =np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[1],skiprows=14,max_rows=9)
                CMS13TeVptdis_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[2],skiprows=14,max_rows=9)
                CMS13TeVptdis_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[3],skiprows=14,max_rows=9)
                CMS13TeVptdis_error1_sys =np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[4],skiprows=14,max_rows=9)
                CMS13TeVptdis_error2_sys =np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[5],skiprows=14,max_rows=9)

                CMS13TeVptdis_error1 = pow(CMS13TeVptdis_error1_stat*CMS13TeVptdis_error1_stat+CMS13TeVptdis_error1_sys*CMS13TeVptdis_error1_sys,0.5)
                CMS13TeVptdis_error2 = pow(CMS13TeVptdis_error2_stat*CMS13TeVptdis_error2_stat+CMS13TeVptdis_error2_sys*CMS13TeVptdis_error2_sys,0.5)

                CMS13TeVpt_low = np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[6],skiprows=14,max_rows=9)
                CMS13TeVpt_high = np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV/HEPData-ins1397173-v1-csv/Table33.csv',delimiter=',',usecols=[7],skiprows=14,max_rows=9)
                CMS13TeVdpt = CMS13TeVpt_high-CMS13TeVpt_low

                norm_CMS=0
                for i in range(len(CMS13TeVdpt)):
                        norm_CMS += CMS13TeVptdis[i]*CMS13TeVdpt[i]
                        # print(Alice13TeVdpt[i])
                norm_CMS = 1/norm_CMS

                CMS13TeVptdis = CMS13TeVptdis*norm_CMS

                alicemse_5=0
                cmsmse_5 = 0

                # print(pt[24], Alice13TeVpt[0], Alice13TeVptdis[0], aliceRidgedis[24])
                # print(pt[34], Alice13TeVpt[1], Alice13TeVptdis[1], aliceRidgedis[34])
                # print(pt[44], Alice13TeVpt[2], Alice13TeVptdis[2], aliceRidgedis[44])
                # print(pt[54], Alice13TeVpt[3], Alice13TeVptdis[3], aliceRidgedis[54])

                alicemse_5 = abs(Alice13TeVptdis[0]-aliceRidgedis[24])/aliceRidgedis[24]+abs(Alice13TeVptdis[1]-aliceRidgedis[34])/aliceRidgedis[34]+abs(Alice13TeVptdis[2]-aliceRidgedis[44])/aliceRidgedis[44]+abs(Alice13TeVptdis[3]-aliceRidgedis[54])/aliceRidgedis[54]
                alicemse_5 = (alicemse_5/4)*100

                print('Y^Ridge')
                # print(Alice13TeVpt)
                # print(pt)
                # print(CMS13TeVpt)
                dpt = pt[1]-pt[0]
                dpt = 10
                j0 = j = 4
                for i in range(len(CMS13TeVpt)-3):
                        # print(CMS13TeVpt[i],pt[i*dpt+j0])
                        cmsmse_5 += abs(CMS13TeVptdis[i]-cmsRidgedis[i*dpt+j0])/cmsRidgedis[i*dpt+j0]
                cmsmse_5 += abs(CMS13TeVptdis[7]-cmsRidgedis[68])/cmsRidgedis[68]+abs(CMS13TeVptdis[7]-cmsRidgedis[98])/cmsRidgedis[98]
                cmsmse_5 = (cmsmse_5/8)*100

                errors_5 = 0
                errors_5 += alicemse_5 + cmsmse_5

                print('ALICE error = ',alicemse_5)
                print('CMS error = ',cmsmse_5)

                wr.writerow(['Y^Ridge','ALICE : ',alicemse_5])
                wr.writerow(['Y^Ridge','CMS : ',cmsmse_5])

                # yerr=(Alice13TeVptdis_error1,abs(Alice13TeVptdis_error2))
                plt.errorbar(Alice13TeVpt, Alice13TeVptdis, yerr=(abs(Alice13TeVptdis_error2)*norm_Alice,Alice13TeVptdis_error1*norm_Alice), color="blue",zorder=1,markersize=35,marker='v',linestyle=' ',linewidth=3,label=r'$pp, \, 13TeV, \, ALICE$',capsize=10)
                # result,\, ALICE, \,1.6<\vert\Delta\eta\vert<1.8
                plt.errorbar(CMS13TeVpt, CMS13TeVptdis, yerr=(abs(CMS13TeVptdis_error2)*norm_CMS,CMS13TeVptdis_error1*norm_CMS), color="black",zorder=1,markersize=35,marker='^',linestyle=' ',linewidth=3,label=r'$pp, \, 13TeV, \, CMS$',capsize=10)
                # result,\, CMS, \,2.0<\vert\Delta\eta\vert<4.0

                # plt.plot(pt,aliceRidgedis, color="blue",linewidth=5,linestyle = '-',label=r'$pp,Ridge,1.6<\vert\Delta\eta\vert<1.8$')
                # plt.scatter(pt,aliceRidgedis, color="skyblue",s=1000,linewidths=6,facecolors='none',marker='D',label=fr'$pp,Ridge,1.6<\vert\Delta\eta\vert<1.8,\, Deviation : {alicemse_5:.3E}$',zorder=2)
                # plt.scatter(pt,cmsRidgedis, color="crimson",s=1000,linewidths=6,facecolors='none',marker='s',label=fr'$pp,Ridge,2.0<\vert\Delta\eta\vert<4.0$',zorder=2)

                # ,\, Error : {cmsmse_5:.3E}$


                plt.plot(pt,aliceRidgedis, color="skyblue",linewidth=6,linestyle = '-',label=fr'$result, \, ALICE$',zorder=2)
                # pp,Ridge,1.6<\vert\Delta\eta\vert<1.8
                # ,\, Error : {alicemse:.3E}

                plt.plot(pt,cmsRidgedis, color="crimson",linewidth=6, linestyle = '-',label=fr'$result, \, CMS$',zorder=2)

                # pp,Ridge,2.0<\vert\Delta\eta\vert<4.0
                # ,\, Error : {cmsmse_5:.3E}$


                # plt.plot(alicept,alicejetdis, color="blue",linewidth=5,linestyle = '--',label=r'$pp,Jet,1.6<\vert\Delta\eta\vert<1.8$')
                # plt.plot(cmspt,cmsRidgedis, color="black",linewidth=5,linestyle = '-',label=r'$pp,Ridge,2.0<\vert\Delta\eta\vert<4.0$')

                # x=np.arange(0.45,10,0.01)
                # y1=0.06782985130264102*np.exp(-(0.9341561431595556/(x**2))-0.984985854370021*x)
                # y2=0.07314538280113343*np.exp((-0.8078484337720946/(x**2))-0.8777989424828683*x)



                #######Fitting 해본 흔적들

                # def func(x,a,b,c,d,e):
                #         # return a*np.exp((b/(x**2))+c*x)+d
                #         # return a*np.exp(b/(x**2)+c*x+d/x)
                #         return a*np.exp(b/(x**2)+c*x)
                #         # return a*d**(b/(x**2))

                # # popt, pcov = curve_fit(f=func,xdata=CMS13TeVpt,ydata=CMS13TeVptdis)
                # # fitx=np.arange(0,10,0.01)
                # # plt.plot(fitx, func(fitx, *popt), color='black', linewidth=7, label='Python Fit')


                # popt, pcov = curve_fit(f=func,xdata=pt,ydata=cmsRidgedis)
                # fitx=np.arange(0.3,10,0.01)
                # plt.plot(fitx, func(fitx, *popt), color='black', linewidth=7, label='exponential, Python Fit')

                # popt, pcov = curve_fit(f=func,xdata=pt,ydata=aliceRidgedis)
                # fitx=np.arange(0.3,10,0.01)
                # # plt.plot(fitx, func(fitx, *popt), color='blue', linewidth=7, label='exponential, Python Fit')

                # fit1 = np.polyfit(pt,cmsRidgedis,17)
                # def nppolyfit(x):
                #         y=0
                #         for i in range(len(fit1)):
                #                 y += fit1[i]*(x**(17-i))
                #         return y

                # plt.plot(fitx, nppolyfit(fitx), color='red', linewidth=7, label='Polynomial Python Fit')




                # plt.plot(x,y1,color = 'skyblue', linewidth=7, linestyle = '-')
                # plt.plot(x,y2,color = 'red', linewidth=7, linestyle = '-')



                plt.xlabel(r'$p_{t}\quad(GeV)$',size=70)
                # plt.ylabel(r'$(1/N_{trig})dN_{ch}/p_{t}dp_{t}{}\quad(GeV^{-2})$',size=70)
                plt.ylabel(r'$Y^{Ridge}$',size=70)
                # plt.plot(pt,jetdis,color="blue", linestyle = '-.', linewidth=7,label=r'$pp,Jet$')
                # plt.plot(pt,Ridgedis+0.632*jetdis, color="black",linewidth=7,label=r'$p+p,Jet+Ridge$')

                plt.minorticks_on()
                # plt.yscale('log')
                # ax.axis([0,4,0.0004,20])

                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

                plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
                plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
                # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


                plt.grid(color='silver',linestyle=':',linewidth=3)
                plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
                plt.tight_layout()

                fig.savefig('pTdis.png')

                print(f'average : {errors_5/9}')


                fig.clear()


        if numbers == 6:
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

                plt.xlabel(r'$p_{T,min}^{Jet}$',size=70)
                plt.ylabel(r'$Y^{near}$',size=70)

                # plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 75)
                plt.title(r'$Y^{near}$', fontsize = 75)
                plt.minorticks_on()
                # plt.yscale('log')
                # ax.axis([-0.1,3.1,0,0.1])

                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

                plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
                plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
                # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


                plt.grid(color='silver',linestyle=':',linewidth=3)
                plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
                # plt.legend(fontsize=45)
                plt.tight_layout()

                fig.savefig('Y^near.png')

                fig.clear()

        if numbers == 7:

                Alice13TeV_jetphi_over10=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=53,max_rows=13)
                Alice13TeV_jetdNdphi_over10=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=53,max_rows=13)

                Alice13TeV_jetdNdphi_over10_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=53,max_rows=13)
                Alice13TeV_jetdNdphi_over10_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=53,max_rows=13)
                Alice13TeV_jetdNdphi_over10_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=53,max_rows=13)
                Alice13TeV_jetdNdphi_over10_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=53,max_rows=13)
                Alice13TeV_jetdNdphi_over10_min=min(Alice13TeV_jetdNdphi_over10)

                Alice13TeV_jetdNdphi_over10_error1 = (Alice13TeV_jetdNdphi_over10_error1_stat**2+Alice13TeV_jetdNdphi_over10_error1_sys**2)**0.5
                Alice13TeV_jetdNdphi_over10_error2 = (Alice13TeV_jetdNdphi_over10_error2_stat**2+Alice13TeV_jetdNdphi_over10_error2_sys**2)**0.5

                resultjet_phi=np.loadtxt('phiCorrelation_pt_jetcut_10.csv',delimiter=',',usecols=[0],skiprows=1)
                resultjet_dNdphi=np.loadtxt('phiCorrelation_pt_jetcut_10.csv',delimiter=',',usecols=[1],skiprows=1)
                resultjet_dNdphi_min=min(resultjet_dNdphi)

                alice_czyam = Alice13TeV_jetdNdphi_over10-Alice13TeV_jetdNdphi_over10_min
                alice_result_czyam = resultjet_dNdphi-resultjet_dNdphi_min
                dphi = resultphi[1]-resultphi[0]
                dphi = 0.198/dphi
                j0 = j = 12
                alicemse_7 = 0
                cmsmse_7 = 0
                for i in range(len(Alice13TeV_jetdNdphi_over10)):
                        # print(resultjet_phi[round(j)], Alice13TeV_jetphi_over10[i])
                        # alicemse_7 += abs(alice_czyam[i]-alice_result_czyam[round(j)])**2

                        # print(Alice13TeV_jetphi_over10[i], resultjet_phi[round(j)], alice_czyam[i], alice_result_czyam[round(j)])

                        alicemse_7 += abs(alice_czyam[i]-alice_result_czyam[round(j)])/alice_result_czyam[round(j)]
                        j = j0+(i+1)*dphi
                alicemse_7 = (alicemse_7/len(Alice13TeV_jetdNdphi_over10))*100
                print('pT,Jet>10')
                print('ALICE error = ', alicemse_7)
                wr.writerow(['pT,Jet>10','ALICE : ',alicemse_7])


                plt.errorbar(Alice13TeV_jetphi_over10,Alice13TeV_jetdNdphi_over10-Alice13TeV_jetdNdphi_over10_min, yerr=(Alice13TeV_jetdNdphi_over10_error1,Alice13TeV_jetdNdphi_over10_error2), color="blue",markersize=20,marker='o',linestyle=' ',linewidth=5,label=r'$pp,13TeV \, ALICE$',capsize=15)
                # ,fillstyle='none'
                # plt.scatter(deltaphi,dNdphi, color="black", s = 80)
                # plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
                plt.plot(resultjet_phi, resultjet_dNdphi-resultjet_dNdphi_min, color = 'blue', linewidth=7, linestyle = '-',label=r'$result$')
                # plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

                plt.xlabel(r'$\Delta\phi$',size=70)
                plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=70)

                # plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 75)
                plt.title(r'$p_{T,Jet}>10 GeV/c$', fontsize = 75)
                plt.minorticks_on()
                # plt.yscale('log')
                # ax.axis([-0.1,3.1,0,0.1])

                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

                plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
                plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
                # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


                plt.grid(color='silver',linestyle=':',linewidth=3)
                plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
                # plt.legend(fontsize=45)
                plt.tight_layout()

                fig.savefig('1D_PhiCorrelation_pt_jetcut_10.png')

                fig.clear()

        if numbers == 8:

                Alice13TeV_jetphi_over20=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[0],skiprows=91,max_rows=13)
                Alice13TeV_jetdNdphi_over20=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[3],skiprows=91,max_rows=13)

                Alice13TeV_jetdNdphi_over20_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[4],skiprows=91,max_rows=13)
                Alice13TeV_jetdNdphi_over20_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[5],skiprows=91,max_rows=13)
                Alice13TeV_jetdNdphi_over20_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[6],skiprows=91,max_rows=13)
                Alice13TeV_jetdNdphi_over20_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/JETBIASED1-NTRIGDN-DPHI.csv',delimiter=',',usecols=[7],skiprows=91,max_rows=13)
                Alice13TeV_jetdNdphi_over20_min=min(Alice13TeV_jetdNdphi_over20)

                Alice13TeV_jetdNdphi_over20_error1 = (Alice13TeV_jetdNdphi_over20_error1_stat**2+Alice13TeV_jetdNdphi_over20_error1_sys**2)**0.5
                Alice13TeV_jetdNdphi_over20_error2 = (Alice13TeV_jetdNdphi_over20_error2_stat**2+Alice13TeV_jetdNdphi_over20_error2_sys**2)**0.5

                resultjet_phi=np.loadtxt('phiCorrelation_pt_jetcut_20.csv',delimiter=',',usecols=[0],skiprows=1)
                resultjet_dNdphi=np.loadtxt('phiCorrelation_pt_jetcut_20.csv',delimiter=',',usecols=[1],skiprows=1)
                resultjet_dNdphi_min=min(resultjet_dNdphi)

                alice_czyam = Alice13TeV_jetdNdphi_over20-Alice13TeV_jetdNdphi_over20_min
                alice_result_czyam = resultjet_dNdphi-resultjet_dNdphi_min
                dphi = resultphi[1]-resultphi[0]
                dphi = 0.198/dphi
                j0 = j = 12
                alicemse_8 = 0
                cmsmse_8 = 0
                for i in range(len(Alice13TeV_jetdNdphi_over20)):
                        # print(resultjet_phi[round(j)], Alice13TeV_jetphi_over10[i])
                        # alicemse_8 += abs(alice_czyam[i]-alice_result_czyam[round(j)])**2
                        alicemse_8 += abs(alice_czyam[i]-alice_result_czyam[round(j)])/alice_result_czyam[round(j)]
                        j = j0+(i+1)*dphi
                alicemse_8 = (alicemse_8/len(Alice13TeV_jetdNdphi_over20))*100
                print('pT,Jet>20')
                print('ALICE error = ', alicemse_8)
                wr.writerow(['pT,Jet>20','ALICE : ',alicemse_8])

                plt.errorbar(Alice13TeV_jetphi_over20,Alice13TeV_jetdNdphi_over20-Alice13TeV_jetdNdphi_over20_min, yerr=(Alice13TeV_jetdNdphi_over20_error1,Alice13TeV_jetdNdphi_over20_error2), color="blue",markersize=20,marker='o',linestyle=' ',linewidth=5,label=r'$pp,13TeV \, ALICE$',capsize=15)
                # ,fillstyle='none'
                # plt.scatter(deltaphi,dNdphi, color="black", s = 80)
                # plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
                plt.plot(resultjet_phi, resultjet_dNdphi-resultjet_dNdphi_min, color = 'blue', linewidth=7, linestyle = '-',label=r'$result$')
                # plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

                plt.xlabel(r'$\Delta\phi$',size=70)
                plt.ylabel(r'$\frac{1}{N_{trig}}\frac{dN^{pair}}{d\Delta\phi}-C_{ZYAM}$',size=70)

                # plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 75)
                plt.title(r'$p_{T,Jet}>20 GeV/c$', fontsize = 75)
                plt.minorticks_on()
                # plt.yscale('log')
                # ax.axis([-0.1,3.1,0,0.1])

                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

                plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
                plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)
                # plt.ticklabel_format(axis='both',style='plain',useOffset=False)


                plt.grid(color='silver',linestyle=':',linewidth=3)
                plt.legend(fontsize=45,framealpha=False, bbox_to_anchor=(1, 1))
                # plt.legend(fontsize=45)
                plt.tight_layout()

                fig.savefig('1D_PhiCorrelation_pt_jetcut_20.png')

                fig.clear()

        if numbers == 9:
                Alice13TeV_Y_ridge_jet_pt=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge,JET}.csv',delimiter=',',usecols=[0],skiprows=12,max_rows=7)
                Alice13TeV_Y_ridge_jet=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge,JET}.csv',delimiter=',',usecols=[1],skiprows=12,max_rows=7)

                Alice13TeV_Y_ridge_jet_error1_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge,JET}.csv',delimiter=',',usecols=[2],skiprows=12,max_rows=7)
                Alice13TeV_Y_ridge_jet_error2_stat=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge,JET}.csv',delimiter=',',usecols=[3],skiprows=12,max_rows=7)
                Alice13TeV_Y_ridge_jet_error1_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge,JET}.csv',delimiter=',',usecols=[4],skiprows=12,max_rows=7)
                Alice13TeV_Y_ridge_jet_error2_sys=np.loadtxt('/home/jaesung/Desktop/Dropbox/Code/WongCode/13TeV-Alice/HEPData-ins1840098-v1-csv/Y^mathrm{ridge,JET}.csv',delimiter=',',usecols=[5],skiprows=12,max_rows=7)
                Alice13TeV_Y_ridge_jet_min=min(Alice13TeV_Y_ridge_jet)

                Alice13TeV_Y_ridge_jet_error1 = (Alice13TeV_Y_ridge_jet_error1_stat**2+Alice13TeV_Y_ridge_jet_error1_sys**2)**0.5
                Alice13TeV_Y_ridge_jet_error2 = (Alice13TeV_Y_ridge_jet_error2_stat**2+Alice13TeV_Y_ridge_jet_error2_sys**2)**0.5

                result_Yridge_jetpt=np.loadtxt('pTdis_jetcut.csv',delimiter=',',usecols=[0],skiprows=1)
                result_Yridge_jet=np.loadtxt('pTdis_jetcut.csv',delimiter=',',usecols=[1],skiprows=1)
                result_Yridge_jet_min=min(result_Yridge_jet)


                plt.errorbar(Alice13TeV_Y_ridge_jet_pt,Alice13TeV_Y_ridge_jet-Alice13TeV_Y_ridge_jet_min, yerr=(Alice13TeV_Y_ridge_jet_error1,Alice13TeV_Y_ridge_jet_error2), color="blue",markersize=20,marker='o',linestyle=' ',linewidth=5,label=r'$pp,13TeV \, ALICE$',capsize=15)
                # ,fillstyle='none'
                # plt.scatter(deltaphi,dNdphi, color="black", s = 80)
                # plt.plot(deltaphi,dNdphi-datamin, color="blue", linestyle = '', linewidth = 13)
                plt.plot(result_Yridge_jetpt, result_Yridge_jet, color = 'blue', linewidth=7, linestyle = '-',label=r'$result$')
                # plt.plot(resultphi, resultdNdphi, color = 'red', linewidth=7, linestyle = '-',label=r'$result$')

                plt.xlabel(r'$p_{T,min}^{Jet}$',size=70)
                # plt.xlabel(r'$q(GeV)$',size=70)
                plt.ylabel(r'$Y^{ridge}$',size=70)

                # plt.title(r'$1.0<p_T<2.0 \quad N_{trk}^{offline}\geq105,\quad C_{ZYAM}=1.27,\quad \Delta\phi_{ZYAM}=1.18$', fontsize = 75)
                plt.title(r'$Y^{ridge}$', fontsize = 75)
                plt.minorticks_on()
                # plt.yscale('log')
                ax.axis([-0.1,41,-0.1,0.45])

                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))

                plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45)
                plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45)


                # plt.ticklabel_format(axis='both',style='plain',useOffset=False)
                # ax2 = ax.twinx()
                # q = 0.05*result_Yridge_jetpt+0.5
                # qy = np.zeros(len(result_Yridge_jet))-5
                # line2 = ax2.plot(q,qy)
                # ax2.set_xlabel('q')
                # ax2.spines.bottom.set_position(("axes", -0.15))
                # ax2.xaxis.set_label_position('bottom')
                # ax2.xaxis.set_ticks_position('bottom')
                # # ax2.axis([0,9,-1,1])
                # ax2.get_xlim(9)

                plt.grid(color='silver',linestyle=':',linewidth=3)
                # plt.legend(fontsize=45,framealpha=False, loc='upper left')
                # ax2.legend()
                # plt.legend(fontsize=45)
                plt.tight_layout()

                fig.savefig('Y^Ridge_Jet.png')


if __name__ == '__main__':
        procs = []

        for i in range(1,11):
                proc = Process(target=functions, args=(i,))
                procs.append(proc)
                proc.start()
        
        for proc in procs:
                proc.join()

