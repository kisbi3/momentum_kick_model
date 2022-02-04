#배열의 경우 https://yunmorning.tistory.com/23 참조
import ctypes
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from scipy.optimize import curve_fit
from pprint import pprint   #배열 출력 이쁘게 해주는 것


param_F=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV_atlas/result/ridge_parameters.csv',delimiter=',',usecols=[1],skiprows=2)
param_G=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV_atlas/result/ridge_parameters.csv',delimiter=',',usecols=[2],skiprows=2)
param_v=np.loadtxt('/home/jaesung/OneDrive/Code/WongCode/13TeV_atlas/result/ridge_parameters.csv',delimiter=',',usecols=[3],skiprows=2)
phi = np.linspace(-1.,1.,301)

path="./function.so"
c_function = ctypes.cdll.LoadLibrary(path)

Aridge = c_function.py_Aridge
Aridge.restype = ctypes.c_double

cal = Aridge()

print(cal)

Ridge = c_function.py_phiCorr
Ridge.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double)
Ridge.restype = ctypes.POINTER(ctypes.c_double*301) #phi에 대해서 300개의 data를 얻어냄

ptr_free = c_function.freeptr
ptr_free.argtype = ctypes.c_void_p
ptr_free.restype = None

#double ptf_st, double ptf_end, double etaf_st, double etaf_end, double etacms_st, double etacms_end, double etaatlas_st, double etaatlas_end, double phif_st, double phif_end, int n, int check2, double ptjetcut, double q



def ridge(phi, G, v):
    return G*(1+2*v*np.cos(2*phi))

def atlas_ridge():
    ridge_integrate = np.zeros(5)
    integrate_phi = np.arange(-1., 1., 0.001)
    for j in range(5):
        ridge_integrate[j] = np.trapz(ridge(integrate_phi, param_G[4-j], param_v[4-j])-min(ridge(integrate_phi, param_G[4-j], param_v[4-j])), x = integrate_phi)



def ols(Aridge):
    Ridgearray_C = Ridge(Aridge, .5, 5., 1.6, 1.8, 2., 4., 2., 5., -1., 1., 300, 1, 0., 0.2)

    Ridgearray_py = [x for x in Ridgearray_C.contents]
    ptr_free(Ridgearray_C)
    Ridgearray=np.array(Ridgearray_py)
    theor_integrate = np.zeros(5)

    #그냥 numpy로 적분하자.
    theor_integrate[0] = np.trapz(Ridgearray-min(Ridgearray), x = phi)
    theor_integrate[1] = np.trapz(Ridgearray-min(Ridgearray), x = phi)
    theor_integrate[2] = np.trapz(Ridgearray-min(Ridgearray), x = phi)
    theor_integrate[3] = np.trapz(Ridgearray-min(Ridgearray), x = phi)
    theor_integrate[4] = np.trapz(Ridgearray-min(Ridgearray), x = phi)
    

    ratio_095 = ridge_integrate[0]/theor_integrate[0]
    ratio_105 = ridge_integrate[1]/theor_integrate[1]
    ratio_115 = ridge_integrate[2]/theor_integrate[2]
    ratio_125 = ridge_integrate[3]/theor_integrate[3]
    ratio_135 = ridge_integrate[4]/theor_integrate[4]

    theory_norm_095 = (Ridgearray-min(Ridgearray))*ratio_095
    theory_norm_105 = (Ridgearray-min(Ridgearray))*ratio_105
    theory_norm_115 = (Ridgearray-min(Ridgearray))*ratio_115
    theory_norm_125 = (Ridgearray-min(Ridgearray))*ratio_125
    theory_norm_135 = (Ridgearray-min(Ridgearray))*ratio_135

    return theory_norm_095, theory_norm_105, theory_norm_115, theory_norm_125, theory_norm_135


def only_ridge(theory_norm_095, theory_norm_105, theory_norm_115, theory_norm_125, theory_norm_135):
    fig = plt.figure()
    ax = plt.axes()
    mpl.rcParams["text.usetex"] = True
    color_ridge = ['grey', 'blue', 'orange', 'green', 'black']
    fig.set_size_inches(25, 20, forward=True)

    for j in range(5):
        graph_ridge = ridge(phi, param_G[4-j], param_v[4-j])
        ax.plot(phi, graph_ridge-min(graph_ridge), linewidth=5, color = color_ridge[j], linestyle='--')


    # print(ratio_095, ratio_105, ratio_115, ratio_125, ratio_135)

    plt.plot(phi, theory_norm_095, color = color_ridge[0], linewidth=5, linestyle = '-',label=r'$90 \,\, \leq N_{ch}^{rec}<100$')
    plt.plot(phi, theory_norm_105, color = color_ridge[1], linewidth=5, linestyle = '-',label=r'$100 \leq N_{ch}^{rec}<110$')
    plt.plot(phi, theory_norm_115, color = color_ridge[2], linewidth=5, linestyle = '-',label=r'$110 \leq N_{ch}^{rec}<120$')
    plt.plot(phi, theory_norm_125, color = color_ridge[3], linewidth=5, linestyle = '-',label=r'$120 \leq N_{ch}^{rec}<130$')
    plt.plot(phi, theory_norm_135, color = color_ridge[4], linewidth=5, linestyle = '-',label=r'$130 \leq N_{ch}^{rec}$')
 

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

    fig.savefig('only_ridges.png')

    fig.clear()

theory_norm_095, theory_norm_105, theory_norm_115, theory_norm_125, theory_norm_135 = ols(cal)
only_ridge(theory_norm_095, theory_norm_105, theory_norm_115, theory_norm_125, theory_norm_135)

