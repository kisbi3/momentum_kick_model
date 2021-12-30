import numpy as np
import cupy as cp
# import scipy.integrate as integrate
from scipy import integrate
import time

time_start = time.time()

md = 1.
Njet = 0.75
Tjet = 0.55
m = 0.13957018
sqrSnn = 200.
mp = 0.938272046
T = 0.5
a = 0.5

# frnkconst = Njet/(Tjet*(m+Tjet)*2*np.pi)


# def frnk(pt):
#     return frnkconst*np.exp(-pt/np.sqrt(md*md+pt*pt))


#numpy
def lightcone(pti, yi):
    yb = np.arccosh(sqrSnn/(2.*mp))
    squareroot = np.sqrt(m*m+pti*pti)
    yiabs = np.fabs(yi)
    
    return (squareroot/m)*np.exp(yiabs-yb)
    

def integralAridge(pti, yi):
    x = lightcone(pti,yi)
    # suqreroot = np.sqrt(m*m+pti*pti)

    if x>=1.:
        return 0.
    
    else:
        # return pti*np.power(1-x,a)*np.exp(-np.sqrt(m*m+pti*pti)/T)/np.sqrt(md*md+pti*pti)
        return pti*np.power(1-x,a)*np.exp(-np.sqrt(m*m+pti*pti)/T)/np.sqrt(md*md+pti*pti)
    


#cupy
def lightcone(pti, yi):
    yb = cp.arccosh(sqrSnn/(2.*mp))
    squareroot = cp.sqrt(m*m+pti*pti)
    yiabs = cp.absolute(yi)
    
    return (squareroot/m)*cp.exp(yiabs-yb)
    

def integralAridge(pti, yi):
    x = lightcone(pti,yi)
    # suqreroot = np.sqrt(m*m+pti*pti)

    if x>=1.:
        return 0.
    
    else:
        # return pti*np.power(1-x,a)*np.exp(-np.sqrt(m*m+pti*pti)/T)/np.sqrt(md*md+pti*pti)
        return pti*cp.power(1-x,a)*cp.exp(-cp.sqrt(m*m+pti*pti)/T)/cp.sqrt(md*md+pti*pti)

# lightcone2 = cp.ElementwiseKernel(
#     'T pti, T yi', #에러가 뜬다면 float -> float32로 바꾸어 보자.
#     'T lightcone2',
#     '''
#         yb = cp.arccosh(sqrSnn/(2.*mp));
#         squareroot = cp.sqrt(m*m+pti*pti);
#         yiabs = cp.absolute(yi);

#         lightcone2 = (squareroot/m)*cp.exp(yiabs-yb);
#     ''',
#     'lightcone2'
# )

# integralAridge2 = cp.ElementwiseKernel(
#     'T pti, T yi',
#     'T integralAridge2',
#     '''
#         x = lightcone2(pti, yi);
#         if x>=1,:
#     '''
# )

#ElementwiseKernel -> Rawkernel 로 바꿔보자.

lightcone = cp.RawKernel(r'''
    extern "C" __device__

    double lightcone(pti, yi){
        double sqrSnn = 200.;
        double mp = 0.938272046;
        double a = 0.5;
        double T = 0.5;
        double md = 1.;
        double m = 0.13957018;

        double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mpi

        double squareroot=sqrt(m*m+pti*pti);
        double yiabs = abs(yi);


        return (squareroot/m)*exp(yiabs-yb);   
    }


''', 'lightcone')

integralAridge2 = cp.RawKernel(r'''
//    extern "C" __global__
    extern "C"

//    double x = 0.;
//    double md = 1.;
//    double Njet = 0.75;
//    double Tjet = 0.55;
//    double m = 0.13957018;
//    double sqrSnn = 200.;
//    double mp = 0.938272046;
//    double T = 0.5;
//    double a = 0.5;

//    # void lightcone(double* pti, double* yi, double* x){
//    #     double sqrSnn = 200.;
//    #     double m = 0.13957018;
//    #     double mp = 0.938272046;
//
//    #     double yb = acosh(sqrSnn/(2.*mp));
//    #     double squareroot=sqrt(m*m+pti*pti);
//    #     double yiabs = std::fabs(yi);

//    #     x = (squareroot/m)*exp(yiabs-yb);
//    # }
    __device__ double lightcone(double* pti, double* yi){
        double sqrSnn = 200.;
        double mp = 0.938272046;
        double m = 0.13957018;

        double yb = acosh(sqrSnn/(2.*mp));    //mN=mbeam, mb = mpi

        double squareroot=sqrt(m*m+pti*pti);
        //double yiabs = fabs(yi);


        return (squareroot/m)*exp(yi-yb);   
    }
    __global__ void integralAridge(double* pti, double* yi, int n, double* sum){
        double sqrSnn = 200.;
        double mp = 0.938272046;
        double md = 1.;
        double m = 0.13957018;
        double T = 0.5;
        double a = 0.5;
        double x = 0.;

//        int index = blockIdx.x*blockDim.x*n+threadIdx.x;

//        int row = blockIdx.x*blockDim.x;
//        int column = threadIdx.x;

        int row = blockIdx.x;
//        int column = blockIdx.y;
        int column = threadIdx.x;


        double yb = acosh(sqrSnn/(2.*mp));
        double squareroot=sqrt(m*m+pti[row]*pti[row]);
        double yiabs = fabs(yi[column]);
//        x = (squareroot/m)*exp(yiabs-yb);

        x = lightcone(pti, yi);

//        # lightcone(pti, yi, x);
        
        
        if(x>=1.){
            sum[row*blockDim.x+column] = 0.;
//            sum[row][column] = 0.;
        }
        else{
            sum[row*blockDim.x+column] = pti[row]*pow(1-x,a)*exp(-sqrt(m*m+pti[row]*pti[row])/T)/sqrt(md*md+pti[row]*pti[row]);
//            sum[row][column] = pti[row]*pow(1-x,a)*exp(-sqrt(m*m+pti[row]*pti[row])/T)/sqrt(md*md+pti[row]*pti[row]);
        }
    }

''', 'integralAridge')

def function(x,y,z,a):
    return np.exp(x**2+y**2+z**2+a**2)


# total = integrate.nquad(function, [[1, 2],[1, 2],[1,2],[1,2]])
# total = integrate.quad( integrate.quad(function, 1, 2), 1, 2)
# total = integrate.nquad(integralAridge, [[0., 10.], [0., 10.]])
# total = integrate.quad(integralAridge, 0.15, 1.0)


# dist = integrate.nquad(integralAridge, [[5.,5.], [0.,10.]])



# dpti = (10.-0.)/200.
dpti = 0.01
pti = 0.
another = 0.
# pti = np.arange(0.,10.,0.01)
pti_gpu = cp.arange(0.,10.,0.0001)
# yi_gpu = cp.arange(0.,10.,0.0001)

# dist = [0.,0.]
# for i in range(int(10/0.01)):
#     dist = integrate.nquad(integralAridge, [[pti[i], pti[i]+dpti], [0.,10.]])
#     # pti += dpti
#     print(pti[i], dist[0])
#     another += 2*dist[0]


#cupy
# sum = 0.
# matrix = cp.zeros(len(pti_gpu), dtype=float)
# cal = cp.empty((len(pti_gpu),2), dtype=float, order='C')
# sum_gpu = cp.zeros((len(pti_gpu),len(yi_gpu)), dtype=float, order='C')

start = 0.
end = 10./100.
totalsum = 0.
# # pti_dist = cp.zeros((10,1000))
# # yi_dist = cp.zeros((10,1000))
# # sum_dist = cp.zeros((1000,1000))
# # print(pti_dist)
sum_dist = cp.zeros((100000, 1000), dtype=float, order='C')
for i in range(int(len(pti_gpu)/1000)):
    # print(pti_dist[i])
    # pti_dist = cp.arange(start, end, 0.001)
    # print('111')
    # print(pti_dist)
    yi_dist = cp.arange(start, end, 0.0001)
    # print(yi_dist)
    # sum_dist = cp.zeros((100000, 1000), dtype=float, order='C')
    integralAridge2((100000,),(1000,), (pti_gpu, yi_dist, 1000, sum_dist))
    # print(sum_dist)
    totalsum += cp.sum(sum_dist)
    start += 0.1
    end += 0.1

    # y initial에 대해서 적분 값이 이상함. 디버깅 필요.

# integralAridge2((len(pti_gpu),), (len(pti_gpu),), (pti_gpu, yi_gpu, len(pti_gpu), sum_gpu))
# integralAridge2((10000,10000), (1000,), (pti_gpu, yi_gpu, len(pti_gpu), sum_gpu))
# print(deviceQuery)
# integralAridge2((len(pti_gpu),len(yi_gpu)), (1,), (pti_gpu, yi_gpu, len(pti_gpu), sum_gpu))
#  ((grid), (block), (arguments))
# print(len(sum_gpu))
# print(len(pti_gpu))
# cp.empty(matrix, dtype=float)
# matrix = 2*integralAridge(pti_gpu,yi_gpu)
# matrix = cupy.zeros((len(pti_gpu),2), dtype = float, )

# print(matrix)

# n = len(pti_gpu)
# for i in range(len(pti_gpu)):
#     print(pti_gpu[i], sum_gpu[i])


# constant = cp.sum(sum_gpu)*4*cp.pi*0.001*0.001
constant = totalsum*4*cp.pi*0.0001*0.0001
Aridge = 1/constant
print(constant, Aridge)

time_end = time.time()
print("time : ",time_end-time_start, "sec")

# for i in range(len(pti_gpu)):
#     # cp.cuda.Device(i).use()
#     # sum = integrate.nquad(integralAridge, [[pti_gpu[i], pti_gpu[i]+0.01], [0.,10.]])
#     for j in range(len(yi_gpu)):
#         dist = 2*integralAridge(pti_gpu[i],yi_gpu[j])*0.01*0.01
#         sum += dist
#         # print(pti_gpu[i], yi_gpu[j], dist, sum)
#     print(pti_gpu[i], dist, sum)


# print(sum*4*cp.pi)


# print(another*2*np.pi)

# print(type(np.pi))
# print(type(total))
# print(type(total[0]))

# Aridge = total[0]*4*np.pi
# rapidity 는 대칭이므로 x2 해주어야 함. 

# print(Aridge)
# print(1/(Aridge))