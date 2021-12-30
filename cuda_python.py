import numpy as np
import cupy as cp
#numba는 단순히 compiler임.
#cupy는 numpy의 gpu버전이라고 생각하면 됨.

#example

#numpy
x_cpu = np.array([1,2,3])
cpu_norm = np.linalg.norm(x_cpu)
#cupy
x_gpu = cp.array([1,2,3,])
gpu_norm = cp.linalg.norm(x_gpu)

#시작할 때에는 Device 번호는 0이다.

#1번 device를 시작하고 싶을 경우는 다음을 써야 한다.
cp.cuda.Device(1).use()

x.device # x를 이용하고 있는 device 번호를 출력하는 명령어

with cp.cuda.Device(1):
    x_on_gpu1 = cp.array([1,2,3,4,5]) # switch to another GPU device

x_gpu = cp.asarray(x_cpu) # move the data to the current device.(array on the current device cpu -> gpu)
#cpu -> gpu    gpu0 -> gpu1 같은 경우도 가능.

#cupy.asarry() does not copy the input array if possible. So, if you put an array of the current device, it returns the tinput ubject itself.
# If we do copy the array in this situation, you can use cupy.array() with copy=True. Actually, cupy.asarray() is equivalent to cupy.array(arr,dtype, copy=False).


Returns = cp.get_array_modue(*args)
# This function is used to implement CPU/GPU generic code.
# If at least one of the arguments is a cupy.ndarray object, the cupy module is returned.

# Parameters
# args – Values to determine whether NumPy or CuPy should be used.
# Returns
# cupy or numpy is returned based on the types of the arguments.
# Return type
# module

cp.asnumpy(x) # -> returns a Numpy array(array on the host gpu -> cpu)
#배열의 덧셈 등은 배열들이 모두 gpu에 있거나 cpu에 있는 경우에만 가능함을 주의하자.


cp.arange(start, stop=None, step=1, dtype=None).reshape() #->또 다른 배열 만드는 법 x배열 만드는 데에 유용해 보임.
cp.empty(dimension, dtype=float,order='C') # 공배열 만들기.dimension : 만드는 배열 차원 설정, dtype : 배열을 int, float등으로 설정, order : 배열 타입을 C/Fortran등으로 설정

