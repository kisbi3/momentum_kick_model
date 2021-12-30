import numpy as np
import cupy as cp

# __global__ : cpu가 gpu에 할당하는 것
# __device__ : gpu에서 요청하여 다른 gpu가 실행하는 것.
# __host__ : cpu에서 요청하여 cpu에서 실행하는 것. -> 일반적으로 생략.

func = cp.RawKernel(r'''
    extern "C" __global__
    void grid(){
        
//        char x[100] = blockIdx.x;
//        char y[100] = blockDim.x;
//        char z[100] = threadIdx.x;
        int x = blockIdx.x;
        int y = blockIdx.y;
        int z = blockIdx.z;
        
        int x1 = blockDim.x;
        int y1 = blockDim.y;
        int z1 = blockDim.z;

        int x2 = threadIdx.x;
        int y2 = threadIdx.y;
        int z2 = threadIdx.z;

        int x3 = gridDim.x;
        int y3 = gridDim.y;
        int z3 = gridDim.z;
        printf("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3);
    }

''', 'grid')

func((100,101,102), (10,10,10), ())
print('complete')
# func((?1,?2), (?3,), ())
# ?3의 경우 Dim은 변하지 않는다. -> 고정된 상수이다. 하지만 이는 thread의 개수인 것 같다.
# blockDim -> 한 block에서 thread가 (x or y or z 축에서) 몇개가 잘려있는지 표현하는 것.
# ex : Dim은 2로 고정일 때 thread는 0, 1이 반복되어서 나옴
# ?1, ?2는 block의 xyz좌표임. 하나 더 추가해서 3차원으로 만들 수 있음.

# 따라서 위의 경우 blockDim은 x, y, z의 경우 모두 10으로 고정
# blockIdx, threadIdx의 경우에는 x, y, z의 경우 모두 0~9까지 나란함.
# thread의 최대 개수는 1024개로 추정. 따라서 총 1000개가 가장 유용할 것으로 추정.
# 이때 thread의 최대 개수는 thread.x .y .z를 모두 곱한 값임.
# gridDim : block 개수


# Maximum number of threads per multiprocessor:   1024
# Maximum number of thread per block:     1024
# Maximum sizes of each dimension of a block:     1024 x 1024 x 64
# Maximum sizes of each dimension of a grid:      2147483647 x 65535 x 65535
# => 2080 Ti    =>  1660 super와도 동일한듯.

# threadIdx.x : 해당 thread block내에서 thread의 index를 가리킴
# blockDim.x : thread block에 몇개의 thread가 들어있는지 그 size를 가리킴
# blockIdx.x : 해당 grid에서 몇번 째 thread block인지 그 index를 가리킴
# gridDim.x : gird의 size를 의미함. 즉 몇개의 thread block이 이 grid에 들어있는지를 가리킴


# https://newsight.tistory.com/135