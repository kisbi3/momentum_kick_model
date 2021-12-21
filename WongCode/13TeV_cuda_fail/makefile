# Options
LINK = --device-link
OUT = --output-file
THREAD = -lpthread
CUDA_COMPILE = --device-c

# Header files
# CUDA_HEADER = cuda_function cuda_integral
# HEADER = integral function

CUD_HEADER  = $(wildcard ./headerfiles/*.cu)
CUH_HEADER = $(wildcard ./headerfiles/*.cuh)
# CU_HEAD = $(addsuffix .o, $(basename $(CUD_HEADER)))
CU_HEAD = $(CUD_HEADER:.cu=.o)
CPP_HEADER = $(wildcard ./headerfiles/*.cpp)
HPP_HEADER = $(wildcard ./headerfiles/*.hpp)
# CP_HEAD = $(addsuffix .o, $(basename $(CPP_HEADER)))
CP_HEAD = $(CPP_HEADER:.cpp=.o)
TARGET = test

# all : $(TARGET).out

all : test.out

$(CP_HEAD) : $(HPP_HEADER) $(CPP_HEADER)
	g++ -c $(CPP_HEADER)
	mv $(notdir $(CPP_HEADER:.cpp=.o)) ./headerfiles


$(CU_HEAD) : $(CU_FILES) $(CUH_FILES)
	nvcc $(CUDA_COMPILE) $(CUD_HEADER)
	mv $(notdir $(CUD_HEADER:.cu=.o)) ./headerfiles

cuda_headers.o : $(CU_HEAD)
	nvcc $(LINK) $(CU_HEAD) $(OUT) ./headerfiles/cuda_headers.o

$(TARGET).o : $(TARGET).cpp
	g++ -c $(TARGET).cpp

$(TARGET).out : $(TARGET).o $(CP_HEAD) $(CU_HEAD) cuda_headers.o
	nvcc $(TARGET).o $(CP_HEAD) $(CU_HEAD) ./headerfiles/cuda_headers.o $(OUT) $(TARGET).out


clean :
	rm -f *.o
	rm -f ./headerfiles/*.o

run : 
	rm -f *.o
	rm -f ./headerfiles/*.o
	./test.out
