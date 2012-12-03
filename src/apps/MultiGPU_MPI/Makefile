################################################################################
#
# Makefile 
# compile app. simply by typing 'make' (with default values for T,X,NSB)
# with non default values type, e.g.,
# 'make T=48 X=96'
#
################################################################################

#default values of lattice size and ThreadBlock size
T=32
X=32
NSB=32

# name of the output program
PROG = MultiGPU_MPI_LandauGaugeFixingSU3_4D



# the lattice size and ThreadBlock size
CCDEFS = -D_T_=$(T) -D_X_=$(X) -D_NSB_=$(NSB)
CUDEFS = $(CCDEFS)
# mpi compiler
MPICC = mpicxx
# nvcc compiler
NVCC = nvcc
# the objects to be compiled with MPICC
CCOBJ = main.o MPI_ProcInfo.o
# the objects to be compiled with NVCC
CUOBJ = Chronotimer.o $(PROG).o MultiGPU_MPI_LandauKernelsSU3.o MultiGPU_MPI_Reduce.o 
# libs for the linking process
LIBS = -lboost_program_options -L/usr/local/cuda/lib64 -lcudart 
# flags for the NVCC compiler
CUFLAGS = --ptxas-options=-v -arch=sm_20 -use_fast_math -Xptxas -dlcm=cg
# cuda include path
INC = -I/usr/local/cuda/include/


$(PROG): $(CCOBJ) $(CUOBJ)
	$(MPICC) -o $@ $(CCOBJ) $(CUOBJ) $(LIBS)

%.o: %.cpp
	$(MPICC) $(INC) -c $< $(CCDEFS)

%.o: %.cu
	$(NVCC) -c $(CUFLAGS) $< $(CUDEFS)

Chronotimer.o: ../../util/timer/Chronotimer.cc
	$(NVCC) -c ../../util/timer/Chronotimer.cc

clean: 
	rm -f *.o $(PROG)