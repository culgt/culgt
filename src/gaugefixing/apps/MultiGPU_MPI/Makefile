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

# double precision?
ifeq ($(PREC),DP)
DPREC = -DDOUBLEPRECISION
endif

# do the updates in DP even if links are stored in SP
ifeq ($(DPUPDATES),true)
DPREC += -DUSE_DP_ORUPDATE -DUSE_DP_MICROUPDATE -DUSE_DP_SAUPDATE # TODO add the flags for SA,SR
endif

# the lattice size and ThreadBlock size
CCDEFS = $(DPREC) -D_T_=$(T) -D_X_=$(X) -D_NSB_=$(NSB)
CUDEFS = $(CCDEFS)
# mpi compiler
MPICC = mpicxx
# nvcc compiler
NVCC = nvcc
# the objects to be compiled with MPICC
CCOBJ = main.o MultiGPU_MPI_AlgorithmOptions.o
# the objects to be compiled with NVCC
CUOBJ = Chronotimer.o MultiGPU_MPI_LandauKernelsSU3.o MultiGPU_MPI_Reduce.o 
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

Chronotimer.o: ../../../util/timer/Chronotimer.cc
	$(NVCC) -c ../../../util/timer/Chronotimer.cc

GlobalConstants.o: ../../GlobalConstants.cu
	$(NVCC) -c $(CUDEFS) $(CUFLAGS) ../../GlobalConstants.cu

clean: 
	rm -f *.o $(PROG)
