################################################################################
#
# Makefile 
# ========
# examples:
#
# make (uses default values and default app.)
# make APP=LandauGaugeFixingSU3_4D (compiles and links LandauGaugeFixingSU3_4D.cu)
# make T=96 X=48 (lattice size)
# make PREC=DP (use double precision)
# make DPUPDATES=true (do the updates in DP even if links are stored in SP)
#
################################################################################

# default values of lattice size and ThreadBlock size
T=32
X=32
NSB=32

# default app
APP = CoulombGaugeFixingSU3_4D

# double precision?
ifeq ($(PREC),DP)
DPREC = -DDOUBLEPRECISION
endif

# do the updates in DP even if links are stored in SP
ifeq ($(DPUPDATES),true)
DPREC += -DUSE_DP_ORUPDATE # TODO add the flags for SA,SR,Micro
endif

# define 
CCDEFS = $(DPREC) -D_T_=$(T) -D_X_=$(X) -D_NSB_=$(NSB)
CUDEFS = $(CCDEFS)
# gcc compiler
CC = gcc
# nvcc compiler
NVCC = nvcc
# the objects to be compiled with NVCC
CUOBJ = $(APP).o
# libs for the linking process
LIBS = -lboost_program_options
# flags for the NVCC compiler
CUFLAGS = --ptxas-options=-v -arch=sm_20 -use_fast_math -Xptxas -dlcm=cg
# flags for the CC compiler
CCFLAGS = -O3

CUFLAGS += $(CCFLAGS)

$(APP): $(CUOBJ) Chronotimer.o
	$(NVCC) -o $@ $(CUOBJ) Chronotimer.o $(LIBS)

%.o: %.cu
	$(NVCC) -c $(CUFLAGS) $< $(CUDEFS) 

%.o: %.cpp
	$(CC) -c $< $(CCDEFS)

Chronotimer.o: ../util/timer/Chronotimer.cc
	$(CC) -c ../util/timer/Chronotimer.cc


clean: 
	rm -f *.o $(APP)


