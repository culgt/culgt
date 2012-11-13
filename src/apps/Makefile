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

# define 
CCDEFS = $(DPREC) -D_T_=$(T) -D_X_=$(X) -D_NSB_=$(NSB)
CUDEFS = $(CCDEFS) -DCUDA
# gcc compiler
CC = gcc
# nvcc compiler
NVCC = nvcc
# the objects to be compiled with NVCC
CUOBJ = $(APP).o
# libs for the linking process
LIBS = -lboost_program_options
# flags for the NVCC compiler
CUFLAGS = --ptxas-options=-v -arch=sm_20 -use_fast_math


$(APP): $(CUOBJ) Chronotimer.o
	$(NVCC) -o $@ $(CUOBJ) Chronotimer.o $(LIBS)

%.o: %.cu
	$(NVCC) -c $(CUFLAGS) $< $(CUDEFS) 

%.o: %.cpp
	$(CC) -c $< $(CCDEFS)

Chronotimer.o: ../util/timer/Chronotimer.cc
	$(CC) -c ../util/timer/Chronotimer.cc


clean: 
	rm *.o $(APP)


