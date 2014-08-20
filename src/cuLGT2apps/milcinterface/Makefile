
# define 
CCDEFS = 
CUDEFS = $(CCDEFS)

# c++ compiler
CC = gcc
# nvcc compiler
NVCC = nvcc

#includes
MILCINCLUDES = -I/home/vogt/milc/milc_qcd-7.7.11/pure_gauge -I/home/vogt/milc/milc_qcd-7.7.11/include
INCLUDES =  -I../../cuLGT2

INCLUDES += $(MILCINCLUDES)

# libs for the linking process
LIBS = 



# flags for the NVCC compiler
CUFLAGS = --ptxas-options=-v -arch=sm_20 -use_fast_math -Xptxas -dlcm=cg
# flags for the CC compiler
CCFLAGS = -O3

CUFLAGS += $(CCFLAGS)

all: milcinterface

milcinterface: MILCcuLGT.a

MILCcuLGT.a: linkMILCcuLGT.o MILCcuLGT.o Chronotimer.o
	ar rcs MILCcuLGT.a linkMILCcuLGT.o MILCcuLGT.o Chronotimer.o

Chronotimer.o: ../../util/timer/Chronotimer.cc
	$(CC) -c $(CCFLAGS) ../../util/timer/Chronotimer.cc

MILCcuLGT.o: MILCcuLGT.cu
	$(NVCC) -dc $(CUDEFS) $(CUFLAGS) $(INCLUDES) MILCcuLGT.cu
	
linkMILCcuLGT.o: MILCcuLGT.o Chronotimer.o
	$(NVCC) -dlink $(CUDEFS) $(CUFLAGS) MILCcuLGT.o -o linkMILCcuLGT.o

clean: 
	rm -f *.o MILCcuLGT.a 
