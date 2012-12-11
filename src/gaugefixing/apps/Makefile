################################################################################
#  Copyright 2012 Mario Schroeck, Hannes Vogt
#
#  This file is part of cuLGT.
#
#  cuLGT is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  any later version.
#
#  cuLGT is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
#
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

# double precision?
ifeq ($(PREC),DP)
DPREC = -DDOUBLEPRECISION
endif

# do the updates in DP even if links are stored in SP
ifeq ($(DPUPDATES),true)
DPREC += -DUSE_DP_ORUPDATE -DUSE_DP_MICROUPDATE -DUSE_DP_SAUPDATE # TODO add the flags for SA,SR
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

$(APP): $(CUOBJ) Chronotimer.o GlobalConstants.o
	$(NVCC) -o $@ $(CUOBJ) Chronotimer.o GlobalConstants.o $(LIBS)

%.o: %.cu
	$(NVCC) -c $(CUFLAGS) $< $(CUDEFS) 

%.o: %.cpp
	$(CC) -c $(CCFLAGS) $< $(CCDEFS)
#%.o: %.cc
#	$(CC) -c $(CCFLAGS) $< $(CCDEFS)

Chronotimer.o: ../../util/timer/Chronotimer.cc
	$(CC) -c $(CCFLAGS) ../../util/timer/Chronotimer.cc

GlobalConstants.o: ../GlobalConstants.cu
	$(NVCC) -c $(CUDEFS) $(CUFLAGS) ../GlobalConstants.cu

clean: 
	rm -f *.o $(APP)


