#possible flags: -DDOUBLEPRECISION, -DOSX
#DOUBLEPRECISION: to use double precission otherwise single precission
#FILEOPPOSITEPREC: load/save the configuration files as opposite precission as the being used

coulomb:
	nvcc --ptxas-options=-v -arch=sm_12 ../util/timer/Chronotimer.cc CoulombGaugeFixingSU3_4D.cu -use_fast_math -o coulomb -D_T_=8 -D_X_=8 -lboost_program_options

landau:
	nvcc --ptxas-options=-v -arch=sm_12 ../util/timer/Chronotimer.cc LandauGaugeFixingSU3_4D.cu -use_fast_math -o landau -D_T_=8 -D_X_=8 -lboost_program_options
SRlandau:
	nvcc --ptxas-options=-v -arch=sm_20 ../util/timer/Chronotimer.cc SRLandauGaugeFixingSU3_4D.cu -use_fast_math -o SRlandau -D_T_=8 -D_X_=8 -lboost_program_options
MAG:
	nvcc --ptxas-options=-v -arch=sm_20 ../util/timer/Chronotimer.cc MAGaugeFixingSU3_4D.cu -use_fast_math -o MAG -D_T_=32 -D_X_=32 -lboost_program_options
MAG_U1xU1:
	nvcc --ptxas-options=-v -arch=sm_20 ../util/timer/Chronotimer.cc MA_U1xU1_GaugeFixingSU3_4D.cu -use_fast_math -o MAG_U1xU1 -D_T_=64 -D_X_=20 -DDOUBLEPRECISION -lboost_program_options

multiGPU:
	nvcc --ptxas-options=-v -arch=sm_12 ../util/timer/Chronotimer.cc multiGPU_GaugeFixingSU3_4D.cu -use_fast_math -o multiGPU -DOSX -D_T_=8 -D_X_=8 -lboost_program_options

clean:
	rm coulomb landau SRlandau MAG
