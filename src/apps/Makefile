coulomb:
	nvcc --ptxas-options=-v -arch=sm_20 ../util/timer/Chronotimer.cc CoulombGaugeFixingSU3_4D.cu -use_fast_math -o coulomb -D_T_=32 -D_X_=32

landau:
	nvcc --ptxas-options=-v -arch=sm_20 ../util/timer/Chronotimer.cc LandauGaugeFixingSU3_4D.cu -use_fast_math -o landau -D_T_=32 -D_X_=32


clean:
	rm coulomb
	rm landau
