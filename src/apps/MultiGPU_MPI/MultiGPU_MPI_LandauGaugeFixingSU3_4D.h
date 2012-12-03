// wrapper functions for kernel calls etc.
//TODO this file is temporary.

#ifndef MULTIGPU_MPI_LANDAUGAUGEFIXINGSU3_4D_H_
#define MULTIGPU_MPI_LANDAUGAUGEFIXINGSU3_4D_H_

// initialize the device per process
void initDevice( const int device );

// set hot gauge field: kernel wrapper
void _set_hot( Real* U, int counter);

// gauge quality per site: kernel wrapper
void _generateGaugeQualityPerSite( Real* UtUp, Real* UtDw, lat_index_t* nnt, double *dGff, double *dA );

// average over the gauge qualities per site: kernel wrapper
void _averageGaugeQuality( double* dGff, double* dA );

// the overrelaxation step: kernel wrapper
void _orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter, cudaStream_t stream=0  );

#endif /* MULTIGPU_MPI_LANDAUGAUGEFIXINGSU3_4D_H_ */