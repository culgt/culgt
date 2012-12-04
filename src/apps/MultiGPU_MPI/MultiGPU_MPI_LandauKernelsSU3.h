/*
 * CoulombKernelsSU3.hxx
 *
 *
 *  Created on: Oct 26, 2012
 *      Author: vogt
 */

#ifndef MULTIGPU_MPI_LANDAUKERNELSSU3_H_
#define MULTIGPU_MPI_LANDAUKERNELSSU3_H_


// #include "./OrUpdate.hxx"
// #include "./MPI_ProcInfo.h"
#include "./MultiGPU_MPI_AlgorithmOptions.h"

// kernels as class members are not supported (even static): wrap the kernel calls and hide the kernels in namespace.

// kernels: 
namespace MPILKSU3
{
static const int Ndim = 4;
static const int Nc = 3;
template<class Algorithm> inline __device__ void applyOneTimeslice( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, Algorithm algorithm  );
__global__ void generateGaugeQualityPerSite( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, double *dGff, double *dA );
__global__ void restoreThirdLine( Real* U, lat_index_t* nnt );
__global__ void randomTrafo( Real* U,lat_index_t* nnt, bool parity, int counter );
__global__ void orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter );
__global__ void microStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity );
__global__ void saStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float temperature, int counter );
}

// wrappers:
class MultiGPU_MPI_LandauKernelsSU3
{
public:
//	__global__ static void heatbathStep( Real* UtDw, Real* Ut, Real* UtUp, lat_index_t* nnt, float beta, bool parity, int counter );

	// constructor
	MultiGPU_MPI_LandauKernelsSU3();

	// TODO remove static and make the init in constructor
	void initCacheConfig();
	
	void applyOneTimeslice( int a, int b, cudaStream_t stream, Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, MultiGPU_MPI_AlgorithmOptions algoOptions  );

	void generateGaugeQualityPerSite( int a, int b, cudaStream_t stream, Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, double *dGff, double *dA );
// 	static double getGaugeQualityPrefactorA();
// 	static double getGaugeQualityPrefactorGff();
// 
// 	static void restoreThirdLine( int a, int b, Real* U, lat_index_t* nnt );
// 	static void randomTrafo( int a, int b, Real* U,lat_index_t* nnt, bool parity, int counter );
// 	static void orStep( int a, int b,  Real* U, lat_index_t* nnt, bool parity, float orParameter );
// 	static void microStep( int a, int b, Real* U, lat_index_t* nnt, bool parity );
// 	static void saStep( int a, int b, Real* U, lat_index_t* nnt, bool parity, float temperature, int counter );
	
private:
	
};



#endif /* COULOMBKERNELSSU2_H_ */
