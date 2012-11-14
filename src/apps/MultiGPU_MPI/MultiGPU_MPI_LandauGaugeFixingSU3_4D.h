// wrapper functions for kernel calls etc.

#ifndef CUDA_WRAPPER_H_
#define CUDA_WRAPPER_H_

enum _cudaMemcpyKind {
	_cudaMemcpyHostToHost,
	_cudaMemcpyHostToDevice,
	_cudaMemcpyDeviceToHost,
	_cudaMemcpyDeviceToDevice,
	_cudaMemcpyDefault
};

void initDevice( const int device );

// void _cudaMalloc( void** devPtr, size_t size );
// 
// void _cudaHostAlloc( void ** pHost, size_t size, unsigned int flags );
// 
// void _cudaMemcpy( void* dst, const void* src, size_t count, enum _cudaMemcpyKind kind );
// 
// void _cudaMemcpyAsync( void * dst, const void* src, size_t count, enum _cudaMemcpyKind kind, int stream = 0 );

void _set_hot( Real* U, int counter);

void _generateGaugeQualityPerSite( Real* UtUp, Real* UtDw, lat_index_t* nnt, double *dGff, double *dA );

void _averageGaugeQuality( double* dGff, double* dA );

void _orStep( Real* UtUp, Real* UtDw, lat_index_t* nnt, bool parity, float orParameter  );

#endif /* CUDA_WRAPPER_H_ */