/************************************************************************
 * Wraps the reduction kernel of the CUDA samples, parts of reduction.cpp are copied...
 */

#ifndef REDUCTION_H_
#define REDUCTION_H_

#include "cudacommon/cuda_host_device.h"
#include "reduction_kernel.cu"

namespace culgt
{



template<class T> class Reduction
{
public:
	Reduction( int size );
	Reduction();
	void init();
	void getNumBlocksAndThreads(int whichKernel, int n, int maxBlocks, int maxThreads, int &blocks, int &threads);
	T reduceAll( T* d_idata );
	T reduceAllDot( T* d_idata, T* d_idata2 );
	T reduceAllDotConjugate( T* d_idata, T* d_idata2 );
	T reduceAllAbs( T* d_idata );
	T reduceAllMax( T* d_idata );
	static inline void assign( T* dest, T a );
private:
	int numBlocks;
	int numThreads;
	int size;
	int cpuFinalThreshold; // const for testing...
	int whichKernel;
	int maxThreads;
	int maxBlocks;
	T *d_odata; // summation array
	T *h_odata;
};

template<class T> culgt::Reduction<T>::Reduction(int size) : size(size)
{
	cpuFinalThreshold = 1; // const for testing...
	whichKernel = 6;
	maxThreads = 256;

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	maxBlocks = deviceProp.maxThreadsPerMultiProcessor/maxThreads*deviceProp.multiProcessorCount;
	init();
}

template<class T> void culgt::Reduction<T>::init()
{
	numBlocks = 0;
	numThreads = 0;
	getNumBlocksAndThreads(whichKernel, size, maxBlocks, maxThreads, numBlocks, numThreads);

	if (numBlocks == 1)
	{
		cpuFinalThreshold = 1;
	}

	cudaMalloc((void **) &d_odata, numBlocks*sizeof(T));
	h_odata = (T *) malloc(numBlocks*sizeof(T));
}

template<class T> void culgt::Reduction<T>::getNumBlocksAndThreads(int whichKernel, int n, int maxBlocks, int maxThreads, int &blocks, int &threads)
{
	if (whichKernel < 3)
	{
		threads = (n < maxThreads) ? reduction::nextPow2(n) : maxThreads;
		blocks = (n + threads - 1) / threads;
	}
	else
	{
		threads = (n < maxThreads*2) ? reduction::nextPow2((n + 1)/ 2) : maxThreads;
		blocks = (n + (threads * 2 - 1)) / (threads * 2);
	}


	if (whichKernel == 6)
	{
		blocks = (maxBlocks < blocks) ? maxBlocks : blocks;
	}
}

template<class T> T culgt::Reduction<T>::reduceAll( T* d_idata )
{
	int s=numBlocks;
	int kernel = whichKernel;
	T gpu_result = 0.;
	bool needReadBack = true;

	reduction::reduce<T>(size, numThreads, numBlocks, whichKernel, d_idata, d_odata);

	while (s > cpuFinalThreshold )
	{
		int threads = 0, blocks = 0;
		getNumBlocksAndThreads(kernel, s, maxBlocks, maxThreads, blocks, threads);

		reduction::reduce<T>(s, threads, blocks, kernel, d_odata, d_odata);

		if (kernel < 3)
		{
			s = (s + threads - 1) / threads;
		}
		else
		{
			s = (s + (threads*2-1)) / (threads*2);
		}
	}

	if (s > 1)
	{
		// copy result from device to host
		cudaMemcpy(h_odata, d_odata, s * sizeof(T), cudaMemcpyDeviceToHost);

		for (int i=0; i < s; i++)
		{
			gpu_result += h_odata[i];
		}

		needReadBack = false;
	}

	if (needReadBack)
	{
		// copy final sum from device to host
		cudaMemcpy(&gpu_result, d_odata, sizeof(T), cudaMemcpyDeviceToHost);
	}
	return gpu_result;
}


template<class T> T culgt::Reduction<T>::reduceAllDot( T* d_idata, T* d_idata2 )
{
	int s=numBlocks;
	int kernel = whichKernel;
	T gpu_result = 0.;
	bool needReadBack = true;

	reduction::reducedot<T,false>(size, numThreads, numBlocks, whichKernel, d_idata, d_idata2, d_odata);

	while (s > cpuFinalThreshold )
	{
		int threads = 0, blocks = 0;
		getNumBlocksAndThreads(kernel, s, maxBlocks, maxThreads, blocks, threads);

		reduction::reduce<T>(s, threads, blocks, kernel, d_odata, d_odata);

		if (kernel < 3)
		{
			s = (s + threads - 1) / threads;
		}
		else
		{
			s = (s + (threads*2-1)) / (threads*2);
		}
	}

	if (s > 1)
	{
		// copy result from device to host
		cudaMemcpy(h_odata, d_odata, s * sizeof(T), cudaMemcpyDeviceToHost);

		for (int i=0; i < s; i++)
		{
			gpu_result += h_odata[i];
		}

		needReadBack = false;
	}

	if (needReadBack)
	{
		// copy final sum from device to host
		cudaMemcpy(&gpu_result, d_odata, sizeof(T), cudaMemcpyDeviceToHost);
	}
	return gpu_result;
}

template<class T> T culgt::Reduction<T>::reduceAllDotConjugate( T* d_idata, T* d_idata2 )
{
	int s=numBlocks;
	int kernel = whichKernel;
	T gpu_result = 0.;
	bool needReadBack = true;

	reduction::reducedot<T,true>(size, numThreads, numBlocks, whichKernel, d_idata, d_idata2, d_odata);

	while (s > cpuFinalThreshold )
	{
		int threads = 0, blocks = 0;
		getNumBlocksAndThreads(kernel, s, maxBlocks, maxThreads, blocks, threads);

		reduction::reduce<T>(s, threads, blocks, kernel, d_odata, d_odata);

		if (kernel < 3)
		{
			s = (s + threads - 1) / threads;
		}
		else
		{
			s = (s + (threads*2-1)) / (threads*2);
		}
	}

	if (s > 1)
	{
		// copy result from device to host
		cudaMemcpy(h_odata, d_odata, s * sizeof(T), cudaMemcpyDeviceToHost);

		for (int i=0; i < s; i++)
		{
			gpu_result += h_odata[i];
		}

		needReadBack = false;
	}

	if (needReadBack)
	{
		// copy final sum from device to host
		cudaMemcpy(&gpu_result, d_odata, sizeof(T), cudaMemcpyDeviceToHost);
	}
	return gpu_result;
}

template<class T> T culgt::Reduction<T>::reduceAllAbs( T* d_idata )
{
	int s=numBlocks;
	int kernel = whichKernel;
	T gpu_result = 0;
	bool needReadBack = true;

	reduction::reduceabs<T>(size, numThreads, numBlocks, whichKernel, d_idata, d_odata);

	while (s > cpuFinalThreshold )
	{
		int threads = 0, blocks = 0;
		getNumBlocksAndThreads(kernel, s, maxBlocks, maxThreads, blocks, threads);

		reduction::reduce<T>(s, threads, blocks, kernel, d_odata, d_odata);

		if (kernel < 3)
		{
			s = (s + threads - 1) / threads;
		}
		else
		{
			s = (s + (threads*2-1)) / (threads*2);
		}
	}

	if (s > 1)
	{
		// copy result from device to host
		cudaMemcpy(h_odata, d_odata, s * sizeof(T), cudaMemcpyDeviceToHost);

		for (int i=0; i < s; i++)
		{
			gpu_result += h_odata[i];
		}

		needReadBack = false;
	}

	if (needReadBack)
	{
		// copy final sum from device to host
		cudaMemcpy(&gpu_result, d_odata, sizeof(T), cudaMemcpyDeviceToHost);
	}
	return gpu_result;
}

template<class T> T culgt::Reduction<T>::reduceAllMax( T* d_idata )
{
	int s=numBlocks;
	int kernel = whichKernel;
	T gpu_result = 0;
	bool needReadBack = true;

	reduction::reducemax<T>(size, numThreads, numBlocks, whichKernel, d_idata, d_odata);

	while (s > cpuFinalThreshold )
	{
		int threads = 0, blocks = 0;
		getNumBlocksAndThreads(kernel, s, maxBlocks, maxThreads, blocks, threads);

		reduction::reducemax<T>(s, threads, blocks, kernel, d_odata, d_odata);

		if (kernel < 3)
		{
			s = (s + threads - 1) / threads;
		}
		else
		{
			s = (s + (threads*2-1)) / (threads*2);
		}
	}

	if (s > 1)
	{
		// copy result from device to host
		cudaMemcpy(h_odata, d_odata, s * sizeof(T), cudaMemcpyDeviceToHost);

		for (int i=0; i < s; i++)
		{
			if( gpu_result < h_odata[i])
				gpu_result = h_odata[i];
		}

		needReadBack = false;
	}

	if (needReadBack)
	{
		// copy final sum from device to host
		cudaMemcpy(&gpu_result, d_odata, sizeof(T), cudaMemcpyDeviceToHost);
	}
	return gpu_result;
}

}

#endif

