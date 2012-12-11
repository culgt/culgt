/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
 * 
 * This class and the kernels in the namespace
 * reduces (sums over) an array Ar of doubles to a single double
 * by calling three levels of tree-like reducers.
 *
 */

#ifndef MULTIGPU_MPI_REDUCE_HXX_
#define MULTIGPU_MPI_REDUCE_HXX_


// kernels: hide in namespace
namespace MPIRED
{
__global__ void reduceStep1( double *Ar );
__global__ void reduceStep2( double *Ar );
__global__ void reduceStep3( double *Ar, int redBlockSize );
__global__ void reduceSerial( double *Ar, int size );
__global__ void setArrayZero( double *Ar );
// __device__ double average_or_max( double a, double b );
}

class Reduce
{
public:
	Reduce( int _size );
	~Reduce();
	double getReducedValue( cudaStream_t stream, double *Ar );
	void setArrayZero( cudaStream_t stream, double *Ar );
private:
	int size;
	int redBlockSize;
	double currentAr;
	void initReductionBlockSize();
	// kernel wrapper
	void reduceStep1( int a, int b, cudaStream_t stream, double *Ar );
	void reduceStep2( int a, int b, cudaStream_t stream, double *Ar );
	void reduceStep3( int a, int b, cudaStream_t stream, double *Ar, int redBlockSize );
	void reduceSerial( int a, int b, cudaStream_t stream, double *Ar );
};


#endif /* MULTIGPU_MPI_REDUCE_HXX_ */
