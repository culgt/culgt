/*
 * GaugeFixingSubgroupStep.hxx
 *
 * GFlop:
 * - summation of local update: 4*2*8 = 64 (Landau), 3*2*8 = 48 (Coulomb)
 * - calculating update (see update classes)
 * - apply transformation: 4*2*3*26 (mu*(up+down)*Nc*...) = 504
 *
 *  Created on: May 10, 2012
 *      Author: vogt
 */

#ifndef GAUGEFIXINGSUBGROUPSTEP_HXX_
#define GAUGEFIXINGSUBGROUPSTEP_HXX_

#include "../lattice/datatype/datatypes.h"
#include "../lattice/datatype/lattice_typedefs.h"
#include "../lattice/Quaternion.hxx"
#include "./GlobalConstants.h"

#define USEATOMIC
#ifdef DOUBLEPRECISION
#undef USEATOMIC
#endif

template<class SUx, class Algorithm, GaugeType lc > class GaugeFixingSubgroupStep
{
public:
	__device__ inline GaugeFixingSubgroupStep( SUx* U, Algorithm algo, const short id, const short mu, const bool updown );
	__device__ inline void subgroup( const  int i, const  int j );
private:
	SUx* U;
	Algorithm algo;
	const short id;
	const short mu;
	const bool updown;
};

template<class SUx, class Algorithm, GaugeType lc > __device__ GaugeFixingSubgroupStep<SUx, Algorithm, lc>::GaugeFixingSubgroupStep( SUx* U, Algorithm algo, const short id, const short mu, const bool updown ): U(U), algo(algo), id(id), mu(mu), updown(updown)
{
}

template<class SUx, class Algorithm, GaugeType lc > __device__ void GaugeFixingSubgroupStep<SUx, Algorithm, lc>::subgroup( const int i, const int j )
{
	Quaternion<Real> q;
#ifdef USEATOMIC
	__shared__ Real shA[4*NSB];
#else
	volatile __shared__ Real shA[4*NSB];
#endif

	__syncthreads();
	if( mu == 0 )
	{
		if( updown == 0 )
		{
			shA[id]	= 0;
			shA[id+NSB] = 0;
			shA[id+2*NSB] = 0;
			shA[id+3*NSB] = 0;
		}
	}
	__syncthreads();


	if( lc == LANDAU || lc == COULOMB )
	{

#ifdef USEATOMIC
		q = U->getSubgroupQuaternion( i, j );
		if( updown == 0 )
		{
			if( lc == LANDAU || mu != 0 )
			{
				atomicAdd( &shA[id], q[0] );
				atomicAdd( &shA[id+NSB], -q[1] );
				atomicAdd( &shA[id+2*NSB], -q[2] );
				atomicAdd( &shA[id+3*NSB], -q[3] );
			}
		}
		else
		{
			if( lc == LANDAU || mu != 0 )
			{
				atomicAdd( &shA[id], q[0] );
				atomicAdd( &shA[id+NSB], q[1] );
				atomicAdd( &shA[id+2*NSB], q[2] );
				atomicAdd( &shA[id+3*NSB], q[3] );
			}
		}


	#else

	// TODO atomic_add and volatile did not work properly... can't say why... induced a race condition in an earlier version (earlier compiler?)

			if( lc == LANDAU )
			{
				if( updown == 0 && mu == 0 )
				{
				q = U->getSubgroupQuaternion( i, j );
					shA[id] += q[0];
					shA[id+NSB] -= q[1];
					shA[id+2*NSB] -= q[2];
					shA[id+3*NSB] -= q[3];
				}
				__syncthreads();
			}

			if( updown == 0 && mu == 1 )
			{
			q = U->getSubgroupQuaternion( i, j );
				shA[id] += q[0];
				shA[id+NSB] -= q[1];
				shA[id+2*NSB] -= q[2];
				shA[id+3*NSB] -= q[3];
			}
			__syncthreads();
			if( updown == 0 && mu == 2 )
			{
			q = U->getSubgroupQuaternion( i, j );
				shA[id] += q[0];
				shA[id+NSB] -= q[1];
				shA[id+2*NSB] -= q[2];
				shA[id+3*NSB] -= q[3];
			}
			__syncthreads();
			if( updown == 0 && mu == 3 )
			{
			q = U->getSubgroupQuaternion( i, j );
				shA[id] += q[0];
				shA[id+NSB] -= q[1];
				shA[id+2*NSB] -= q[2];
				shA[id+3*NSB] -= q[3];
			}
			__syncthreads();


			if( lc == LANDAU )
			{
				if( updown == 1 && mu == 0 )
				{
				q = U->getSubgroupQuaternion( i, j );
					shA[id] += q[0];
					shA[id+NSB] += q[1];
					shA[id+2*NSB] += q[2];
					shA[id+3*NSB] += q[3];
				}
				__syncthreads();
			}

			if( updown == 1 && mu == 1 )
			{
			q = U->getSubgroupQuaternion( i, j );
				shA[id] += q[0];
				shA[id+NSB] += q[1];
				shA[id+2*NSB] += q[2];
				shA[id+3*NSB] += q[3];
			}
			__syncthreads();
			if( updown == 1 && mu == 2 )
			{
			q = U->getSubgroupQuaternion( i, j );
				shA[id] += q[0];
				shA[id+NSB] += q[1];
				shA[id+2*NSB] += q[2];
				shA[id+3*NSB] += q[3];
			}
			__syncthreads();
			if( updown == 1 && mu == 3 )
			{
			q = U->getSubgroupQuaternion( i, j );
				shA[id] += q[0];
				shA[id+NSB] += q[1];
				shA[id+2*NSB] += q[2];
				shA[id+3*NSB] += q[3];
			}
#endif


			__syncthreads();
			if( mu == 0 )
			{
				if( updown == 0 )
				{
					algo.calculateUpdate( shA, id );
				}
			}
		__syncthreads();

		q[0] = shA[id];
		q[1] = shA[id+NSB];
		q[2] = shA[id+2*NSB];
		q[3] = shA[id+3*NSB];


//	q[0] = 1./2.;
//	q[1] = 1./2.;
//	q[2] = 1./2.;
//	q[3] = 1./2.;

//	q[0] = 1.;
//	q[1] = 0.;
//	q[2] = 0.;
//	q[3] = 0.;


//			if( id == 0 ) printf( "update: %f, %f\n", q[0], q[1] );

//

//	for( int i = 0; i < 3; i++ )
//	{
//		for( int j = 0; j < 3; j++ )
//		{
//			U->set(i,j,complex(0,0));
//		}
//	}

//	if( q.det().x > 1.1 || q.det().x < .9 )
//	{
//		printf( "q.det: %f; id=%d; mu=%d; updown=%d; site=%d; speicherstelle:%d\n", q.det().x, id, mu, updown, blockIdx.x * blockDim.x/8 + (threadIdx.x % NSB4) % NSB, &q[0]);
//	}


//	(*q)[0] = 1./2.;
//	(*q)[1] = 1./2.;
//	(*q)[2] = 1./2.;
//	(*q)[3] = 1./2.;
//	(*q)[0] = 1.;
//	(*q)[1] = .0;
//	(*q)[2] = .0;
//	(*q)[3] = .0;



		if( updown == 0 )
		{
			U->leftSubgroupMult( i, j, &q );
		}
		else
		{
			q.hermitian();
			U->rightSubgroupMult( i, j, &q );
		}
//		__syncthreads();

//	if( U->det().x > 1.1 || U->det().x < .9 )
//	{
//		printf( "U.det: %f; id=%d; mu=%d; updown=%d; site=%d; speicherstelle:%d\n", U->det().x, id, mu, updown, blockIdx.x * blockDim.x/8 + (threadIdx.x % NSB4) % NSB, &q[0]);
//	}






//	if( updown == 0 )
//	{
//		Real q2[4];
//		q2[0] = shA[id];
//		q2[1] = shA[id+NSB];
//		q2[2] = shA[id+2*NSB];
//		q2[3] = shA[id+3*NSB];
//		U->leftSubgroupMult( i, j, q2 );
//		Real temp = q2[0]*q2[0]-q2[3]*q2[3]+q2[2]*q2[2]-q2[1]*q2[1];
//
//		if( temp > 1.1 || temp < .9 )
//		{
//			printf( "q.det: %f; id=%d; mu=%d; updown=%d; site=%d;\n", temp, id, mu, updown, blockIdx.x * blockDim.x/8 + (threadIdx.x % NSB4) % NSB );
//		}
//	}
//	else
//	{
//		Real q2[4];
//		q2[0] = shA[id];
//		q2[1] = -shA[id+NSB];
//		q2[2] = -shA[id+2*NSB];
//		q2[3] = -shA[id+3*NSB];
//		Real temp = q2[0]*q2[0]-q2[3]*q2[3]+q2[2]*q2[2]-q2[1]*q2[1];
//
//		if( temp > 1.1 || temp < .9 )
//		{
//			printf( "q.det: %f; id=%d; mu=%d; updown=%d; site=%d;\n", temp, id, mu, updown, blockIdx.x * blockDim.x/8 + (threadIdx.x % NSB4) % NSB );
//		}
//
//		U->rightSubgroupMult( i, j, q2 );
//	}
//
//	if( U->det().x > 1.01 )
//	{
//		printf( "det > 1: %f\n", U->det().x );
//	}
	}
	else if( lc == U1xU1 )
	{
		if( updown == 0 && mu == 0 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] -= q[1];
// 			shA[id+2*NSB] -= q[2];
			shA[id+3*NSB] -= q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 1 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] -= q[1];
// 			shA[id+2*NSB] -= q[2];
			shA[id+3*NSB] -= q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 2 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] -= q[1];
// 			shA[id+2*NSB] -= q[2];
			shA[id+3*NSB] -= q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 3 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] -= q[1];
// 			shA[id+2*NSB] -= q[2];
			shA[id+3*NSB] -= q[3];
		}
		__syncthreads();


		if( updown == 1 && mu == 0 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] += q[1];
// 			shA[id+2*NSB] += q[2];
			shA[id+3*NSB] += q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 1 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] += q[1];
// 			shA[id+2*NSB] += q[2];
			shA[id+3*NSB] += q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 2 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] += q[1];
// 			shA[id+2*NSB] += q[2];
			shA[id+3*NSB] += q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 3 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
// 			shA[id+NSB] += q[1];
// 			shA[id+2*NSB] += q[2];
			shA[id+3*NSB] += q[3];
		}
		__syncthreads();

				
		if( mu == 0 )
		{
			if( updown == 0 )
			{
				shA[id+NSB]=0.0;
				shA[id+2*NSB]=0.0;
				algo.calculateUpdate( shA, id );
			}
		}
		__syncthreads();

		q[0] = shA[id];
		q[1] = 0.0;//shA[id+NSB];
		q[2] = 0.0;//shA[id+2*NSB];
		q[3] = shA[id+3*NSB];


		if( updown == 0 )
		{
			U->leftSubgroupMult( i, j, &q );
		}
		else
		{
			q.hermitian();
			U->rightSubgroupMult( i, j, &q );
		}


	}
	else if( lc == MAG )
	{
		if( updown == 0 && mu == 0 )
		{
			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
			shA[id]    +=  q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2]; //=D
			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3]; //=E
			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3]; //=F
		}
		__syncthreads();
		if( updown == 0 && mu == 1 )
		{
			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
			shA[id]    +=  q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 2 )
		{
			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
			shA[id]    +=  q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 3 )
		{
			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
			shA[id]    +=  q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();


		if( updown == 1 && mu == 0 )
		{
			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
			shA[id]    += q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 1 )
		{
			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
			shA[id]    += q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 2 )
		{
			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
			shA[id]    += q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 3 )
		{
			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
			shA[id]    += q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
		}
		__syncthreads();


		if( mu == 0 )
		{
			if( updown == 0 )
			{
				//g0 -> D+sqrt(D*D+E*E+F*F):
//				shA[id] *= .5;
//				shA[id] -= 1./(Real)4.0;
				shA[id+NSB] *= 2.;
				shA[id+2*NSB] *= 2.;
				shA[id] = shA[id]+sqrt(shA[id]*shA[id]+shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]);

				algo.calculateUpdate( shA, id );
			}
		}
		__syncthreads();

		q[0] = shA[id];
		q[1] = shA[id+NSB];
		q[2] = shA[id+2*NSB];
		q[3] = 0; // should be zero here (because it was zero before) but it is not (if we use shA[id+3*NSB]! why?)

		if( updown == 0 )
		{
			U->leftSubgroupMult( i, j, &q );
		}
		else
		{
			q.hermitian();
			U->rightSubgroupMult( i, j, &q );
		}
		//U->projectSU3withoutThirdRow();





//		if( updown == 0 && mu == 0 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3]; //=D
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3]; //=E
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3]; //=F
//		}
//		__syncthreads();
//		if( updown == 0 && mu == 1 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//		if( updown == 0 && mu == 2 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//		if( updown == 0 && mu == 3 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    +=  q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += -q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += -q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//
//
//		if( updown == 1 && mu == 0 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    += q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//		if( updown == 1 && mu == 1 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    += q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//		if( updown == 1 && mu == 2 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    += q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//		if( updown == 1 && mu == 3 )
//		{
//			q = U->getSubgroupQuaternion( i, j );
//			q.projectSU2();
//			shA[id]    += q[0]*q[0]+q[3]*q[3];
//			shA[id+NSB] += q[0]*q[1]-q[2]*q[3];
//			shA[id+2*NSB] += q[0]*q[2]+q[1]*q[3];
//		}
//		__syncthreads();
//
//
//		if( mu == 0 )
//		{
//			if( updown == 0 )
//			{
//				//g0 -> D+sqrt(D*D+E*E+F*F):
//				shA[id] /= 2.;
//				shA[id] -= 1./(Real)4.0;
//				shA[id] = shA[id]+sqrt(shA[id]*shA[id]+shA[id+NSB]*shA[id+NSB]+shA[id+2*NSB]*shA[id+2*NSB]);
////				shA[id+NSB] *= 2.;
////				shA[id+2*NSB] *= 2.;
//				algo.calculateUpdate( shA, id );
//			}
//		}
//		__syncthreads();
//
//		q[0] = shA[id];
//		q[1] = shA[id+NSB];
//		q[2] = shA[id+2*NSB];
//		q[3] = 0.0;
//
//		if( updown == 0 )
//		{
//			U->leftSubgroupMult( i, j, &q );
//		}
//		else
//		{
//			q.hermitian();
//			U->rightSubgroupMult( i, j, &q );
//		}


	}

}





#endif /* GAUGEFIXINGSUBGROUPSTEP_HXX_ */
