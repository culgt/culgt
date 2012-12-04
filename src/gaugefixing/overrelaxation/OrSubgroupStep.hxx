/*
 *  Created on: Apr 20, 2012
 *      Author: vogt
 */

#ifndef ORSUBGROUPSTEP_HXX_
#define ORSUBGROUPSTEP_HXX_

#error "Use GaugeFixingSubgroupStep!"


template<class SUx> class OrSubgroupStep
{
public:
	__device__ inline OrSubgroupStep( SUx* U, const short id, const short mu, const bool updown );
	__device__ inline void subgroup( const  int i, const  int j );
private:
	SUx* U;
//	volatile Real* shA;
	const short id;
	const short mu;
	const bool updown;
//	Quaternion<Real>* q;
};


template<class SUx> __device__ OrSubgroupStep<SUx>::OrSubgroupStep( SUx* U, const short id, const short mu, const bool updown ): U(U), id(id), mu(mu), updown(updown)
{
}

template<class SUx> __device__ void OrSubgroupStep<SUx>::subgroup( const int i, const int j )
{
	volatile __shared__ Real shA[128];
	Quaternion<Real> q;

	__syncthreads();
	if( mu == 0 )
	{
		if( updown == 0 )
		{
			shA[id]	= 0;
			shA[id+32] = 0;
			shA[id+64] = 0;
			shA[id+96] = 0;
		}
	}
	__syncthreads();



// TODO atomic_add and volatile did not work properly... can't say why... induced a race condition
		if( updown == 0 && mu == 1 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+32] -= q[1];
			shA[id+64] -= q[2];
			shA[id+96] -= q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 2 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+32] -= q[1];
			shA[id+64] -= q[2];
			shA[id+96] -= q[3];
		}
		__syncthreads();
		if( updown == 0 && mu == 3 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+32] -= q[1];
			shA[id+64] -= q[2];
			shA[id+96] -= q[3];
		}
		__syncthreads();



		if( updown == 1 && mu == 1 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+32] += q[1];
			shA[id+64] += q[2];
			shA[id+96] += q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 2 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+32] += q[1];
			shA[id+64] += q[2];
			shA[id+96] += q[3];
		}
		__syncthreads();
		if( updown == 1 && mu == 3 )
		{
		q = U->getSubgroupQuaternion( i, j );
			shA[id] += q[0];
			shA[id+32] += q[1];
			shA[id+64] += q[2];
			shA[id+96] += q[3];
		}
		__syncthreads();


//	if( mu != 0) // COULOMB: maybe put a switch here for Landau gauge
//	{
//		q = U->getSubgroupQuaternion( i, j );
//		if( id ==31 && blockIdx.x == 0 ) printf( "%f %f %f %f\n", q[0], q[1], q[2], q[3] );
//			if( id==0 && blockIdx.x == 0 ) printf( "rein: %f %f %f %f\n", shA[id], shA[id+32], shA[id+64], shA[id+96] );
//		if( updown == 0 )
//		{
//			Real temp = -q[0];
//			atomicAdd( &shA[id], temp );
//			temp = -q[1];
//			atomicAdd( &shA[id+32], temp );
//			temp = -q[2];
//			atomicAdd( &shA[id+64], temp);
//			temp = -q[3];
//			atomicAdd( &shA[id+96], temp );
////			shA[id] -= q[0];
////			shA[id+32] -= q[1];
////			shA[id+64] -= q[2];
////			shA[id+96] -= q[3];
//		}
//		else
//		{
//			Real temp = q[0];
//			atomicAdd( &shA[id], temp );
//			temp = q[1];
//			atomicAdd( &shA[id+32], temp );
//			temp = q[2];
//			atomicAdd( &shA[id+64], temp );
//			temp = q[3];
//			atomicAdd( &shA[id+96], temp );
////			shA[id] += q[0];
////			shA[id+32] += q[1];
////			shA[id+64] += q[2];
////			shA[id+96] += q[3];
//		}
//		if( threadIdx.x == 32 && blockIdx.x == 0 ) printf( "nach add: %f %f %f %f\n", shA[id], shA[id+32], shA[id+64], shA[id+96] );
//			if( id==0 && blockIdx.x == 0 ) printf( "raus: %f %f %f %f\n", shA[id], shA[id+32], shA[id+64], shA[id+96] );
//	}

	__syncthreads();

	if( mu == 0 )
	{
		if( updown == 0 )
		{
			// first order
			Real orParameter = 1.7;

			Real ai_sq = shA[id+32]*shA[id+32]+shA[id+64]*shA[id+64]+shA[id+96]*shA[id+96];
			Real a0_sq = shA[id]*shA[id];

			Real b=(orParameter*a0_sq+ai_sq)/(a0_sq+ai_sq);
			Real c=rsqrt(a0_sq+b*b*ai_sq);

			shA[id]*=c;
			shA[id+32]*=b*c;
			shA[id+64]*=b*c;
			shA[id+96]*=b*c;






			// second order
//			Real orParameter = 1.7;
//			Real ai_sq=shA[id+32]*shA[id+32]+shA[id+64]*shA[id+64]+shA[id+96]*shA[id+96];
//			Real a0_sq=shA[id+0]*shA[id+0];
//
//
//			Real sq = sqrt( a0_sq + ai_sq );
//			shA[id+0] /= sq;
//			shA[id+32] /= sq;
//			shA[id+64] /= sq;
//			shA[id+96] /= sq;
//
//			Real tmp = 1 + orParameter*(shA[id+0]-1) + .5*orParameter*(orParameter-1)*((shA[id+0]-1)*(shA[id+0]-1)-shA[id+96]*shA[id+96]-shA[id+64]*shA[id+64]-shA[id+32]*shA[id+32]);
//
//			shA[id+32] = orParameter*shA[id+32] + .5*orParameter*(orParameter-1)*2*(shA[id+0]-1)*shA[id+32];
//			shA[id+64] = orParameter*shA[id+64] + .5*orParameter*(orParameter-1)*2*(shA[id+0]-1)*shA[id+64];
//			shA[id+96] = orParameter*shA[id+96] + .5*orParameter*(orParameter-1)*2*(shA[id+0]-1)*shA[id+96];
//
//			shA[id+0] = tmp;
//
//			ai_sq=shA[id+32]*shA[id+32]+shA[id+64]*shA[id+64]+shA[id+96]*shA[id+96];
//			a0_sq=shA[id+0]*shA[id+0];
//
//			sq = sqrt( a0_sq + ai_sq );
//			shA[id+0] /= sq;
//			shA[id+32] /= sq;
//			shA[id+64] /= sq;
//			shA[id+96] /= sq;








//			if( id==0 && blockIdx.x == 0 ) printf( "in calc: %f %f %f %f\n", shA[id], shA[id+32], shA[id+64], shA[id+96] );
		}
	}

	__syncthreads();

	q[0] = shA[id];
	q[1] = shA[id+32];
	q[2] = shA[id+64];
	q[3] = shA[id+96];

//	if( id == 0 && blockIdx.x == 0 ) printf( "%f %f %f %f\n", q[0], q[1], q[2], q[3] );
//	assert( false);



//	q[0] = 1./2.;
//	q[1] = 1./2.;
//	q[2] = 1./2.;
//	q[3] = 1./2.;

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
//		printf( "q.det: %f; id=%d; mu=%d; updown=%d; site=%d; speicherstelle:%d\n", q.det().x, id, mu, updown, blockIdx.x * blockDim.x/8 + (threadIdx.x % 128) % 32, &q[0]);
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








//	if( updown == 0 )
//	{
//		Real q2[4];
//		q2[0] = shA[id];
//		q2[1] = shA[id+32];
//		q2[2] = shA[id+64];
//		q2[3] = shA[id+96];
//		U->leftSubgroupMult( i, j, q2 );
//		Real temp = q2[0]*q2[0]-q2[3]*q2[3]+q2[2]*q2[2]-q2[1]*q2[1];
//
//		if( temp > 1.1 || temp < .9 )
//		{
//			printf( "q.det: %f; id=%d; mu=%d; updown=%d; site=%d;\n", temp, id, mu, updown, blockIdx.x * blockDim.x/8 + (threadIdx.x % 128) % 32 );
//		}
//	}
//	else
//	{
//		Real q2[4];
//		q2[0] = shA[id];
//		q2[1] = -shA[id+32];
//		q2[2] = -shA[id+64];
//		q2[3] = -shA[id+96];
//		Real temp = q2[0]*q2[0]-q2[3]*q2[3]+q2[2]*q2[2]-q2[1]*q2[1];
//
//		if( temp > 1.1 || temp < .9 )
//		{
//			printf( "q.det: %f; id=%d; mu=%d; updown=%d; site=%d;\n", temp, id, mu, updown, blockIdx.x * blockDim.x/8 + (threadIdx.x % 128) % 32 );
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



#endif /* ORSUBGROUPSTEP_HXX_ */
