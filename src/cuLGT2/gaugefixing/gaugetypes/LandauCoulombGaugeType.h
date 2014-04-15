/**
 *
 *  Created on: Mar 24, 2014
 *      Author: vogt
 */

#ifndef LANDAUCOULOMBGAUGETYPE_H_
#define LANDAUCOULOMBGAUGETYPE_H_

enum GaugeTypeName{LANDAU,COULOMB};

template<GaugeTypeName gaugetype> class LandauCoulombGaugeType
{
public:
	static const int LinksInvolved = (gaugetype==LANDAU)?(8):(6);
	static const int SharedArraySize = 4;

	template<typename T> __device__ static inline void gatherInfo( T* shA, typename Real4<T>::VECTORTYPE& q, const short id, const short mu, const bool updown, const int NSB )
	{
		if( gaugetype==LANDAU || mu > 0 )
		{
			if( updown == 0 )
			{
				atomicAdd( &shA[id], q.x );
				atomicAdd( &shA[id+NSB], -q.y );
				atomicAdd( &shA[id+2*NSB], -q.z );
				atomicAdd( &shA[id+3*NSB], -q.w );
			}
			else
			{
				atomicAdd( &shA[id], q.x );
				atomicAdd( &shA[id+NSB], q.y );
				atomicAdd( &shA[id+2*NSB], q.z );
				atomicAdd( &shA[id+3*NSB], q.w );

			}
		}
	}
	template<typename T> __device__ static inline void gatherInfo( T* shA, typename Real4<T>::VECTORTYPE& qUp, typename Real4<T>::VECTORTYPE& qDown, const short id, const short mu, const int NSB )
	{
		if( gaugetype==LANDAU || mu > 0 )
		{
			atomicAdd( &shA[id], qDown.x+qUp.x );
			atomicAdd( &shA[id+NSB], qDown.y-qUp.y );
			atomicAdd( &shA[id+2*NSB], qDown.z-qUp.z );
			atomicAdd( &shA[id+3*NSB], qDown.w-qUp.w );
		}
	}

	__device__ static inline void gatherInfo( double* shA, typename Real4<double>::VECTORTYPE& q, const short id, const short mu, const bool updown, const int NSB )
	{
		if( updown == 0 && mu == 0)
		{
			if( gaugetype== LANDAU)
				add( shA, q, id, mu, updown, NSB );
		}
		__syncthreads();
		if( updown == 0 && mu == 1)
		{
			add( shA, q, id, mu, updown, NSB );
		}
		__syncthreads();
		if( updown == 0 && mu == 2)
		{
			add( shA, q, id, mu, updown, NSB );
		}
		__syncthreads();
		if( updown == 0 && mu == 3)
		{
			add( shA, q, id, mu, updown, NSB );
		}
		__syncthreads();
		if( updown == 1 && mu == 0)
		{
			if( gaugetype== LANDAU)
				add( shA, q, id, mu, updown, NSB );
		}
		__syncthreads();
		if( updown == 1 && mu == 1)
		{
			add( shA, q, id, mu, updown, NSB );
		}
		__syncthreads();
		if( updown == 1 && mu == 2)
		{
			add( shA, q, id, mu, updown, NSB );
		}
		__syncthreads();
		if( updown == 1 && mu == 3)
		{
			add( shA, q, id, mu, updown, NSB );
		}
	}

	__device__ static inline void gatherInfo( double* shA, typename Real4<double>::VECTORTYPE& qUp, typename Real4<double>::VECTORTYPE& qDown, const short id, const short mu, const int NSB )
	{
		if( mu == 0)
		{
			if( gaugetype== LANDAU)
				add( shA, qUp, qDown, id, mu, NSB );
		}
		__syncthreads();
		if( mu == 1)
		{
			add( shA, qUp, qDown, id, mu, NSB );
		}
		__syncthreads();
		if( mu == 2)
		{
			add( shA, qUp, qDown, id, mu, NSB );
		}
		__syncthreads();
		if( mu == 3)
		{
			add( shA, qUp, qDown, id, mu, NSB );
		}
	}


	template<typename T> __device__ static inline void prepareUpdate( T* shA, const short id, const int NSB )
	{
	}

	template<typename T> __device__ static inline void collectUpdate( T* shA, typename Real4<T>::VECTORTYPE& q, const short id, const short mu, const bool updown, const int NSB )
	{
		q.x = shA[id];
		q.y = shA[id+NSB];
		q.z = shA[id+2*NSB];
		q.w = shA[id+3*NSB];
	}

private:
	template<typename T> __device__ static inline void add( T* shA, typename Real4<T>::VECTORTYPE& q, const short id, const short mu, const bool updown, const int NSB )
	{
		if( updown == 0 )
		{
			shA[id] += q.x;
			shA[id+NSB] -= q.y;
			shA[id+2*NSB] -= q.z;
			shA[id+3*NSB] -=q.w;
		}
		else
		{
			shA[id] += q.x;
			shA[id+NSB] += q.y;
			shA[id+2*NSB] += q.z;
			shA[id+3*NSB] +=q.w;
		}
	}

	__device__ static inline void add( double* shA, typename Real4<double>::VECTORTYPE& qUp, typename Real4<double>::VECTORTYPE& qDown, const short id, const short mu, const int NSB )
	{
		shA[id] += qDown.x +qUp.x;
		shA[id+NSB] +=qDown.y +-qUp.y;
		shA[id+2*NSB] +=qDown.z +-qUp.z;
		shA[id+3*NSB] +=qDown.w +-qUp.w;
	}

};


#endif /* LANDAUCOULOMBGAUGETYPE_H_ */
