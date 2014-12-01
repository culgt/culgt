/**
 *
 *  Created on: Mar 24, 2014
 *      Author: vogt
 */

#ifndef LANDAUCOULOMBLOGARITHMICGAUGETYPE_H_
#define LANDAUCOULOMBLOGARITHMICGAUGETYPE_H_

//enum GaugeTypeName{LANDAU,COULOMB};

template<GaugeTypeName gaugetype> class LandauCoulombLogarithmicGaugeType
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
//				T rnorm = q.y*q.y+q.z*q.z+q.w*q.w;
//				if( rnorm <= 0 )
//				{
//					rnorm = 1.;
//				}
//				else
//				{
//					rnorm = rsqrt( rnorm )*acos(q.x);
//				}
//				atomicAdd( &shA[id+NSB], 2*q.y*rnorm );
//				atomicAdd( &shA[id+2*NSB], 2*q.z*rnorm );
//				atomicAdd( &shA[id+3*NSB], 2*q.w*rnorm );

				//linear def
				atomicAdd( &shA[id+NSB], 2*q.y );
				atomicAdd( &shA[id+2*NSB], 2*q.z );
				atomicAdd( &shA[id+3*NSB], 2*q.w );
			}
			else
			{
//				T rnorm = q.y*q.y+q.z*q.z+q.w*q.w;
//				if( rnorm <= 0 )
//				{
//					rnorm = 1.;
//				}
//				else
//				{
//					rnorm = rsqrt( rnorm )*acos(q.x);
//				}
//				atomicAdd( &shA[id+NSB], -2*q.y*rnorm );
//				atomicAdd( &shA[id+2*NSB], -2*q.z*rnorm );
//				atomicAdd( &shA[id+3*NSB], -2*q.w*rnorm );

				// linear def
				atomicAdd( &shA[id+NSB], -2*q.y );
				atomicAdd( &shA[id+2*NSB], -2*q.z );
				atomicAdd( &shA[id+3*NSB], -2*q.w );
			}
		}
	}
	template<typename T> __device__ static inline void gatherInfo( T* shA, typename Real4<T>::VECTORTYPE& qUp, typename Real4<T>::VECTORTYPE& qDown, const short id, const short mu, const int NSB )
	{
		if( gaugetype==LANDAU || mu > 0 )
		{
//			T acosUp = (qUp.x>1)?(acos(1.)):acos(qUp.x);
//			T acosDown = (qDown.x>1)?(acos(1.)):acos(qDown.x);
//
//			T normSqUp = qUp.y*qUp.y+qUp.z*qUp.z+qUp.w*qUp.w;
//			T rnormUp;
//			if( normSqUp <= 0 )
//			{
//				rnormUp = 1.;
//			}
//			else
//			{
//				rnormUp = rsqrt( normSqUp )*acosUp;
//			}
//
//			T normSqDown = qDown.y*qDown.y+qDown.z*qDown.z+qDown.w*qDown.w;
//			T rnormDown;
//			if( normSqDown <= 0)
//			{
//				rnormDown = 1.;
//			}
//			else
//			{
//				rnormDown = rsqrt( normSqDown )*acosDown;
//			}
//			atomicAdd( &shA[id+NSB],  2*qUp.y*rnormUp-2*qDown.y*rnormDown );
//			atomicAdd( &shA[id+2*NSB], 2*qUp.z*rnormUp-2*qDown.z*rnormDown );
//			atomicAdd( &shA[id+3*NSB], 2*qUp.w*rnormUp-2*qDown.w*rnormDown );

			// linear def
			atomicAdd( &shA[id+NSB],  2*qUp.y-2*qDown.y );
			atomicAdd( &shA[id+2*NSB], 2*qUp.z-2*qDown.z );
			atomicAdd( &shA[id+3*NSB], 2*qUp.w-2*qDown.w );
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
//		if( mu == 0)
//		{
//			if( gaugetype== LANDAU)
//				add( shA, qUp, qDown, id, mu, NSB );
//		}
//		__syncthreads();
//		if( mu == 1)
//		{
//			add( shA, qUp, qDown, id, mu, NSB );
//		}
//		__syncthreads();
//		if( mu == 2)
//		{
//			add( shA, qUp, qDown, id, mu, NSB );
//		}
//		__syncthreads();
//		if( mu == 3)
//		{
//			add( shA, qUp, qDown, id, mu, NSB );
//		}
	}


	template<typename T> __device__ static inline void prepareUpdate( T* shA, const short id, const int NSB )
	{
//		if( id == 0 )
//		{
//			printf( "%f\t%f\t%f\n", shA[id+NSB] ,shA[id+2*NSB],shA[id+3*NSB]);
//		}
	}

	template<typename T> __device__ static inline void collectUpdate( T* shA, typename Real4<T>::VECTORTYPE& q, const short id, const short mu, const bool updown, const int NSB )
	{
//		if( id == 0 )
//		{
//			printf( "%f\t%f\t%f\t%f\n", shA[id],shA[id+NSB] ,shA[id+2*NSB],shA[id+3*NSB]);
//		}
		q.x = shA[id];
		q.y = -shA[id+NSB];
		q.z = -shA[id+2*NSB];
		q.w = -shA[id+3*NSB];
	}

private:
	__device__ static inline void add( double* shA, typename Real4<double>::VECTORTYPE& q, const short id, const short mu, const bool updown, const int NSB )
	{
		if( updown == 0 )
		{
//			double rnorm = q.y*q.y+q.z*q.z+q.w*q.w;
//			if( rnorm <= 0 )
//			{
//				rnorm = 1.;
//			}
//			else
//			{
//				rnorm = rsqrt( rnorm )*acos(q.x);
//			}
//			shA[id+NSB] += 2*q.y*rnorm;
//			shA[id+2*NSB] += 2*q.z*rnorm;
//			shA[id+3*NSB] += 2*q.w*rnorm;


			shA[id+NSB] += 2*q.y;
			shA[id+2*NSB] += 2*q.z;
			shA[id+3*NSB] += 2*q.w;
		}
		else
		{
//			double rnorm = q.y*q.y+q.z*q.z+q.w*q.w;
//			if( rnorm <= 0 )
//			{
//				rnorm = 1.;
//			}
//			else
//			{
//				rnorm = rsqrt( rnorm )*acos(q.x);
//			}
//			shA[id+NSB] -= 2*q.y*rnorm;
//			shA[id+2*NSB] -= 2*q.z*rnorm;
//			shA[id+3*NSB] -= 2*q.w*rnorm;

			shA[id+NSB] -= 2*q.y;
			shA[id+2*NSB] -= 2*q.z;
			shA[id+3*NSB] -= 2*q.w;
		}
	}

	__device__ static inline void add( double* shA, typename Real4<double>::VECTORTYPE& qUp, typename Real4<double>::VECTORTYPE& qDown, const short id, const short mu, const int NSB )
	{
//		double acosUp = (qUp.x>1)?(acos(1.)):acos(qUp.x);
//		double acosDown = (qDown.x>1)?(acos(1.)):acos(qDown.x);
//
//		double normSqUp = qUp.y*qUp.y+qUp.z*qUp.z+qUp.w*qUp.w;
//		double rnormUp;
//		if( normSqUp <= 0 )
//		{
//			rnormUp = 1.;
//		}
//		else
//		{
//			rnormUp = rsqrt( normSqUp )*acosUp;
//		}
//
//		double normSqDown = qDown.y*qDown.y+qDown.z*qDown.z+qDown.w*qDown.w;
//		double rnormDown;
//		if( normSqDown <= 0)
//		{
//			rnormDown = 1.;
//		}
//		else
//		{
//			rnormDown = rsqrt( normSqDown )*acosDown;
//		}
//
//
//		shA[id+NSB] +=2*qUp.y*rnormUp-2*qDown.y*rnormDown;
//		shA[id+2*NSB] +=2*qUp.z*rnormUp-2*qDown.z*rnormDown;
//		shA[id+3*NSB] +=2*qUp.w*rnormUp-2*qDown.w*rnormDown;



		shA[id+NSB] +=2*qUp.y-2*qDown.y;
		shA[id+2*NSB] +=2*qUp.z-2*qDown.z;
		shA[id+3*NSB] +=2*qUp.w-2*qDown.w;
	}

};


#endif /* LANDAUCOULOMBGAUGETYPE_H_ */
