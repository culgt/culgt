/**
 *
 *  Created on: Mar 24, 2014
 *      Author: vogt
 */

#ifndef RANDOMTRAFOTYPE_H_
#define RANDOMTRAFOTYPE_H_

class RandomTrafoType
{
public:
	static const int LinksInvolved = -1;
	static const int SharedArraySize = 4;

	template<typename T, typename LocalLinkType> __device__ static inline void gatherInfo( T* shA, LocalLinkType& q, lat_group_index_t iSub, lat_group_index_t jSub, const short id, const short mu, const bool updown, const int NSB )
	{
	}
	template<typename T, typename LocalLinkType> __device__ static inline void gatherInfo( T* shA, LocalLinkType& qUp, LocalLinkType& qDown, lat_group_index_t iSub, lat_group_index_t jSub, const short id, const short mu, const int NSB )
	{
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
};


#endif /* LANDAUCOULOMBGAUGETYPE_H_ */
