/**
 * Plaquette.h
 *
 *  Created on: Mar 13, 2014
 *      Author: vogt
 */

#ifndef PLAQUETTE_H_
#define PLAQUETTE_H_

template<typename PatternType> class Plaquette
{
public:
	typedef typename PatternType::PARAMTYPE::TYPE T;

private:
	T* U;
};


#endif /* PLAQUETTE_H_ */
