/**
 * GaugeConfigurationHelper.h
 *
 *  Created on: Feb 26, 2014
 *      Author: vogt
 */

#ifndef GAUGECONFIGURATIONHELPER_H_
#define GAUGECONFIGURATIONHELPER_H_

template<typename T> class GaugeConfigurationHelper
{
public:
	static void allocateMemory( T** pointer, size_t size );
	static void freeMemory( T* pointer );
	static void setElement( T* pointer, int i, const T val );
	static T getElement( T* pointer, int i );
	static void copyToDevice( T* dest, T* src, size_t size );
	static void copyToHost( T* dest, T* src, size_t size );
};



#endif /* GAUGECONFIGURATIONHELPER_H_ */
