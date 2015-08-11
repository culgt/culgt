/*
 * compiler.h
 *
 *  Created on: Aug 11, 2015
 *      Author: vogt
 */

#ifndef COMPILER_H_
#define COMPILER_H_


namespace culgt
{
#define CULGT_MAKE_STRING(s) CULGT_MAKE_STRING_TMP(s)
#define CULGT_MAKE_STRING_TMP(s) #s

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 \
		+ __GNUC_MINOR__ * 100 \
        + __GNUC_PATCHLEVEL__)
#define CULGT_HOST_COMPILER "GCC " CULGT_MAKE_STRING(__GNUC__) "." CULGT_MAKE_STRING(__GNUC_MINOR__) "." CULGT_MAKE_STRING(__GNUC_PATCHLEVEL__)
#endif

#include "nvcc_version.h"
#define CULGT_NVCC_COMPILER "NVCC " CULGT_MAKE_STRING(CULGT_CUDA_VERSION)

}


#endif /* COMPILER_H_ */
