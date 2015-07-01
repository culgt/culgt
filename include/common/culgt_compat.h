/*
 * culgt_compat.h
 *
 *  Created on: Jul 1, 2015
 *      Author: vogt
 */

#ifndef CULGT_COMPAT_H_
#define CULGT_COMPAT_H_

#if __cplusplus >= 201103L
#define CULGT_OVERRIDE override
#define CULGT_FINAL final
#else
#define CULGT_OVERRIDE
#define CULGT_FINAL
#endif


#endif /* CULGT_COMPAT_H_ */
