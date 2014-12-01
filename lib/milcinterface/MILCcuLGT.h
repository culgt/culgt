/*
 * MILCcuLGT.h
 *
 *  Created on: Aug 18, 2014
 *      Author: vogt
 */

#ifndef MILCCULGT_H_
#define MILCCULGT_H_


typedef float REAL;

#ifdef __cplusplus
   extern "C" {
#endif

void cuLGTinitLandau( int nx, int ny, int nz, int nt );
void cuLGTfixLandau( int nx, int ny, int nz, int nt );


#ifdef __cplusplus
   }
#endif

#endif /* MILCCULGT_H_ */
