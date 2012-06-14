/*
 * lattice_typedefs.h
 *
 *  Created on: Apr 26, 2012
 *      Author: vogt
 */

#ifndef LATTICE_TYPEDEFS_H_
#define LATTICE_TYPEDEFS_H_

// type for lattic dimension, i.e. 3 or 4, ...
typedef short lat_dim_t;

// type for the group indices: SU(N)
typedef short lat_group_dim_t;

// datatype for lattice coordinates, i.e. the x0, x1, ... coordinates, as wells as the corresponding lattice sizes
typedef short lat_coord_t;

// type for the lattice index, i.e. the index built from the coordinates
typedef int lat_index_t;

// type for index of the global link array
typedef int lat_array_index_t;

// type for the kind of gauge, i.e. Landau, Coulomb, ...
enum GaugeType {LANDAU, COULOMB};

#endif /* LATTICE_TYPEDEFS_H_ */
