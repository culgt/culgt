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

// enum for the kind of gauge, i.e. Landau, Coulomb, Maximally Abelian, U(1)_3 x U(1)_8, ...
enum GaugeType {LANDAU, COULOMB, MAG, U1xU1};

// enum for the two different stopping criteria: stop if max(precision)<\eps or average(precision)<\eps
enum StoppingCrit {MAX, AVERAGE};

#endif /* LATTICE_TYPEDEFS_H_ */
