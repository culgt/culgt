/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
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

// TODO define somewhere in "/gaugefixing"
// enum for the kind of gauge, i.e. Landau, Coulomb, Maximally Abelian, U(1)_3 x U(1)_8, ...
enum GaugeType {LANDAU, COULOMB, MAG, U1xU1};

// enum for the two different stopping criteria: stop if max(precision)<\eps or average(precision)<\eps
enum StoppingCrit {MAX, AVERAGE};

// parity ordering in Site
enum ParityType { NO_SPLIT, FULL_SPLIT, TIMESLICE_SPLIT };

#endif /* LATTICE_TYPEDEFS_H_ */
