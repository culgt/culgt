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
 *
 * This class contains a copy of these data members
 * of the class ProgramOptions which are needed by
 * the kernels. The corresponding object will be 
 * passed to the MultiGPU_MPI_Communicator object
 * when applying an algorithm.
 */

#ifndef MULTIGPU_MPI_ALGORITHMOPTIONS_H_
#define MULTIGPU_MPI_ALGORITHMOPTIONS_H_


// enum of the algorithms for MPI comm.
// OR = overrelaxation
// SR = stochastic relaxation
// SA = simulated annealing
// MS = micro step
// RT = random transformation
enum AlgoType { OR, SR, SA, MS, RT };


class MultiGPU_MPI_AlgorithmOptions
{
public:
	// constructor
	MultiGPU_MPI_AlgorithmOptions();
	// get methods
	float getOrParameter() const;
	float getSrParameter() const;
	float getSaMax() const;
	float getTemperature() const;
	float getTempStep() const;
	int getSaMicroupdates() const;
	float getSaMin() const;
	int getSaSteps() const;
	enum AlgoType getAlgorithm() const;
	long getSeed() const;
	// set methods
	void setOrParameter( float );
	void setSrParameter( float );
	void setSaMax( float );
	void setTemperature( float );
	void setTempStep( float );
	void setSaMicroupdates( int );
	void setSaMin( float );
	void setSaSteps( int );
	// set an algorithm to active
	void setAlgorithm( enum AlgoType );
	void setSeed( long );
	// temperature -= tempStep:
	void decreaseTemperature();
private:
	int saSteps;
	float saMin;
	float saMax;
	float temperature;
	float tempStep;
	int saMicroupdates;
	float orParameter;
	float srParameter;
	enum AlgoType algo;
	long seed;
};



#endif /* MULTIGPU_MPI_ALGORITHMOPTIONS_H_ */
