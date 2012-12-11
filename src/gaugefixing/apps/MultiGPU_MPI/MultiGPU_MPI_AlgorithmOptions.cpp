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
 *
 */

#include "./MultiGPU_MPI_AlgorithmOptions.h"

MultiGPU_MPI_AlgorithmOptions::MultiGPU_MPI_AlgorithmOptions() {
}

float MultiGPU_MPI_AlgorithmOptions::getOrParameter() const {
	return orParameter;
}

float MultiGPU_MPI_AlgorithmOptions::getSrParameter() const {
	return srParameter;
}

float MultiGPU_MPI_AlgorithmOptions::getSaMax() const {
	return saMax;
}

float MultiGPU_MPI_AlgorithmOptions::getTemperature() const {
	return temperature;
}

float MultiGPU_MPI_AlgorithmOptions::getTempStep() const {
	return tempStep;
}

int MultiGPU_MPI_AlgorithmOptions::getSaMicroupdates() const {
	return saMicroupdates;
}

float MultiGPU_MPI_AlgorithmOptions::getSaMin() const {
	return saMin;
}

int MultiGPU_MPI_AlgorithmOptions::getSaSteps() const {
	return saSteps;
}

enum AlgoType MultiGPU_MPI_AlgorithmOptions::getAlgorithm() const {
	return algo;
}

long MultiGPU_MPI_AlgorithmOptions::getSeed() const {
	return seed;
}

void MultiGPU_MPI_AlgorithmOptions::setOrParameter( float _orParameter ) {
	orParameter = _orParameter;
}

void MultiGPU_MPI_AlgorithmOptions::setSrParameter( float _srParameter ) {
	srParameter = _srParameter;
}

void MultiGPU_MPI_AlgorithmOptions::setSaMax( float _saMax ) {
	saMax = _saMax;
}

void MultiGPU_MPI_AlgorithmOptions::setTemperature( float _temperature ) {
	temperature = _temperature;
}

void MultiGPU_MPI_AlgorithmOptions::setTempStep( float _tempStep ) {
	tempStep = _tempStep;
}

void MultiGPU_MPI_AlgorithmOptions::setSaMicroupdates( int _saMicroupdates ) {
	saMicroupdates = _saMicroupdates;
}

void MultiGPU_MPI_AlgorithmOptions::setSaMin( float _saMin ) {
	saMin = _saMin;
}

void MultiGPU_MPI_AlgorithmOptions::setSaSteps( int _saSteps ) {
	saSteps = _saSteps;
}

void MultiGPU_MPI_AlgorithmOptions::setAlgorithm( enum AlgoType _algo ) {
	algo = _algo;
}

void MultiGPU_MPI_AlgorithmOptions::setSeed( long _seed ) {
	seed = _seed;
}

void MultiGPU_MPI_AlgorithmOptions::decreaseTemperature() {
	temperature -= tempStep;
}

