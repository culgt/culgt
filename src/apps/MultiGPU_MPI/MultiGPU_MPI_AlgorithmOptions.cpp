/*
 *  MultiGPU_MPI_AlgorithmOptions.cpp
 *
 *  Created on: Nov 30, 2012
 *      Author: schroeck
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

