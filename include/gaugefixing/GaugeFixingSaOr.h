/**
 * Abstraction for the gauge fixing logic that is same in Landau and Coulomb gauge.
 *
 *  Created on: Apr 23, 2014
 *      Author: vogt
 */

#ifndef GAUGEFIXINGSAOR_H_
#define GAUGEFIXINGSAOR_H_

#include "RunInfo.h"
#include "GaugeSettings.h"
#include "GaugeStats.h"
#include "util/timer/Chronotimer.h"

namespace culgt
{

enum GaugeFieldDefinition { GAUGEFIELD_STANDARD, GAUGEFIELD_LOGARITHMIC };

class GaugeFixingSaOr
{
public:
	GaugeFixingSaOr( lat_index_t size ): fullTimeSa(0), fullTimeOr(0), iterSa(0), iterOr(0), iterMicro(0), copyMemoryIsAllocated( false )
	{
		CUDA_SAFE_CALL( cudaMalloc( &dA, size*sizeof(double) ), "malloc dA");
		CUDA_SAFE_CALL( cudaMalloc( &dGff, size*sizeof(double) ), "malloc dGff");
	}

	virtual ~GaugeFixingSaOr()
	{
		cudaFree( dA );
		cudaFree( dGff );
	}

	virtual void randomTrafo() = 0;
	virtual void reproject() = 0;
	virtual void runOverrelaxation( float orParameter ) = 0;
	virtual RunInfo getRunInfoOverrelaxation( double time, long iter ) = 0;
	virtual void runSimulatedAnnealing( float temperature ) = 0;
	virtual RunInfo getRunInfoSimulatedAnnealing( double time, long iterSa, long iterMicro ) = 0;
	virtual void runMicrocanonical() = 0;
	virtual GaugeStats getGaugeStats( GaugeFieldDefinition definition = GAUGEFIELD_STANDARD ) = 0;

	virtual void allocateCopyMemory() = 0;
	virtual void freeCopyMemory() = 0;
	virtual void saveCopy() = 0;
	virtual void writeBackCopy() = 0;
	virtual void storeCleanCopy() = 0;
	virtual void takeCleanCopy() = 0;

	void fix( GaugeSettings settings )
	{
		this->settings = settings;

		GaugeStats stats = getGaugeStats();
		double bestGff = stats.getGff();

		bool useCopyMemory = ( settings.getGaugeCopies() > 1 );

		if( useCopyMemory && !copyMemoryIsAllocated )
		{
			allocateCopyMemory();
			copyMemoryIsAllocated = true;
		}

		if( useCopyMemory )
		{
			storeCleanCopy();
			saveCopy();
		}

		for( int copy = 0; copy < settings.getGaugeCopies(); copy++ )
		{
			if( useCopyMemory ) takeCleanCopy();

			if( settings.isRandomTrafo() )
			{
				if( settings.isPrintStats() ) std::cout << "Applying random trafo" << std::endl;
				randomTrafo();
			}

			if( settings.isPrintStats() )
			{
				checkPrecision(0);
			}

			fixSimulatedAnnealing();

			fixOverrelaxation();

			stats = getGaugeStats();

			if( stats.getGff() > bestGff )
			{
				if( settings.isPrintStats() )std::cout << "Found BETTER Copy! (" << std::setprecision(16) << stats.getGff() << ")"  << std::endl;
				bestGff = stats.getGff();

				// this might be an infinite loop
//				while( !(stats.getPrecision() < settings.getPrecision() ) )
//				{
//					// we have the best GFF but the precision is not reached
//					fixOverrelaxation();
//					stats = getGaugeStats();
//				}

				if( useCopyMemory ) saveCopy();
			}
			else
			{
				if( settings.isPrintStats() )std::cout << "DID NOT FIND BETTER COPY! (" << std::setprecision(16) << stats.getGff() << ")" << std::endl;
			}
		}
		if( useCopyMemory ) writeBackCopy();

		if( copyMemoryIsAllocated )
		{
			freeCopyMemory();
			copyMemoryIsAllocated = false;
		}
	}


protected:
	bool copyMemoryIsAllocated;

	double* dA;
	double* dGff;

	GaugeSettings settings;

	Chronotimer timerSa;
	Chronotimer timerOr;

	double fullTimeSa;
	double fullTimeOr;

	long iterSa;
	long iterOr;
	long iterMicro;



	bool check( int iteration )
	{
		if( (iteration) % settings.getReproject() == 0 )
		{
			reproject();
		}

		if( (iteration) % settings.getCheckPrecision() == 0 )
		{
			return checkPrecision( iteration );
		}

		return false;
	}

	bool checkPrecision( int iteration )
	{
		GaugeStats stats = getGaugeStats();
		if( settings.isPrintStats() ) std::cout << std::setw(14) << iteration << std::fixed << std::setprecision( 14) << std::setw(24) << stats.getGff() << std::setprecision( 10 ) << std::setw(20) << std::scientific << stats.getPrecision() << std::endl;
		return stats.getPrecision() < settings.getPrecision();
	}

	void fixSimulatedAnnealing()
	{
		if( settings.getSaSteps() > 0 )
		{
			if( settings.isPrintStats() )std::cout << "Simulated Annealing" << std::endl;

			float tStep = ( settings.getSaMax()-settings.getSaMin() )/(float)settings.getSaSteps();
			float temperature = settings.getSaMax();
			cudaDeviceSynchronize();
			timerSa.reset();
			timerSa.start();
			for( int i = 0; i < settings.getSaSteps(); i++ )
			{
				runSimulatedAnnealing( temperature );
				iterSa++;
				temperature -= tStep;

				for( int j = 0; j < settings.getMicroiter(); j++ )
				{
					runMicrocanonical();
					iterMicro++;
				}
				check( i+1 );
			}
			cudaDeviceSynchronize();
			timerSa.stop();

			fullTimeSa += timerSa.getTime();

			if( settings.isPrintStats() ) std::cout << std::fixed << std::setprecision( 2 ) << "Time for Simulated Annealing: " << timerSa.getTime() << " s" << std::endl;
			if( settings.isPrintStats() ) getRunInfoSimulatedAnnealing( fullTimeSa, iterSa, iterMicro ).print();
		}
	}

	void fixOverrelaxation()
	{
		if( settings.getOrMaxIter() > 0 )
		{
			cudaDeviceSynchronize();
			timerOr.reset();
			timerOr.start();
			if( settings.isPrintStats() )std::cout << "Overrelaxation" << std::endl;
			for( int i = 0; i < settings.getOrMaxIter(); i++ )
			{
				runOverrelaxation( settings.getOrParameter() );
				CUDA_LAST_ERROR( "fixOverrelaxation()runOverrelaxtion()" );

				iterOr++;

				if( check( i+1 ) ) break;
			}
			cudaDeviceSynchronize();
			timerOr.stop();

			fullTimeOr += timerOr.getTime();

			if( settings.isPrintStats() ) std::cout << std::fixed << std::setprecision( 2 ) << "Time for Overrelaxation: " << timerOr.getTime() << " s" << std::endl;
			if( settings.isPrintStats() ) getRunInfoOverrelaxation( fullTimeOr, iterOr ).print();
		}
	}

};

}

#endif /* GAUGEFIXINGSAOR_H_ */
