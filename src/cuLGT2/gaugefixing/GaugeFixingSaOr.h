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

namespace culgt
{

class GaugeFixingSaOr
{
public:
	GaugeFixingSaOr( lat_index_t size )
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
	virtual RunInfo runOverrelaxation( float orParameter, int iter, int id = -1 ) = 0;
	virtual RunInfo runSimulatedAnnealing( float saMax, float saMin, int steps, int id = -1 ) = 0;
	virtual GaugeStats getGaugeStats() = 0;

	void fix( GaugeSettings settings )
	{


		double lastGff = 0;

		for( int copy = 0; copy < settings.getGaugeCopies(); copy++ )
		{
			if( settings.isRandomTrafo() )
			{
				std::cout << "fix: random Trafo" << std::endl;
				randomTrafo();
			}

			if( settings.isPrintStats() )
			{
				// TODO print stats
			}

			if( settings.getSaSteps() > 0 )
			{
				std::cout << "fix: Simulated Annealing" << std::endl;
				// TODO pass checkPrecision to saSteps if (isPrintStats)
				// TODO pass reproject to saSteps
				runSimulatedAnnealing( settings.getSaMax(), settings.getSaMin(), settings.getSaSteps() );
			}

			if( settings.getOrMaxIter() > 0 )
			{
				std::cout << "fix: Overrelaxation" << std::endl;
				// TODO pass reproject to orSteps
				runOverrelaxation( settings.getOrParameter(), settings.getOrMaxIter() );
			}

			GaugeStats stats = getGaugeStats();

			if( stats.getGff() > lastGff )
			{
				// TODO check if precision is reached
				// TODO saveConfig
			}
		}
	}

protected:

	double* dA;
	double* dGff;




//	template<typename GaugeFixingType> static void fix( GaugeFixingType& gaugefixing, GaugeSettings settings )
//	{
//		for( int copy = 0; copy < settings.getGaugeCopies(); copy++ )
//		{
//			if( settings.isRandomTrafo() )
//			{
//				gaugefixing.randomTrafo();
//			}
//
//			if( settings.isPrintStats() )
//			{
//				// TODO print stats
//			}
//
//			if( settings.getSaSteps() > 0 )
//			{
//				// TODO pass checkPrecision to saSteps if (isPrintStats)
//				// TODO pass reproject to saSteps
//				gaugefixing.sasteps( settings.getSaMax(), settings.getSaMin(), settings.getSaSteps() );
//			}
//
//			if( settings.getOrMaxIter() > 0 )
//			{
//				// TODO pass reproject to orSteps
//				gaugefixing.orSteps( settings.getOrParameter(), settings.getOrMaxIter() );
//			}
//
//			GaugeStats stats = gaugefixing.getGaugeStats();
//
//			if( stats.getGff() > lastGff )
//			{
//				// TODO check if precision is reached
//				// TODO saveConfig
//			}
//		}
//	}
};

}

#endif /* GAUGEFIXINGSAOR_H_ */
