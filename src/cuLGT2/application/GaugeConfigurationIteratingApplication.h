/**
 * GaugeConfigurationIteratingApplication.h
 *
 *  Created on: Apr 15, 2014
 *      Author: vogt
 */

#ifndef GAUGECONFIGURATIONITERATINGAPPLICATION_H_
#define GAUGECONFIGURATIONITERATINGAPPLICATION_H_

#include "../lattice/GaugeConfiguration.h"
#include "../lattice/LinkFile.h"
#include "../lattice/LatticeDimension.h"
#include "FileIterator.h"
#include "ProgramOptions.h"

namespace culgt
{


template<typename PatternType, typename LinkFileType> class GaugeConfigurationIteratingApplication
{
public:
	GaugeConfigurationIteratingApplication( const LatticeDimension<PatternType::SITETYPE::Ndim> dimension, FileIterator fileiterator, ProgramOptions* programOptions ): dimension(dimension), linkFile(  dimension ), configuration(dimension), fileiterator(fileiterator), programOptions(programOptions)
	{
		configuration.allocateMemory();
	}

	virtual ~GaugeConfigurationIteratingApplication()
	{
		configuration.freeMemory();
	}

	virtual void iterate() = 0;

	void run()
	{
		for( fileiterator.reset(); fileiterator.hasElement(); fileiterator.next() )
		{
			iterate();
		}
	}

	bool load()
	{
		linkFile.setFilename( fileiterator.getFilename() );
		std::cout << fileiterator.getFilename() << std::endl;
		linkFile.load( configuration.getHostPointer() );
		return false;
	}

	template<typename ConcreteApplicationType> static ConcreteApplicationType* init( const int argc, const char* argv[] )
	{
		// do the command line interpretation here.
		ProgramOptions* po = new ProgramOptions( argc, argv );

		int fileNumberEnd = po->getFileNumberStart()+(po->getNConf()-1)*po->getFileNumberStep();
		FileIterator fileiterator( po->getFileBasename(), po->getFileEnding(), po->getFileNumberformat(), po->getFileNumberStart(), fileNumberEnd );

		APP = new ConcreteApplicationType( LatticeDimension<PatternType::SITETYPE::Ndim>(32,32,32,32), fileiterator, po );

		return dynamic_cast<ConcreteApplicationType*>(APP);
	}

	static void destroy()
	{
		delete APP;
	}

	LinkFileType& getLinkFile()
	{
		return linkFile;
	}

	const LatticeDimension<PatternType::SITETYPE::Ndim>& getDimension() const
	{
		return dimension;
	}

	template<typename ConcreteApplicationType> static void main( const int argc, const char* argv[] )
	{
		init<ConcreteApplicationType>( argc, argv );
		APP->run();
		APP->destroy();
	}

protected:
	LatticeDimension<PatternType::SITETYPE::Ndim> dimension;
	LinkFileType linkFile;
	GaugeConfiguration<PatternType> configuration;
	FileIterator fileiterator;
	ProgramOptions* programOptions;

private:
	static GaugeConfigurationIteratingApplication<PatternType, LinkFileType>* APP;
};

template<typename PatternType, typename LinkFileType> GaugeConfigurationIteratingApplication<PatternType,LinkFileType>* GaugeConfigurationIteratingApplication<PatternType,LinkFileType>::APP = NULL;

}

#endif /* GAUGECONFIGURATIONITERATINGAPPLICATION_H_ */
