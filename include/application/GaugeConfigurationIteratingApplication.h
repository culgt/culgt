/**
 * GaugeConfigurationIteratingApplication.h
 *
 *  Created on: Apr 15, 2014
 *      Author: vogt
 */

#ifndef GAUGECONFIGURATIONITERATINGAPPLICATION_H_
#define GAUGECONFIGURATIONITERATINGAPPLICATION_H_

#include "lattice/GaugeConfiguration.h"
#include "lattice/LatticeDimension.h"
#include "FileIterator.h"
#include "ProgramOptions.h"
#include <typeinfo>

using std::type_info;

namespace culgt
{


template<typename PatternType, typename LinkFileType, typename LinkFileTypeOut=LinkFileType> class GaugeConfigurationIteratingApplication
{
public:
	GaugeConfigurationIteratingApplication( const LatticeDimension<PatternType::SITETYPE::NDIM> dimension, FileIterator fileiterator, ProgramOptions* programOptions ): dimension(dimension), configuration(dimension), fileiterator(fileiterator), programOptions(programOptions)
	{
		configuration.allocateMemory();
	}

	virtual ~GaugeConfigurationIteratingApplication()
	{
		configuration.freeMemory();
	}

	virtual void iterate() = 0;
	virtual void setup() = 0;
	virtual void teardown() = 0;

//	virtual void addProgramOptions() = 0;

	void run()
	{
		setup();
		for( fileiterator.reset(); fileiterator.hasElement(); fileiterator.next() )
		{
			iterate();
		}
		teardown();
	}

	bool load()
	{
		linkFile->setFilename( fileiterator.getFilename() );
		std::cout << fileiterator.getFilename() ;
		try
		{
			linkFile->load( configuration.getHostPointer() );
			std::cout << " loaded!" << std::endl;
			return true;
		}
		catch( FileNotFoundException& e )
		{
			std::cout << " not found! Skipping..." << std::endl;
			return false;
		}
	}

	bool loadToDevice()
	{
		bool isLoaded = load();
		if( isLoaded )
		{
			configuration.copyToDevice();
			std::cout << "Copied to device!" << std::endl;
		}
		return isLoaded;
	}

	void save( string appendix )
	{
		linkFileOut->setFilename( fileiterator.getFilename( appendix ) );
		std::cout << fileiterator.getFilename( appendix )  << " saved!" << std::endl;
		linkFileOut->save( configuration.getHostPointer() );
	}

	void saveFromDevice( string appendix )
	{
		configuration.copyToHost();
		save( appendix );
	}

#if BOOST_VERSION < 105300
	template<typename ConcreteApplicationType> static ConcreteApplicationType* init( int argc, char* argv[] )
#else
	template<typename ConcreteApplicationType> static ConcreteApplicationType* init( const int argc, const char* argv[] )
#endif
	{
		// do the command line interpretation here.
		ProgramOptions* po = new ProgramOptions();
//		ConcreteApplicationType::addProgramOptions( po );
		po->parseOptions( argc, argv, true );

		if( po->getDevicenumber() >= 0 )
		{
			cudaSetDevice( po->getDevicenumber() );
			CUDA_LAST_ERROR( "Set device" );
		}

		std::cout << "Using device: " << DeviceProperties::getName() << "(" << DeviceProperties::getDeviceNumber() << ")" << std::endl;

		int fileNumberEnd = po->getFileNumberStart()+(po->getNConf()-1)*po->getFileNumberStep();
		FileIterator fileiterator( po->getFileBasename(), po->getFileEnding(), po->getFileNumberformat(), po->getFileNumberStart(), fileNumberEnd, po->getFileNumberStep() );


		LatticeDimension<PatternType::SITETYPE::NDIM> dimtest(po->getNt(),po->getNx(),po->getNy(),po->getNz());

		// TODO in principle we could read the sizes from the gaugeconfig file!
		APP = new ConcreteApplicationType( LatticeDimension<PatternType::SITETYPE::NDIM>(po->getNt(),po->getNx(),po->getNy(),po->getNz()), fileiterator, po );

		APP->linkFile = new LinkFileType( APP->dimension, po->getReinterpretReal() );
		if( typeid(LinkFileType) == typeid(LinkFileTypeOut) )
			APP->linkFileOut = APP->linkFile;
		else
			APP->linkFileOut = new LinkFileTypeOut( APP->dimension, po->getReinterpretReal() );

		// parse again for options that where added in the constructor
		po->parseOptions( argc, argv );

		return dynamic_cast<ConcreteApplicationType*>(APP);
	}

	static void destroy()
	{
		delete APP;
	}

	LinkFileType& getLinkFile()
	{
		return *linkFile;
	}

	const LatticeDimension<PatternType::SITETYPE::NDIM>& getDimension() const
	{
		return dimension;
	}

#if BOOST_VERSION < 105300
	template<typename ConcreteApplicationType> static void main( int argc, char* argv[] )
#else
	template<typename ConcreteApplicationType> static void main( const int argc, const char* argv[] )
#endif
	{
		init<ConcreteApplicationType>( argc, argv );
		APP->run();
		APP->destroy();
	}

protected:
	LatticeDimension<PatternType::SITETYPE::NDIM> dimension;
	LinkFileType* linkFile;
	LinkFileTypeOut* linkFileOut;
	GaugeConfiguration<PatternType> configuration;
	FileIterator fileiterator;
	ProgramOptions* programOptions;

private:
	static GaugeConfigurationIteratingApplication<PatternType, LinkFileType, LinkFileTypeOut>* APP;
};

template<typename PatternType, typename LinkFileType, typename LinkFileTypeOut> GaugeConfigurationIteratingApplication<PatternType,LinkFileType,LinkFileTypeOut>* GaugeConfigurationIteratingApplication<PatternType,LinkFileType,LinkFileTypeOut>::APP = NULL;

}

#endif /* GAUGECONFIGURATIONITERATINGAPPLICATION_H_ */
