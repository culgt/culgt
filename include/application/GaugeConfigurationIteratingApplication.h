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
#include "cudacommon/DeviceHelper.h"
#include "lattice/filetypes/LinkFileManager.h"

using std::type_info;

namespace culgt
{


template<typename PatternType> class GaugeConfigurationIteratingApplication
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
		fileiterator.setFileExtension( makeFileExtension( programOptions->getFileExtension(), linkFileManager->getLinkFile()->getPreferredExtension() ) );
		linkFileManager->getLinkFile()->setFilename( fileiterator.getFilename() );
		std::cout << fileiterator.getFilename() ;
		try
		{
			linkFileManager->getLinkFile()->load( configuration.getHostPointer() );
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
		std::cout << fileiterator.getFilenameWithExtension( appendix, makeFileExtension( programOptions->getFileExtensionOut(), linkFileManagerOut->getLinkFile()->getPreferredExtension() ) )  << " saved!" << std::endl;
		linkFileManagerOut->getLinkFile()->setFilename( fileiterator.getFilenameWithExtension( appendix, makeFileExtension( programOptions->getFileExtensionOut(), linkFileManagerOut->getLinkFile()->getPreferredExtension() ) ) );
		linkFileManagerOut->getLinkFile()->save( configuration.getHostPointer() );
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
		po->parseOptions( argc, argv, true );

		if( po->getDevicenumber() >= 0 )
		{
			cudaSetDevice( po->getDevicenumber() );
			CUDA_LAST_ERROR( "Set device" );
		}
		else
		{
			int selectedDevice = DeviceHelper::selectAvailableDevice();
			if( selectedDevice == -1 )
			{
				std::cout << "No device available" << std::endl;
				exit(1);
			}
			else
			{
				std::cout << "Device " << selectedDevice << " is available." << std::endl;
			}
		}

		std::cout << "Using device: " << DeviceProperties::getName() << "(" << DeviceProperties::getDeviceNumber() << ")" << std::endl;

		int fileNumberEnd = po->getFileNumberStart()+(po->getNConf()-1)*po->getFileNumberStep();


		FileIterator fileiterator( po->getFileBasename(), po->getFileExtension(), po->getFileNumberformat(), po->getFileNumberStart(), fileNumberEnd, po->getFileNumberStep() );
		// TODO in principle we could read the sizes from the gaugeconfig file!
		APP = new ConcreteApplicationType( LatticeDimension<PatternType::SITETYPE::NDIM>(po->getNt(),po->getNx(),po->getNy(),po->getNz()), fileiterator, po );

		initLinkFiles( po );

		// parse again for options that where added in the constructor
		po->parseOptions( argc, argv );

		return dynamic_cast<ConcreteApplicationType*>(APP);
	}

	static void initLinkFiles( ProgramOptions* po )
	{
		APP->linkFileManager = new LinkFileManager<PatternType>( po->getFiletype(), APP->dimension, po->getReinterpretReal() );
		if( po->getFiletypeOut() == LinkFileType::DEFAULT || po->getFiletype() == po->getFiletypeOut() )
			APP->linkFileManagerOut = APP->linkFileManager;
		else
			APP->linkFileManagerOut = new LinkFileManager<PatternType>( po->getFiletypeOut(), APP->dimension, po->getReinterpretReal() );
	}

	static void destroy()
	{
		delete APP;
	}

	LinkFileType::FileType& getLinkFile()
	{
		return *linkFileManager->getLinkFile();
	}

	const LatticeDimension<PatternType::SITETYPE::NDIM>& getDimension() const
	{
		return dimension;
	}

	std::string makeFileExtension( string poFileExtension, string filetypeFileExtension )
	{
		if( poFileExtension.size() > 0 )
		{
			return poFileExtension;
		}
		else
		{
			return filetypeFileExtension;
		}
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
	GaugeConfiguration<PatternType> configuration;
	LinkFileManager<PatternType>* linkFileManager;
	LinkFileManager<PatternType>* linkFileManagerOut;
	FileIterator fileiterator;
	ProgramOptions* programOptions;

private:
	static GaugeConfigurationIteratingApplication<PatternType>* APP;
};

template<typename PatternType> GaugeConfigurationIteratingApplication<PatternType>* GaugeConfigurationIteratingApplication<PatternType>::APP = NULL;

}

#endif /* GAUGECONFIGURATIONITERATINGAPPLICATION_H_ */
