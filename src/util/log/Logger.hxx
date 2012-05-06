/**
 *
 * @author Hannes Vogt (hannes@havogt.de) Universitaet Tuebingen - Institut fuer Theoretische Physik
 * @date 2012-04-13
 */

#ifndef LOGGER_H_
#define LOGGER_H_

#include <string>
#include <stdio.h>
#include "../cuda/cuda_host_device.h"


namespace util
{
enum LOGLEVEL { DEBUG, WARN, ERROR, FATAL };

CUDA_DEVICE LOGLEVEL curLevel = DEBUG;

class Logger
{
public:
	CUDA_HOST_DEVICE static void log( LOGLEVEL level, const char* text );
	CUDA_HOST_DEVICE static void setLevel( LOGLEVEL level );
	Logger(){};
private:
};


/**
 * set log-level
 * DEBUG = All messages
 * WARN = warn, error, fatal messages
 * ERROR = error, fatal messages
 * FATAL = fatal messages
 */
void Logger::setLevel( LOGLEVEL level )
{
	util::curLevel = level;
}

/**
 * prints log via printf (for compatibility reasons to CUDA)
 */
void Logger::log( LOGLEVEL level, const char* text )
{

#ifdef LOGON
	if( level >= util::curLevel )
	{
		printf( "%s\n", text );
	}
#endif
}

}


#endif /* LOGGER_H_ */
