/**
 */

#include "Chronotimer.h"


Chronotimer::Chronotimer() :running(false), resetted(true), savedTime(0)
{
}

void Chronotimer::start()
{
	if( !running )
	{
		gettimeofday(&begin, NULL);
		running = true;
	}
}

void Chronotimer::stop()
{
	gettimeofday(&end, NULL);
	savedTime += getElapsedSeconds();
	running = false;
}

void Chronotimer::reset()
{
	resetted = true;
	running = false;
	savedTime = 0;
}

double Chronotimer::getTime()
{
	if (running)
	{
		gettimeofday(&end, NULL);
		return savedTime + getElapsedSeconds();
	}
	else
	{
		return savedTime;
	}
}

double Chronotimer::getElapsedSeconds()
{
	return (double) (end.tv_sec - begin.tv_sec)
			+ (double) (end.tv_usec - begin.tv_usec) / (double) 1000000;
}
