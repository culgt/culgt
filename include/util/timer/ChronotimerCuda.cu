/**
 *
 * This is more ore less a hack to allow easy switching between both implementations...
 *
 */



#include "Chronotimer.h"

Chronotimer::Chronotimer()
{
	cudaEventCreate(&cudastart);
	cudaEventCreate(&cudastop);
	resetted = true;
	running = false;
}

void Chronotimer::start()
{
	if( !running )
	{
		cudaEventRecord(cudastart, 0);
		running = true;
	}
}

void Chronotimer::stop()
{
	cudaEventRecord(cudastop, 0);
	cudaEventSynchronize(cudastop);
	running = false;
}

void Chronotimer::reset()
{
	resetted = true;
	running = false;
}

double Chronotimer::getTime()
{
	cudaEventElapsedTime(&cudatime, cudastart, cudastop);
	return (double)cudatime/1000.;
}

double Chronotimer::getElapsedSeconds()
{
	return 0;
}
