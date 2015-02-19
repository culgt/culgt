/**
 *
 */

#ifndef CHRONOTIMER_H_
#define CHRONOTIMER_H_

#include <string.h>
#include <sys/time.h>

class Chronotimer {
	public:
		Chronotimer();
		void start();
		void stop();
		void reset();
		double getTime();
	private:
		bool running;
		bool resetted;
		timeval begin;
		timeval end;
		double savedTime;

#ifdef __CUDACC__
		cudaEvent_t cudastart; // this is for the ChronotimerCuda hack
		cudaEvent_t cudastop; // this is for the ChronotimerCuda hack
		float cudatime;
#endif
		double getElapsedSeconds();
};

#endif /* CHRONOTIMER_H_ */
