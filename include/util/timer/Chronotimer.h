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

		double getElapsedSeconds();
};

#endif /* CHRONOTIMER_H_ */
