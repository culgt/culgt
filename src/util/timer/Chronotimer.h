/*
 * chronotimer.h
 *
 *  Created on: Feb 24, 2012
 *      Author: vogt
 */

#ifndef CHRONOTIMER_H_
#define CHRONOTIMER_H_

#include <sys/time.h>
#include <string.h>

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
};

#endif /* CHRONOTIMER_H_ */
