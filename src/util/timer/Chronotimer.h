/************************************************************************
 *
 *  Copyright 2012 Mario Schroeck, Hannes Vogt
 *
 *  This file is part of cuLGT.
 *
 *  cuLGT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  cuLGT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with cuLGT.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************
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
