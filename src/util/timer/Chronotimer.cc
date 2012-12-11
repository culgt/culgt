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

#include "Chronotimer.h"

Chronotimer::Chronotimer()
{
	resetted = true;
	running = false;
}

void Chronotimer::start()
{
	running = true;
	if (resetted)
		gettimeofday(&begin, NULL);
	else
	{
		__suseconds_t tmp_us = end.tv_usec - begin.tv_usec;
		__time_t tmp_s = end.tv_sec - begin.tv_sec;
		if (tmp_us < 0)
		{
			tmp_us += 1000000;
			tmp_s -= 1;
		}

		gettimeofday(&begin, NULL);

		if (begin.tv_usec - tmp_us < 0)
		{
			begin.tv_sec += -tmp_s - 1;
			begin.tv_usec += -tmp_us + 1000000;
		}
		else
		{
			begin.tv_sec += -tmp_s;
			begin.tv_usec += -tmp_us;
		}
	}
}

void Chronotimer::stop()
{
	gettimeofday(&end, NULL);
	running = false;
}

void Chronotimer::reset()
{
	resetted = true;
	running = false;
}

double Chronotimer::getTime()
{
	if (running)
	{
		gettimeofday(&end, NULL);
	}

	return (double) (end.tv_sec - begin.tv_sec)
			+ (double) (end.tv_usec - begin.tv_usec) / (double) 1000000;
}
