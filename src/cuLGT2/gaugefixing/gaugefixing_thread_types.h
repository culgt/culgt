/**
 * gaugefixing_thread_types.h
 *
 *  Created on: Apr 8, 2014
 *      Author: vogt
 */

#ifndef GAUGEFIXING_THREAD_TYPES_H_
#define GAUGEFIXING_THREAD_TYPES_H_

namespace culgt
{
	enum GaugeFixingThreadsPerSite{ SINGLE_THREAD_PER_SITE, FOUR_THREAD_PER_SITE, EIGHT_THREAD_PER_SITE, TIMESLICE };
}

#endif /* GAUGEFIXING_THREAD_TYPES_H_ */
