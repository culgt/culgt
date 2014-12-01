/*
 *  Created on: Nov 7, 2014
 *      Author: vogt
 */

#ifndef GAUGEFIXINGLAUNCHBOUNDS_H_
#define GAUGEFIXINGLAUNCHBOUNDS_H_

template<int TsitesPerBlock, int TminBlocksPerMultiprocessor> struct GaugefixingLaunchBounds
{
	static const int SitesPerBlock = TsitesPerBlock;
	static const int MinBlocksPerMultiprocessor = TminBlocksPerMultiprocessor;
};


#endif
