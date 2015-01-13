/*
 *  Created on: Nov 7, 2014
 *      Author: vogt
 */

#ifndef GAUGEFIXINGLAUNCHBOUNDS_H_
#define GAUGEFIXINGLAUNCHBOUNDS_H_

#include <iostream>
#include <string>
using std::ostream;
using std::string;

template<int TsitesPerBlock, int TminBlocksPerMultiprocessor> struct GaugefixingLaunchBounds
{
	static const int SitesPerBlock = TsitesPerBlock;
	static const int MinBlocksPerMultiprocessor = TminBlocksPerMultiprocessor;

	static GaugefixingLaunchBounds<TsitesPerBlock,TminBlocksPerMultiprocessor> printable()
	{
		return *(new GaugefixingLaunchBounds<TsitesPerBlock,TminBlocksPerMultiprocessor>);
	}
};

template<int TsitesPerBlock, int TminBlocksPerMultiprocessor> ostream& operator<<(ostream& out, const GaugefixingLaunchBounds<TsitesPerBlock,TminBlocksPerMultiprocessor>& lb)
{
	out << "GFLB(" << lb.SitesPerBlock << "," << lb.MinBlocksPerMultiprocessor << ")";
    return out;
}


#endif
