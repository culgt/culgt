/*
 * MILCConverter.h
 *
 *  Created on: Aug 18, 2014
 *      Author: vogt
 */

#ifndef MILCCONVERTER_H_
#define MILCCONVERTER_H_
#include "../../lattice/LatticeDimension.h"
#include "../../lattice/LocalLink.h"
#include "../../lattice/parameterization_types/SUNRealFull.h"
#include "../../cuLGT1legacy/SiteCoord.hxx"

#include "su3.h" // milc
#include "lattice.h"


#ifdef __cplusplus
   extern "C" {
#endif

int node_index( int x, int y, int z, int t ); // external function defined in MILC code

#ifdef __cplusplus
   }
#endif


namespace culgt
{


template<typename Pattern, typename MILCRealType> class MILCConverter
{
public:
	typedef LocalLink<SUNRealFull<3,typename Pattern::PARAMTYPE::REALTYPE> > LocalLinkType;
	typedef GlobalLink<Pattern> GlobalLinkType;
	static void convertToMILC( typename Pattern::PARAMTYPE::TYPE* src, int nx, int ny, int nz, int nt )
	{
		LatticeDimension<4> dim( nt, nx, ny, nz );
		typename Pattern::SITETYPE s( dim, DO_NOT_USE_NEIGHBOURS );

		SiteCoord<Pattern::SITETYPE::NDIM, Pattern::SITETYPE::PARITYTYPE > sitecoord( dim );

		for(int t=0;t<nt;t++)
		{
			sitecoord[0] = t;
			for(int z=0;z<nz;z++)
			{
				sitecoord[3] = z;
				for(int y=0;y<ny;y++)
				{
					sitecoord[2] = y;
					for( int x=0;x<nx;x++)
					{
						sitecoord[1] = x;

						int milcindex = node_index( x, y, z, t );
						s.setIndex( sitecoord.getIndex() );
						for( int mu = 0; mu < 4; mu++ )
						{
							LocalLinkType temp;
							GlobalLinkType link( src, s, mu );

							temp = link;

							int milcmu;
							if( mu == 0 )  milcmu = 3;
							else
							{
								milcmu = mu -1;
							}

							for( int i = 0; i < 3; i++ )
							for( int j = 0; j < 3; j++ )
							{
								lattice[milcindex].link[milcmu].e[i][j].real = temp.get( (i*3+j)*2 );
								lattice[milcindex].link[milcmu].e[i][j].imag = temp.get( (i*3+j)*2 +1 );
							}

						}
					}
				}
			}
		}
	}

	static void convertFromMILC( typename Pattern::PARAMTYPE::TYPE* dest, int nx, int ny, int nz, int nt )
	{
		LatticeDimension<4> dim( nt, nx, ny, nz );
		typename Pattern::SITETYPE s( dim, DO_NOT_USE_NEIGHBOURS );

		SiteCoord<Pattern::SITETYPE::NDIM, Pattern::SITETYPE::PARITYTYPE > sitecoord( dim );

		for(int t=0;t<nt;t++)
		{
			sitecoord[0] = t;
			for(int z=0;z<nz;z++)
			{
				sitecoord[3] = z;
				for(int y=0;y<ny;y++)
				{
					sitecoord[2] = y;
					for(int x=0;x<nx;x++)
					{
						sitecoord[1] = x;

						int milcindex = node_index( x, y, z, t );
						s.setIndex( sitecoord.getIndex() );
						for( int mu = 0; mu < 4; mu++ )
						{
							LocalLinkType temp;

							int milcmu;
							if( mu == 0 )  milcmu = 3;
							else
							{
								milcmu = mu -1;
							}

							for( int i = 0; i < 3; i++ )
							for( int j = 0; j < 3; j++ )
							{
								temp.set( (i*3+j)*2, lattice[milcindex].link[milcmu].e[i][j].real );
								temp.set( (i*3+j)*2+1, lattice[milcindex].link[milcmu].e[i][j].imag );
							}

							GlobalLinkType link( dest, s, mu );

							link = temp;
						}
					}
				}
			}
		}
	}
};

}

#endif /* MILCCONVERTER_H_ */
