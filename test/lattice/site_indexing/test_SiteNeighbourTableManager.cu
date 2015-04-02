
#include "gmock/gmock.h"
#include "lattice/site_indexing/SiteIndex.h"
#include "lattice/site_indexing/SiteNeighbourTableManager.h"
#include "cudacommon/DeviceCommunicator.h"

using namespace testing;
using namespace culgt;


class ASiteNeighbourTableManager: public Test
{
public:
	typedef SiteIndex<4,NO_SPLIT> MySite;
	LatticeDimension<4> dim;

	ASiteNeighbourTableManager() : dim(4,4,4,4)
	{

	}
};

TEST_F( ASiteNeighbourTableManager, IsNotAvailableIfNotCreated )
{
	bool result = SiteNeighbourTableManager<MySite>::isAvailable( dim );

	ASSERT_FALSE( result );
}

TEST_F( ASiteNeighbourTableManager, IsAvailableIfCreated )
{
	SiteNeighbourTableManager<MySite>::generate( dim );

	bool result = SiteNeighbourTableManager<MySite>::isAvailable( dim );

	ASSERT_TRUE( result );
}

TEST_F( ASiteNeighbourTableManager, NeighbourIndexIsCorrect )
{
	MySite site( dim, SiteNeighbourTableManager<MySite>::getHostPointer( dim ) );

	site.setIndex( 0 );
	site.setNeighbour( 3, true );

	ASSERT_EQ( 1, site.getIndex() );
}

__global__ void kernelTestNeighbourTable( LatticeDimension<4> dim, lat_index_t* nn, lat_index_t* var )
{
	typedef SiteIndex<4,NO_SPLIT> MySite;
	MySite site( dim, nn );

	site.setIndex( 0 );
	site.setNeighbour( 3, true );

	var[0] = site.getIndex();
}

TEST( ASiteNeighbourTableManagerOnDevice, NeighbourIndexIsCorrect )
{
	typedef SiteIndex<4,NO_SPLIT> MySite;
	LatticeDimension<4> dim(8,8,8,8);

	lat_index_t* deviceVar;
	cudaMalloc( &deviceVar, sizeof(lat_index_t) );

	kernelTestNeighbourTable<<<1,1>>>( dim, SiteNeighbourTableManager<MySite>::getDevicePointer( dim ), deviceVar );
	CUDA_LAST_ERROR( "test kernel" );

	ASSERT_EQ( 1, DeviceCommunicator<lat_index_t>::getValue( deviceVar, 0 ) );
}
