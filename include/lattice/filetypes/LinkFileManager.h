/*
 * LinkFileManager.h
 *
 *  Created on: Jun 22, 2015
 *      Author: vogt
 */

#ifndef LINKFILEMANAGER_H_
#define LINKFILEMANAGER_H_

#include "filetypes.h"
#include "LinkFile.h"
#include "LinkFileILDG.h"
#include "LinkFileVogt.h"
#include "LinkFileHirep.h"
#include "LinkFileHeaderOnly.h"

namespace culgt
{

template<typename MemoryPattern> class LinkFileManager
{
private:
	static const lat_dim_t memoryNdim = MemoryPattern::SITETYPE::NDIM;
	LinkFile<MemoryPattern>* linkFile;
public:
	LinkFileManager( LinkFileType filetype, const LatticeDimension<memoryNdim> size, ReinterpretReal reinterpret )
	{
		switch( filetype )
		{
		case HEADERONLY:
			linkFile = new LinkFileHeaderOnly<MemoryPattern>( size, reinterpret );
			break;
		case HIREP:
			linkFile = new LinkFileHirep<MemoryPattern>( size, reinterpret );
			break;
		case VOGT:
			linkFile = new LinkFileVogt<MemoryPattern>( size, reinterpret );
			break;
		case ILDG:
			linkFile = new LinkFileILDG<MemoryPattern>( size, reinterpret );
			break;
		}
	}

	LinkFile<MemoryPattern>* getLinkFile()
	{
		return linkFile;
	}
};

}

#endif /* LINKFILEMANAGER_H_ */
