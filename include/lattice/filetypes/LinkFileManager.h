/*
 * LinkFileManager.h
 *
 *  Created on: Jun 22, 2015
 *      Author: vogt
 */

#ifndef LINKFILEMANAGER_H_
#define LINKFILEMANAGER_H_

#include "filetype_config.h"
#include "filetypes.h"
#include "LinkFile.h"
#include "LinkFileVogt.h"
#include "LinkFileHirep.h"
#include "LinkFileHeaderOnly.h"
#include "LinkFileNERSC.h"

#ifdef CULGT_HAVE_LINKFILE_ILDG
#include "LinkFileILDG.h"
#endif

using namespace culgt::LinkFileType;

namespace culgt
{

template<typename MemoryPattern> class LinkFileManager
{
private:
	static const lat_dim_t memoryNdim = MemoryPattern::SITETYPE::NDIM;
	LinkFile<MemoryPattern>* linkFile;
public:
	LinkFileManager( FileType filetype, const LatticeDimension<memoryNdim> size, ReinterpretReal reinterpret )
	{
		switch( filetype )
		{
		case DEFAULT:
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
#ifdef CULGT_HAVE_LINKFILE_ILDG
			linkFile = new LinkFileILDG<MemoryPattern>( size, reinterpret );
			break;
#else
			throw LinkFileException( "ILDG file format not supported (recompile with required libraries)");
#endif
		case NERSC:
			if( MemoryPattern::PARAMTYPE::NC == 3 )
			{
				linkFile = new LinkFileNERSC<MemoryPattern>( size, reinterpret );
			}
			else
				throw LinkFileException( "NERSC supported for SU(3) only" );
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
