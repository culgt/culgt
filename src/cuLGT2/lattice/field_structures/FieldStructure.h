/**
 * FieldStructure.h
 *
 *  Created on: Jun 23, 2014
 *      Author: vogt
 */

#ifndef FIELDSTRUCTURE_H_
#define FIELDSTRUCTURE_H_
#include "../../common/culgt_typedefs.h"
#include "../../cudacommon/cuda_host_device.h"

namespace culgt
{

class FieldStructure
{
public:
	CUDA_HOST_DEVICE virtual lat_array_index_t getIndex() const = 0;
	CUDA_HOST_DEVICE virtual lat_array_index_t getSize() const = 0;
	CUDA_HOST_DEVICE virtual ~FieldStructure(){};
};

}

#endif /* FIELDSTRUCTURE_H_ */
