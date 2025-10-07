#ifndef __LBM_INITIALIZATION_CUH
#define __LBM_INITIALIZATION_CUH

#include <string>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include <builtin_types.h>
#include "errorDef.h"
#include "var.h"
#include "nodeTypeMap.h"
#include "globalFunctions.h"
#include "globalStructs.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include CASE_BC
#include COLREC

/*
 *   @brief Initializes moments with equilibrium population, with density
 *          and velocity defined in the function itself
 *   @param fMom: moments to be inialized to be initialized
 */
__global__ void gpuInitialization_mom(
	dfloat* fMom);

/*
 *   @brief Initializes populations in the intefaces based on the moments
 *          defined in the gpuInitialization_mom
 *   @param fMom: moments used to initialize the interface populations
 *   @param ghostInterface interface block transfer information
 */
__global__ void gpuInitialization_pop(
	dfloat* fMom, ghostInterfaceData ghostInterface);

/*
 *   @brief Initialize the boundary condition node type
 *   @param nodeType: node type ID
 */
__global__ void gpuInitialization_nodeType(
	unsigned int* dNodeType);

/*
 *   @brief Initialize the boundary condition node type
 *   @param nodeType: node type ID
 */
__host__ void hostInitialization_nodeType_bulk(
	unsigned int* hNodeType);

/*
 *   @brief Initialize the boundary condition node type
 *   @param nodeType: node type ID
 */
__host__ void hostInitialization_nodeType(
	unsigned int* hNodeType);


	/*
 *   @brief Initialize the boundary condition node type
 *   @param nodeType: node type ID
 */
__host__ void hostInitialization_innerNodes(
	unsigned int* hNodeType,
	dfloat* D_MAX,
	cylinderProperties** cylinder_properties,
	size_t* contour_counter
);

#endif // !__LBM_INITIALIZATION_CUH