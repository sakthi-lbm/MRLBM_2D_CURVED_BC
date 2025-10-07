// main.cuh
#ifndef MAIN_CUH
#define MAIN_CUH

#include <stdio.h>
#include <stdlib.h>

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

// FILE INCLUDES
#include "var.h"
#include "globalStructs.h"
#include "errorDef.h"
#include "lbmInitialization.cuh"
#include "mlbm.cuh"
#include "saveData.cuh"
#include "treat_data.cuh"

/*
 *   @brief Swaps the pointers of two dfloat variables.
 *   @param pt1: reference to the first dfloat pointer to be swapped
 *   @param pt2: reference to the second dfloat pointer to be swapped
 */
__host__ __device__ void interfaceSwap(dfloat *&pt1, dfloat *&pt2)
{
	dfloat *temp = pt1;
	pt1 = pt2;
	pt2 = temp;
}

void initializeCudaEvents(cudaEvent_t &start, cudaEvent_t &stop, cudaEvent_t &start_step, cudaEvent_t &stop_step)
{
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	checkCudaErrors(cudaEventCreate(&start));
	checkCudaErrors(cudaEventCreate(&stop));
	checkCudaErrors(cudaEventCreate(&start_step));
	checkCudaErrors(cudaEventCreate(&stop_step));

	checkCudaErrors(cudaEventRecord(start, 0));
	checkCudaErrors(cudaEventRecord(start_step, 0));
}

dfloat recordElapsedTime(cudaEvent_t &start_step, cudaEvent_t &stop_step, int step)
{
	checkCudaErrors(cudaEventRecord(stop_step, 0));
	checkCudaErrors(cudaEventSynchronize(stop_step));

	float elapsedTime;
	checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start_step, stop_step));
	elapsedTime *= 0.001;

	size_t nodesUpdatedSync = step * NUMBER_LBM_NODES;
	dfloat MLUPS = (nodesUpdatedSync / 1e6) / elapsedTime;
	return MLUPS;
}

/*
 *   @brief Frees the memory allocated for the ghost interface data.
 *   @param ghostInterface: reference to the ghost interface data structure
 */
__host__ void interfaceFree(ghostInterfaceData &ghostInterface)
{
	cudaFree(ghostInterface.fGhost.X_0);
	cudaFree(ghostInterface.fGhost.X_1);
	cudaFree(ghostInterface.fGhost.Y_0);
	cudaFree(ghostInterface.fGhost.Y_1);

	cudaFree(ghostInterface.gGhost.X_0);
	cudaFree(ghostInterface.gGhost.X_1);
	cudaFree(ghostInterface.gGhost.Y_0);
	cudaFree(ghostInterface.gGhost.Y_1);
}

/*
 *   @brief Performs a CUDA memory copy for ghost interface data between source and destination.
 *   @param ghostInterface: reference to the ghost interface data structure
 *   @param dst: destination ghost data structure
 *   @param src: source ghost data structure
 *   @param kind: type of memory copy (e.g., cudaMemcpyHostToDevice)
 *   @param Q: number of quantities in the ghost data that are transfered
 */
__host__ void interfaceCudaMemcpy(GhostInterfaceData &ghostInterface, ghostData &dst, const ghostData &src, cudaMemcpyKind kind, int Q)
{
	struct MemcpyPair
	{
		dfloat *dst;
		const dfloat *src;
		size_t size;
	};

	MemcpyPair memcpyPairs[] = {
		{dst.X_0, src.X_0, sizeof(dfloat) * NUMBER_GHOST_FACE_X * Q},
		{dst.X_1, src.X_1, sizeof(dfloat) * NUMBER_GHOST_FACE_X * Q},
		{dst.Y_0, src.Y_0, sizeof(dfloat) * NUMBER_GHOST_FACE_Y * Q},
		{dst.Y_1, src.Y_1, sizeof(dfloat) * NUMBER_GHOST_FACE_Y * Q},
	};

	checkCudaErrors(cudaDeviceSynchronize());
	for (const auto &pair : memcpyPairs)
	{
		checkCudaErrors(cudaMemcpy(pair.dst, pair.src, pair.size, kind));
	}
}
/*
 *   @brief Swaps the ghost interfaces.
 *   @param ghostInterface: reference to the ghost interface data structure
 */
__host__ void swapGhostInterfaces(GhostInterfaceData &ghostInterface)
{
	// Synchronize device before performing swaps
	checkCudaErrors(cudaDeviceSynchronize());

	// Swap interface pointers for fGhost and gGhost
	interfaceSwap(ghostInterface.fGhost.X_0, ghostInterface.gGhost.X_0);
	interfaceSwap(ghostInterface.fGhost.X_1, ghostInterface.gGhost.X_1);
	interfaceSwap(ghostInterface.fGhost.Y_0, ghostInterface.gGhost.Y_0);
	interfaceSwap(ghostInterface.fGhost.Y_1, ghostInterface.gGhost.Y_1);
}

/*
 *   @brief Allocates memory for the ghost interface data.
 *   @param ghostInterface: reference to the ghost interface data structure
 */
__host__ void interfaceMalloc(ghostInterfaceData &ghostInterface)
{
	cudaMalloc((void **)&(ghostInterface.fGhost.X_0), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
	cudaMalloc((void **)&(ghostInterface.fGhost.X_1), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
	cudaMalloc((void **)&(ghostInterface.fGhost.Y_0), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);
	cudaMalloc((void **)&(ghostInterface.fGhost.Y_1), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);

	cudaMalloc((void **)&(ghostInterface.gGhost.X_0), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
	cudaMalloc((void **)&(ghostInterface.gGhost.X_1), sizeof(dfloat) * NUMBER_GHOST_FACE_X * QF);
	cudaMalloc((void **)&(ghostInterface.gGhost.Y_0), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);
	cudaMalloc((void **)&(ghostInterface.gGhost.Y_1), sizeof(dfloat) * NUMBER_GHOST_FACE_Y * QF);
}

__host__ void allocateHostMemory(
	dfloat **h_fMom, dfloat **rho, dfloat **ux, dfloat **uy, dfloat **ux_center, dfloat **uy_center)
{
	checkCudaErrors(cudaMallocHost((void **)h_fMom, MEM_SIZE_MOM));
	checkCudaErrors(cudaMallocHost((void **)rho, MEM_SIZE_SCALAR));
	checkCudaErrors(cudaMallocHost((void **)ux, MEM_SIZE_SCALAR));
	checkCudaErrors(cudaMallocHost((void **)uy, MEM_SIZE_SCALAR));
	checkCudaErrors(cudaMallocHost((void **)ux_center, L_back * sizeof(dfloat)));
	checkCudaErrors(cudaMallocHost((void **)uy_center, L_back * sizeof(dfloat)));
}

__host__ void allocateDeviceMemory(
	dfloat **d_fMom, unsigned int **dNodeType, GhostInterfaceData *ghostInterface, dfloat **ux_center, dfloat **uy_center)
{
	cudaMalloc((void **)d_fMom, MEM_SIZE_MOM);
	cudaMalloc((void **)dNodeType, sizeof(int) * NUMBER_LBM_NODES);
	interfaceMalloc(*ghostInterface);
	checkCudaErrors(cudaMalloc((void **)ux_center, L_back * sizeof(dfloat)));
	checkCudaErrors(cudaMalloc((void **)uy_center, L_back * sizeof(dfloat)));
}

__host__ void initializeDomain(
	GhostInterfaceData &ghostInterface,
	dfloat *&d_fMom, dfloat *&h_fMom,
	unsigned int *&hNodeType, unsigned int *&dNodeType, int *step,
	dim3 gridBlock, dim3 threadBlock,
#ifdef CYLINDER
	dfloat *D_max, cylinderProperties **h_cylinder_properties,
	cylinderProperties *&d_cylinder_properties,
	size_t *contour_count
#endif
)
{
	// LBM Initialization
	gpuInitialization_mom<<<gridBlock, threadBlock>>>(d_fMom);
	gpuInitialization_pop<<<gridBlock, threadBlock>>>(d_fMom, ghostInterface);

	// Node type initialization
	checkCudaErrors(cudaMallocHost((void **)&hNodeType, sizeof(unsigned int) * NUMBER_LBM_NODES));

	hostInitialization_nodeType_bulk(hNodeType);
	hostInitialization_nodeType(hNodeType);

#ifdef CYLINDER
	hostInitialization_innerNodes(hNodeType, D_max, h_cylinder_properties, contour_count);
	checkCudaErrors(cudaMalloc((void **)&d_cylinder_properties, sizeof(cylinderProperties) * (*contour_count)));
	checkCudaErrors(cudaMemcpy(d_cylinder_properties, *h_cylinder_properties, sizeof(cylinderProperties) * (*contour_count), cudaMemcpyHostToDevice));
#endif

	checkCudaErrors(cudaMemcpy(dNodeType, hNodeType, sizeof(unsigned int) * NUMBER_LBM_NODES, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaDeviceSynchronize());

	// Interface population initialization
	interfaceCudaMemcpy(ghostInterface, ghostInterface.gGhost, ghostInterface.fGhost, cudaMemcpyDeviceToDevice, QF);

	// Synchronize after all initializations
	checkCudaErrors(cudaDeviceSynchronize());

	// Synchronize and transfer data back to host if needed
	checkCudaErrors(cudaDeviceSynchronize());
	checkCudaErrors(cudaMemcpy(h_fMom, d_fMom, sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaDeviceSynchronize());
}

#endif // MAIN_CUH
