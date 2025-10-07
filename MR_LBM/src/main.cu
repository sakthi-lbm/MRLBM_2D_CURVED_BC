
#include "main.cuh"
#include <iostream>
#include <chrono>
#include "saveData.cuh"

#include <string>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include <builtin_types.h>
#include "var.h"

using namespace std;

int main()
{
	printf("BLOCK_NX: %d, BLOCK_NY: %d\n", BLOCK_NX, BLOCK_NY);

	folderSetup();

	// set cuda device
	checkCudaErrors(cudaSetDevice(GPU_INDEX));

	// variable declaration
	dfloat *d_fMom;
	ghostInterfaceData ghostInterface;
	cylinderProperties *d_cylinder_properties;
	cylinderProperties *h_cylinder_properties;

	unsigned int *dNodeType;
	unsigned int *hNodeType;

	dfloat D_Max;
	size_t countor_count;

	dfloat *h_fMom;

	dfloat *rho;

	dfloat *ux;
	dfloat *uy;

	dfloat *d_ux_center;
	dfloat *d_uy_center;

	dfloat *h_ux_center;
	dfloat *h_uy_center;

	dfloat rho_infty = 0.0;

	/* ----------------- GRID AND THREADS DEFINITION FOR LBM ---------------- */
	dim3 threadBlock(BLOCK_NX, BLOCK_NY);
	dim3 gridBlock(NUM_BLOCK_X, NUM_BLOCK_Y);

	/* ------------------------- ALLOCATION FOR CPU ------------------------- */
	int step = 0;

	allocateHostMemory(&h_fMom, &rho, &ux, &uy, &h_ux_center, &h_uy_center);

	/* -------------- ALLOCATION FOR GPU ------------- */
	allocateDeviceMemory(&d_fMom, &dNodeType, &ghostInterface, &d_ux_center, &d_uy_center);

	// Setup Streams
	cudaStream_t streamsLBM[1];
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	checkCudaErrors(cudaStreamCreate(&streamsLBM[0]));
	checkCudaErrors(cudaDeviceSynchronize());

	initializeDomain(ghostInterface, d_fMom, h_fMom, hNodeType, dNodeType,
					 &step, gridBlock, threadBlock,
#ifdef CYLINDER
					 &D_Max, &h_cylinder_properties,
					 d_cylinder_properties, &countor_count
#endif
	);

	printf("final_time: %d, begin_stat: %d\n", N_STEPS, STAT_BEGIN_TIME);
	printf("count: %zu, d_max:%f\n", countor_count, D_Max);

	const dfloat VISC = U_MAX * D_Max / RE;
	const dfloat TAU = 0.5 + 3.0 * VISC; // relaxation time
	const dfloat OMEGA = 1.0 / TAU;		 // (tau)^-1

	int avg_blockSize = 256; // Otimizado para ocupação
	int avg_gridSize = (L_back + avg_blockSize - 1) / avg_blockSize;


	/* ------------------------------ TIMER EVENTS  ------------------------------ */
	checkCudaErrors(cudaSetDevice(GPU_INDEX));
	cudaEvent_t start, stop, start_step, stop_step;
	initializeCudaEvents(start, stop, start_step, stop_step);
	/* ------------------------------ LBM LOOP ------------------------------ */
	saveSimInfo(step, 0.0, D, D_Max, countor_count, rho_infty);

	/* --------------------------------------------------------------------- */
	/* ---------------------------- BEGIN LOOP ------------------------------ */
	/* --------------------------------------------------------------------- */
	for (step = INI_STEP; step < N_STEPS; step++)
	{
#ifdef CYLINDER
		streamingAndMom<<<gridBlock, threadBlock>>>(d_fMom, OMEGA, countor_count, dNodeType, ghostInterface, d_cylinder_properties, step);
		checkCudaErrors(cudaDeviceSynchronize());
		updateInnerBoundaries<<<1, countor_count>>>(d_fMom, d_cylinder_properties, OMEGA, step);

		checkCudaErrors(cudaDeviceSynchronize());
		boundaryAndCollision<<<gridBlock, threadBlock>>>(d_fMom, countor_count, OMEGA, dNodeType, ghostInterface, d_cylinder_properties, step);
#else
		gpuMomCollisionStream<<<gridBlock, threadBlock>>>(d_fMom, dNodeType, ghostInterface, step);
#endif

		// swap interface pointers
		checkCudaErrors(cudaDeviceSynchronize());
		swapGhostInterfaces(ghostInterface);

#ifdef CYLINDER
		if (step >= STAT_BEGIN_TIME && step <= STAT_END_TIME)
		{
			checkCudaErrors(cudaDeviceSynchronize());
			checkCudaErrors(cudaMemcpy(h_cylinder_properties, d_cylinder_properties, sizeof(cylinderProperties) * countor_count, cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(h_fMom, d_fMom, sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS, cudaMemcpyDeviceToHost));

			velocity_on_centerline_average<<<avg_gridSize, avg_blockSize>>>(d_fMom, d_ux_center, d_uy_center, step);

			calculate_forces(h_cylinder_properties, countor_count, step);
			calculate_pressure(h_cylinder_properties, countor_count, step);
			calculate_inlet_density(h_fMom, step, &rho_infty);
		}
#endif // CYLINDER

		if (MACR_SAVE != 0 && step % MACR_SAVE == 0)
		{
			printf("\n----------------------------------- %d -----------------------------------\n", step);

			checkCudaErrors(cudaDeviceSynchronize());
			checkCudaErrors(cudaMemcpy(h_fMom, d_fMom, sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS, cudaMemcpyDeviceToHost));
			saveMacr(h_fMom, rho, ux, uy, step);
		}
	}

	/* --------------------------------------------------------------------- */
	/* ------------------------------ END LOOP ------------------------------ */
	/* --------------------------------------------------------------------- */

	checkCudaErrors(cudaDeviceSynchronize());

	// saving the ux and uy center velocities
	checkCudaErrors(cudaMemcpy(h_ux_center, d_ux_center, L_back * sizeof(dfloat), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_uy_center, d_uy_center, L_back * sizeof(dfloat), cudaMemcpyDeviceToHost));

	saving_centerline_data(h_ux_center, h_uy_center);

	// Calculate MLUPS

	dfloat MLUPS = recordElapsedTime(start_step, stop_step, step);
	printf("\n--------------------------- Last Time Step %06d ---------------------------\n", step);
	printf("MLUPS: %f\n", MLUPS);

	/* ------------------------------ POST ------------------------------ */
	checkCudaErrors(cudaMemcpy(h_fMom, d_fMom, sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS, cudaMemcpyDeviceToHost));
	// save info file
	saveSimInfo(step, MLUPS, D, D_Max, countor_count, rho_infty);

	/* ------------------------------ FREE ------------------------------ */
	cudaFree(d_fMom);
	cudaFree(dNodeType);
	cudaFree(hNodeType);
	cudaFree(hNodeType);
	cudaFree(h_fMom);
	cudaFree(rho);
	cudaFree(ux);
	cudaFree(uy);
	interfaceFree(ghostInterface);
	return 0;
}