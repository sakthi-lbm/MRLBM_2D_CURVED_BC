#ifndef __MLBM_H
#define __MLBM_H

#include <string>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <cuda_runtime.h>
#include <builtin_types.h>
#include "var.h"
#include "includeFiles/interface.h"
#include "boundaryCondition.cuh"

#include COLREC
#include CASE_BC
#include "includeFiles/interface_handling.cuh"

#ifdef CYLINDER
#include "cases/cylinder/aux_functions.cuh"
#endif


#include "globalStructs.h"

/*
 *   @brief Updates macroscopics and then performs collision and streaming
 *   @param fMom: macroscopics moments
 *   @param ghostInterface interface block transfer information
 *   @param d_mean_rho: mean density, used for density correction
 *   @param d_BC_Fx: boundary condition force x
 *   @param d_BC_Fy: boundary condition force x
 *   @param d_BC_Fz: boundary condition force x
 *   @param step: current time step
 *   @param save: if is necessary save some data
 */
__global__ void gpuMomCollisionStream(
	dfloat* fMom, unsigned int* dNodeType, ghostInterfaceData ghostInterface, unsigned int step
);

/*
 *   @brief Updates macroscopics and then performs collision and streaming
 *   @param fMom: macroscopics moments
 *   @param ghostInterface interface block transfer information
 *   @param d_mean_rho: mean density, used for density correction
 *   @param d_BC_Fx: boundary condition force x
 *   @param d_BC_Fy: boundary condition force x
 *   @param d_BC_Fz: boundary condition force x
 *   @param step: current time step
 *   @param save: if is necessary save some data
 */
__global__ void streamingAndMom(
	dfloat* fMom, dfloat OMEGA, size_t cylinder_counter, unsigned int* dNodeType,
	ghostInterfaceData ghostInterface, cylinderProperties* cylinder_properties, unsigned int step
) ;

/*
 *   @brief Updates macroscopics and then performs collision and streaming
 *   @param fMom: macroscopics moments
 *   @param ghostInterface interface block transfer information
 *   @param d_mean_rho: mean density, used for density correction
 *   @param d_BC_Fx: boundary condition force x
 *   @param d_BC_Fy: boundary condition force x
 *   @param d_BC_Fz: boundary condition force x
 *   @param step: current time step
 *   @param save: if is necessary save some data
 */
__global__ void updateInnerBoundaries(dfloat* fMom, cylinderProperties* cylinder_properties, dfloat OMEGA, unsigned int step);


/*
 *   @brief Updates macroscopics and then performs collision and streaming
 *   @param fMom: macroscopics moments
 *   @param ghostInterface interface block transfer information
 *   @param d_mean_rho: mean density, used for density correction
 *   @param d_BC_Fx: boundary condition force x
 *   @param d_BC_Fy: boundary condition force x
 *   @param d_BC_Fz: boundary condition force x
 *   @param step: current time step
 *   @param save: if is necessary save some data
 */
__global__ void boundaryAndCollision(
    dfloat *fMom, size_t cylinder_count, dfloat OMEGA, unsigned int *dNodeType,
    ghostInterfaceData ghostInterface, cylinderProperties *cylinder_properties, unsigned int step);


#endif