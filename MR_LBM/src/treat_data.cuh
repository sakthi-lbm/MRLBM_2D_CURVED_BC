#ifndef TREAT_DATA_CUH
#define TREAT_DATA_CUH

#include <stdio.h>
#include <stdlib.h>

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream> // std::cout, std::fixed
#include <iomanip>  // std::setprecision

#include "var.h"
#include "globalStructs.h"
#include "globalFunctions.h"

__host__ void calculate_pressure(cylinderProperties *h_cylinder_properties, unsigned int count, unsigned int step);

__host__ void calculate_forces(cylinderProperties *h_cylinder_properties, unsigned int count, unsigned int step);

__host__ void calculate_inlet_density(dfloat *h_fMom, unsigned int step, dfloat *rho_infty);

__global__ void domain_avg(dfloat *fMom, dfloat *ux_mean, dfloat *uy_mean, unsigned int step);

__global__ void velocity_on_centerline_average(dfloat *fMom, dfloat *ux_center, dfloat *uy_center, unsigned int step);

#endif // !TREAT_DATA_CUH