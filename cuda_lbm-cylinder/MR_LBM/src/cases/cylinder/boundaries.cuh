#ifndef BOUNDARIES_CUH
#define BOUNDARIES_CUH

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "../../var.h"
#include "../../globalFunctions.h"
#include "../../nodeTypeMap.h"
#include "numerical_solutions.cuh"
#include "aux_functions.cuh"
#include "interpolation_utilities.cuh"

__host__ __device__ inline void boundary_definition(unsigned int *nodeType, unsigned int x, unsigned int y)
{
	if (x == 0 && y == 0)
	{
#ifdef BC_Y_PERIODIC
		*nodeType = WEST;
#else
		*nodeType = SOUTH_WEST;
#endif
	}
	else if (x == 0 && y == (NY - 1))
	{
#ifdef BC_Y_PERIODIC
		*nodeType = WEST;
#else
		*nodeType = NORTH_WEST;
#endif
	}
	else if (x == (NX - 1) && y == 0)
	{
#ifdef BC_Y_PERIODIC
		*nodeType = EAST;
#else
		*nodeType = SOUTH_EAST;
#endif
	}
	else if (x == (NX - 1) && y == (NY - 1))
	{
#ifdef BC_Y_PERIODIC
		*nodeType = EAST;
#else
		*nodeType = NORTH_EAST;

#endif
	}
	else if (y == 0)
	{
#ifdef BC_Y_PERIODIC
		*nodeType = BULK;
#else
		*nodeType = SOUTH;
#endif
	}
	else if (y == (NY - 1))
	{
#ifdef BC_Y_PERIODIC
		*nodeType = BULK;
#else
		*nodeType = NORTH;
#endif
	}
	else if (x == 0)
	{
		*nodeType = WEST;
	}
	else if (x == (NX - 1))
	{
		*nodeType = EAST;
	}
	else
	{
		*nodeType = BULK;
	}
}

__device__ inline void boundary_calculation(unsigned int nodeType, dfloat *rhoVar, dfloat *ux, dfloat *uy, dfloat *mxx, dfloat *myy, dfloat *mxy, dfloat *pop, dfloat *fMom, int x, int y, dfloat OMEGA)
{
	const dfloat pop_0 = pop[0] + W0;
	
	const dfloat pop_1 = pop[1] + W1;
	const dfloat pop_2 = pop[2] + W1;
	const dfloat pop_3 = pop[3] + W1;
	const dfloat pop_4 = pop[4] + W1;

	const dfloat pop_5 = pop[5] + W2;
	const dfloat pop_6 = pop[6] + W2;
	const dfloat pop_7 = pop[7] + W2;
	const dfloat pop_8 = pop[8] + W2;

	switch (nodeType)
	{
	case NORTH:
	{
		const dfloat rhoIn = pop_0 + pop_1 + pop_2 + pop_3 + pop_5 + pop_6;
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop_1 + pop_3 + pop_5 + pop_6) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop_5 - pop_6) * inv_rhoIn;
		const dfloat myyIn = (pop_2 + pop_5 + pop_6) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		*rhoVar = 3.0 * rhoIn * (4.0 + 3.0 * (1.0 - OMEGA) * myyIn) / (9.0 + OMEGA);

		*mxx = 6.0 * rhoIn * mxxIn / (5.0 * (*rhoVar));
		*mxy = 2.0 * rhoIn * mxyIn / (*rhoVar);
		*myy = (*rhoVar + 9.0 * rhoIn * myyIn) / (6.0 * (*rhoVar));

		break;
	}
	case SOUTH:
	{
		const dfloat rhoIn = pop_0 + pop_1 + pop_3 + pop_4 + pop_7 + pop_8;
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop_1 + pop_3 + pop_7 + pop_8) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop_7 - pop_8) * inv_rhoIn;
		const dfloat myyIn = (pop_3 + pop_7 + pop_8) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		*rhoVar = 3.0 * rhoIn * (4.0 + 3.0 * (1.0 - OMEGA) * myyIn) / (9.0 + OMEGA);

		*mxx = 6.0 * rhoIn * mxxIn / (5.0 * (*rhoVar));
		*mxy = 2.0 * rhoIn * mxyIn / (*rhoVar);
		*myy = (*rhoVar + 9.0 * rhoIn * myyIn) / (6.0 * (*rhoVar));

		break;
	}
	case WEST:
	{
		const dfloat rhoIn = pop_0 + pop_2 + pop_3 + pop_4 + pop_6 + pop_7;
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop_3 + pop_6 + pop_7) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop_7 - pop_6) * inv_rhoIn;
		const dfloat myyIn = (pop_2 + pop_4 + pop_6 + pop_7) * inv_rhoIn - cs2;

		*uy = 0.0;
		*ux = U_MAX;

		const dfloat rho = (4.0 * rhoIn + 3.0 * rhoIn * mxxIn) / (3.0 - 3.0 * (*ux));	// weak conservation
		// const dfloat rho = (-6.0 * rhoIn) / (-5.0 + 3.0 * (*ux) + 3.0 * (*ux) * (*ux));	//rhoeq conservation
		*mxx = (rho + 9.0 * rhoIn * mxxIn + 3.0 * rho * (*ux)) / (6.0 * rho);
		*mxy = 2.0 * rhoIn * mxyIn / rho;
		*myy = 6.0 * rhoIn * myyIn / (5.0 * rho);

		*rhoVar = rho;

		break;
	}
	case EAST:
	{
		const dfloat rhoIn = pop_0 + pop_1 + pop_2 + pop_4 + pop_5 + pop_8;
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop_1 + pop_5 + pop_8) * inv_rhoIn - cs2;
		const dfloat mxyIn = (pop_5 - pop_8) * inv_rhoIn;
		const dfloat myyIn = (pop_2 + pop_4 + pop_5 + pop_8) * inv_rhoIn - cs2;

		const dfloat rho = RHO_0;// + fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)];
		*ux = fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)] / F_M_I_SCALE;
		*uy = fMom[idxMom(threadIdx.x - 1, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)] / F_M_I_SCALE;

		*rhoVar = rho;

		*mxx = (rho + 9.0 * rhoIn * mxxIn - 3.0 * rho * (*ux)) / (6.0 * rho);
		*mxy = (6.0 * rhoIn * mxyIn - rho * (*uy)) / (3.0 * rho);
		*myy = 6.0 * rhoIn * myyIn / (5.0 * rho);

		break;
	}

	case SOUTH_WEST:
	{
		const dfloat rhoIn = pop_0 + pop_3 + pop_4 + pop_7;
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop_3 + pop_7) * inv_rhoIn - cs2;
		const dfloat mxyIn = pop_7 * inv_rhoIn;
		const dfloat myyIn = (pop_4 + pop_7) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC)
		{
			*rhoVar = 36.0 * (rhoIn - mxyIn * rhoIn + mxyIn * OMEGA * rhoIn) / (24.0 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn - (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else
		{
			*rhoVar = -12.0 * rhoIn * (-3.0 - 3.0 * mxxIn + 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA - 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0 + 9.0 * OMEGA);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn - 6.0 * rhoIn * mxyIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn + 18.0 * rhoIn * myyIn + 7.0 * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = -2.0 * (6.0 * rhoIn * mxyIn - 9.0 * rhoIn * myyIn - (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	case SOUTH_EAST:
	{
		const dfloat rhoIn = pop_0 + pop_1 + pop_4 + pop_8;
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop_1 + pop_8) * inv_rhoIn - cs2;
		const dfloat mxyIn = -pop_8 * inv_rhoIn;
		const dfloat myyIn = (pop_4 + pop_8) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC)
		{
			*rhoVar = -36.0 * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else
		{
			*rhoVar = -12.0 * rhoIn * (-3.0 - 3.0 * mxxIn - 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA + 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0 + 9.0 * OMEGA);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn + 6.0 * rhoIn * mxyIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(-18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn - 18.0 * rhoIn * myyIn - 7.0 * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = 2.0 * (6.0 * rhoIn * mxyIn + 9.0 * rhoIn * myyIn + (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	case NORTH_WEST:
	{
		const dfloat rhoIn = pop_0 + pop_2 + pop_3 + pop_6;
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop_3 + pop_6) * inv_rhoIn - cs2;
		const dfloat mxyIn = -pop_6 * inv_rhoIn;
		const dfloat myyIn = (pop_2 + pop_6) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC)
		{
			*rhoVar = -36.0 * (mxyIn * OMEGA * rhoIn - rhoIn - mxyIn * rhoIn) / (24 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else
		{
			*rhoVar = -12.0 * rhoIn * (-3.0 - 3.0 * mxxIn - 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA + 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0 + 9.0 * OMEGA);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn + 6.0 * rhoIn * mxyIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(-18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn - 18.0 * rhoIn * myyIn - 7.0 * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = 2.0 * (6.0 * rhoIn * mxyIn + 9.0 * rhoIn * myyIn + (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	case NORTH_EAST:
	{
		const dfloat rhoIn = pop_0 + pop_1 + pop_2 + pop_5;
		const dfloat inv_rhoIn = 1.0 / rhoIn;

		const dfloat mxxIn = (pop_1 + pop_5) * inv_rhoIn - cs2;
		const dfloat mxyIn = pop_5 * inv_rhoIn;
		const dfloat myyIn = (pop_2 + pop_5) * inv_rhoIn - cs2;

		*ux = 0.0;
		*uy = 0.0;

		if (IRBC)
		{
			*rhoVar = 36.0 * (rhoIn - mxyIn * rhoIn + mxyIn * OMEGA * rhoIn) / (24.0 + OMEGA);

			*mxx = 0.0;
			*mxy = (36.0 * mxyIn * rhoIn - (*rhoVar)) / (9.0 * (*rhoVar));
			*myy = 0.0;
		}
		else
		{
			*rhoVar = -12.0 * rhoIn * (-3.0 - 3.0 * mxxIn + 7.0 * mxyIn - 3.0 * myyIn + 3.0 * mxxIn * OMEGA - 7.0 * mxyIn * OMEGA + 3 * myyIn * OMEGA) / (16.0 + 9.0 * OMEGA);

			*mxx = 2.0 * (9.0 * rhoIn * mxxIn - 6.0 * rhoIn * mxyIn + (*rhoVar)) / (9.0 * (*rhoVar));
			*mxy = -(18.0 * rhoIn * mxxIn - 132.0 * rhoIn * mxyIn + 18.0 * rhoIn * myyIn + 7.0 * (*rhoVar)) / (27.0 * (*rhoVar));
			*myy = -2.0 * (6.0 * rhoIn * mxyIn - 9.0 * rhoIn * myyIn - (*rhoVar)) / (9.0 * (*rhoVar));
		}

		break;
	}
	default:
		break;
	}
}

#endif // BOUNDARY_FUNCTIONS_CUH
