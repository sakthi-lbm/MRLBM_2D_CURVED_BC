#include "mlbm.cuh"
#include "globalStructs.h"
#include "globalFunctions.h"

__global__ void streamingAndMom(
	dfloat *fMom, dfloat OMEGA, size_t cylinder_counter, unsigned int *dNodeType,
	ghostInterfaceData ghostInterface, cylinderProperties *cylinder_properties, unsigned int step)
{
	const int x = threadIdx.x + blockDim.x * blockIdx.x;
	const int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;

	dfloat pop[Q];

	__shared__ dfloat s_pop[BLOCK_LBM_SIZE * (Q - 1)];

	// Load moments from global memory

	// rho'
	unsigned int nodeType = dNodeType[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)];
	if (nodeType == 0b11111111)
		return;
	dfloat rhoVar = RHO_0 + fMom[idxMom(threadIdx.x, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)];
	dfloat ux_t30 = fMom[idxMom(threadIdx.x, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)];
	dfloat uy_t30 = fMom[idxMom(threadIdx.x, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_xx_t45 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MXX_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_xy_t90 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MXY_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_yy_t45 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MYY_INDEX, blockIdx.x, blockIdx.y)];

	pop_reconstruction(rhoVar, ux_t30, uy_t30, m_xx_t45, m_yy_t45, m_xy_t90, pop);

	const unsigned short int xp1 = (threadIdx.x + 1 + BLOCK_NX) % BLOCK_NX;
	const unsigned short int xm1 = (threadIdx.x - 1 + BLOCK_NX) % BLOCK_NX;

	const unsigned short int yp1 = (threadIdx.y + 1 + BLOCK_NY) % BLOCK_NY;
	const unsigned short int ym1 = (threadIdx.y - 1 + BLOCK_NY) % BLOCK_NY;

	// save populations in shared memory
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 0)] = pop[1];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 1)] = pop[2];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 2)] = pop[3];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 3)] = pop[4];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 4)] = pop[5];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 5)] = pop[6];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 6)] = pop[7];
	s_pop[idxPopBlock(threadIdx.x, threadIdx.y, 7)] = pop[8];

	// sync threads of the block so all populations are saved
	__syncthreads();

	pop[1] = s_pop[idxPopBlock(xm1, threadIdx.y, 0)];
	pop[2] = s_pop[idxPopBlock(threadIdx.x, ym1, 1)];
	pop[3] = s_pop[idxPopBlock(xp1, threadIdx.y, 2)];
	pop[4] = s_pop[idxPopBlock(threadIdx.x, yp1, 3)];
	pop[5] = s_pop[idxPopBlock(xm1, ym1, 4)];
	pop[6] = s_pop[idxPopBlock(xp1, ym1, 5)];
	pop[7] = s_pop[idxPopBlock(xp1, yp1, 6)];
	pop[8] = s_pop[idxPopBlock(xm1, yp1, 7)];

	/* load pop from global in cover nodes */

	pop_load(ghostInterface, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, pop);

	dfloat invRho;

	if (nodeType != BULK)
	{
		if (nodeType > 100 && nodeType < 115)
		{
			cylinderProperties *bc_property = findCylindeProperty(cylinder_properties, cylinder_counter, x, y);

			immersedBoundaryLoop((*bc_property).is, pop, &rhoVar, &m_xx_t45, &m_yy_t45, &m_xy_t90, x, y);

			if (step >= STAT_BEGIN_TIME && step <= STAT_END_TIME && CALCULATE_FORCES)
			{
				incoming_forces(bc_property, pop);
			}
		}
		else
		{
			boundary_calculation(nodeType, &rhoVar, &ux_t30, &uy_t30, &m_xx_t45, &m_yy_t45, &m_xy_t90, pop, fMom, x, y, OMEGA);
		}
	}
	else
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

		rhoVar = pop_0 + pop_1 + pop_2 + pop_3 + pop_4 + pop_5 + pop_6 + pop_7 + pop_8;
		invRho = 1 / rhoVar;

		ux_t30 = ((pop_1 + pop_5 + pop_8) - (pop_3 + pop_6 + pop_7)) * invRho;
		uy_t30 = ((pop_2 + pop_5 + pop_6) - (pop_4 + pop_7 + pop_8)) * invRho;

		m_xx_t45 = (pop_1 + pop_3 + pop_5 + pop_6 + pop_7 + pop_8) * invRho - cs2;
		m_xy_t90 = ((pop_5 + pop_7) - (pop_6 + pop_8)) * invRho;
		m_yy_t45 = (pop_2 + pop_4 + pop_5 + pop_6 + pop_7 + pop_8) * invRho - cs2;
	}

	fMom[idxMom(threadIdx.x, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)] = rhoVar - RHO_0;

	fMom[idxMom(threadIdx.x, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)] = ux_t30;
	fMom[idxMom(threadIdx.x, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)] = uy_t30;

	fMom[idxMom(threadIdx.x, threadIdx.y, M_MXX_INDEX, blockIdx.x, blockIdx.y)] = m_xx_t45;
	fMom[idxMom(threadIdx.x, threadIdx.y, M_MXY_INDEX, blockIdx.x, blockIdx.y)] = m_xy_t90;
	fMom[idxMom(threadIdx.x, threadIdx.y, M_MYY_INDEX, blockIdx.x, blockIdx.y)] = m_yy_t45;
}

__global__ void updateInnerBoundaries(dfloat *fMom, cylinderProperties *cylinder_properties, dfloat OMEGA, unsigned int step)
{
	cylinderProperties property = cylinder_properties[threadIdx.x];

	const int xb = (int)property.xb;
	const int yb = (int)property.yb;

	const int tx = xb % BLOCK_NX;
	const int ty = yb % BLOCK_NY;

	const int bx = xb / BLOCK_NX;
	const int by = yb / BLOCK_NY;

	dfloat rhoVar = RHO_0 + fMom[idxMom(tx, ty, M_RHO_INDEX, bx, by)];

	dfloat m_xx_t45 = fMom[idxMom(tx, ty, M_MXX_INDEX, bx, by)];
	dfloat m_xy_t90 = fMom[idxMom(tx, ty, M_MXY_INDEX, bx, by)];
	dfloat m_yy_t45 = fMom[idxMom(tx, ty, M_MYY_INDEX, bx, by)];

	// for first point

	dfloat ux1;
	dfloat uy1;

	bilinear_velocity_interpolation(property.x1, property.y1, int(property.x1), int(property.y1),
									int(property.x1) + 1, int(property.y1) + 1, fMom, M_UX_INDEX, &ux1);
	bilinear_velocity_interpolation(property.x1, property.y1, int(property.x1), int(property.y1),
									int(property.x1) + 1, int(property.y1) + 1, fMom, M_UY_INDEX, &uy1);

	// for second point

	dfloat ux2;
	dfloat uy2;

	bilinear_velocity_interpolation(property.x2, property.y2, int(property.x2), int(property.y2),
									int(property.x2) + 1, int(property.y2) + 1, fMom, M_UX_INDEX, &ux2);
	bilinear_velocity_interpolation(property.x2, property.y2, int(property.x2), int(property.y2),
									int(property.x2) + 1, int(property.y2) + 1, fMom, M_UY_INDEX, &uy2);

	// moment interpolation to first point
	dfloat mxx1 = 0.0;
	dfloat myy1 = 0.0;
	dfloat mxx2 = 0.0;
	dfloat myy2 = 0.0;

	if (ROTATIONAL_COORDINATES)
	{
		bilinear_moment_interpolation(property.x1, property.y1, int(property.x1), int(property.y1),
									  int(property.x1) + 1, int(property.y1) + 1, fMom, &mxx1, &myy1);
		bilinear_moment_interpolation(property.x2, property.y2, int(property.x2), int(property.y2),
									  int(property.x2) + 1, int(property.y2) + 1, fMom, &mxx2, &myy2);
	}

	if (CALCULATE_PRESSURE && step >= STAT_BEGIN_TIME && step <= STAT_END_TIME)
	{
		dfloat rho1;
		dfloat rho2;
		dfloat rho3;

		bilinear_density_interpolation(property.x1, property.y1, int(property.x1), int(property.y1), int(property.x1) + 1, int(property.y1) + 1, fMom, M_RHO_INDEX, &rho1);
		bilinear_density_interpolation(property.x2, property.y2, int(property.x2), int(property.y2), int(property.x2) + 1, int(property.y2) + 1, fMom, M_RHO_INDEX, &rho2);
		bilinear_density_interpolation(property.x3, property.y3, int(property.x3), int(property.y3), int(property.x3) + 1, int(property.y3) + 1, fMom, M_RHO_INDEX, &rho3);

		pressure_extrapolation(property.xw, property.yw, property.x1, property.y1, property.x2, property.y2, property.x3, property.y3, rho1, rho2, rho3, &(cylinder_properties[threadIdx.x].ps));
	}

	const dfloat delta = property.dr;

	const dfloat ux_t30 = extrapolation(delta, ux1, ux2);
	const dfloat uy_t30 = extrapolation(delta, uy1, uy2);

	const dfloat m_xx_int = extrapolation(delta, mxx1, mxx2);
	const dfloat m_yy_int = extrapolation(delta, myy1, myy2);

	if (ROTATIONAL_COORDINATES)
	{
		if (RHO_STRONG)
		{
			numericalSolution_rotation(&rhoVar, ux_t30, uy_t30, &m_xx_t45, &m_xy_t90, &m_yy_t45, m_xx_int, m_yy_int, property.is, property.os, OMEGA, xb, yb);
		}
		if (RHO_EQ)
		{
			numericalSolution_rotation_rhoeq(&rhoVar, ux_t30, uy_t30, &m_xx_t45, &m_xy_t90, &m_yy_t45, m_xx_int, m_yy_int, property.is, property.os, OMEGA, xb, yb);
		}
	}
	else
	{
		numericalSolution(&rhoVar, ux_t30, uy_t30, &m_xx_t45, &m_xy_t90, &m_yy_t45, property.is, property.os, OMEGA);
	}

	fMom[idxMom(tx, ty, M_RHO_INDEX, bx, by)] = rhoVar - RHO_0;

	fMom[idxMom(tx, ty, M_UX_INDEX, bx, by)] = ux_t30;
	fMom[idxMom(tx, ty, M_UY_INDEX, bx, by)] = uy_t30;

	fMom[idxMom(tx, ty, M_MXX_INDEX, bx, by)] = m_xx_t45;
	fMom[idxMom(tx, ty, M_MXY_INDEX, bx, by)] = m_xy_t90;
	fMom[idxMom(tx, ty, M_MYY_INDEX, bx, by)] = m_yy_t45;
}

__global__ void boundaryAndCollision(
	dfloat *fMom, size_t cylinder_count, dfloat OMEGA, unsigned int *dNodeType,
	ghostInterfaceData ghostInterface, cylinderProperties *cylinder_properties, unsigned int step)
{
	const int x = threadIdx.x + blockDim.x * blockIdx.x;
	const int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;
	dfloat pop[Q];

	// Load moments from global memory

	// rho'
	unsigned int nodeType = dNodeType[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)];
	if (nodeType == 0b11111111)
		return;
	dfloat rhoVar = RHO_0 + fMom[idxMom(threadIdx.x, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)];
	dfloat ux_t30 = fMom[idxMom(threadIdx.x, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)];
	dfloat uy_t30 = fMom[idxMom(threadIdx.x, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_xx_t45 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MXX_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_xy_t90 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MXY_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_yy_t45 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MYY_INDEX, blockIdx.x, blockIdx.y)];

	ux_t30 = F_M_I_SCALE * ux_t30;
	uy_t30 = F_M_I_SCALE * uy_t30;

	m_xx_t45 = F_M_II_SCALE * (m_xx_t45);
	m_xy_t90 = F_M_IJ_SCALE * (m_xy_t90);
	m_yy_t45 = F_M_II_SCALE * (m_yy_t45);

	moment_collision(ux_t30, uy_t30, &m_xx_t45, &m_yy_t45, &m_xy_t90, OMEGA);

	pop_reconstruction(rhoVar, ux_t30, uy_t30, m_xx_t45, m_yy_t45, m_xy_t90, pop);

	if (nodeType > 100 && step >= STAT_BEGIN_TIME && step <= STAT_END_TIME && CALCULATE_FORCES)
	{
		cylinderProperties *bc_property = findCylindeProperty(cylinder_properties, cylinder_count, x, y);

		outgoing_forces(bc_property, cylinder_count, pop);
	}

	fMom[idxMom(threadIdx.x, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)] = rhoVar - RHO_0;

	fMom[idxMom(threadIdx.x, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)] = ux_t30;
	fMom[idxMom(threadIdx.x, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)] = uy_t30;

	fMom[idxMom(threadIdx.x, threadIdx.y, M_MXX_INDEX, blockIdx.x, blockIdx.y)] = m_xx_t45;
	fMom[idxMom(threadIdx.x, threadIdx.y, M_MXY_INDEX, blockIdx.x, blockIdx.y)] = m_xy_t90;
	fMom[idxMom(threadIdx.x, threadIdx.y, M_MYY_INDEX, blockIdx.x, blockIdx.y)] = m_yy_t45;

	pop_save(ghostInterface, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, x, y, pop);
}