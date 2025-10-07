#include "lbmInitialization.cuh"
#include <cmath>

__global__ void gpuInitialization_mom(
	dfloat *fMom)
{
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y;
	if (x >= NX || y >= NY)
		return;

	// first moments
	dfloat rho, ux, uy;

	rho = RHO_0;
	ux = U_0_X;
	uy = U_0_Y;

	// zeroth moment
	fMom[idxMom(threadIdx.x, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)] = rho - RHO_0;
	fMom[idxMom(threadIdx.x, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)] = F_M_I_SCALE * ux;
	fMom[idxMom(threadIdx.x, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)] = F_M_I_SCALE * uy;

	// second moments
	// define equilibrium populations
	dfloat pop[Q];
	for (int i = 0; i < Q; i++)
	{
		pop[i] = w[i] * RHO_0 * (1.0 + 3.0 * (ux * cx[i] + uy * cy[i]) + 4.5 * (ux * ux * (cx[i] * cx[i] - cs2) + uy * uy * (cx[i] * cx[i] - cs2)) + 9 * ux * uy * cx[i] * cy[i]);
	}

	dfloat invRho = 1.0 / rho;
	dfloat pixx = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;
	dfloat pixy = ((pop[5] + pop[7]) - (pop[6] + pop[8])) * invRho;
	dfloat piyy = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7] + pop[8]) * invRho - cs2;

	fMom[idxMom(threadIdx.x, threadIdx.y, M_MXX_INDEX, blockIdx.x, blockIdx.y)] = F_M_II_SCALE * pixx;
	fMom[idxMom(threadIdx.x, threadIdx.y, M_MXY_INDEX, blockIdx.x, blockIdx.y)] = F_M_IJ_SCALE * pixy;
	fMom[idxMom(threadIdx.x, threadIdx.y, M_MYY_INDEX, blockIdx.x, blockIdx.y)] = F_M_II_SCALE * piyy;
}

__global__ void gpuInitialization_pop(
	dfloat *fMom, ghostInterfaceData ghostInterface)
{
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y;
	if (x >= NX || y >= NY)
		return;

	// zeroth moment
	dfloat rhoVar = RHO_0 + fMom[idxMom(threadIdx.x, threadIdx.y, M_RHO_INDEX, blockIdx.x, blockIdx.y)];
	dfloat ux_t30 = fMom[idxMom(threadIdx.x, threadIdx.y, M_UX_INDEX, blockIdx.x, blockIdx.y)];
	dfloat uy_t30 = fMom[idxMom(threadIdx.x, threadIdx.y, M_UY_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_xx_t45 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MXX_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_xy_t90 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MXY_INDEX, blockIdx.x, blockIdx.y)];
	dfloat m_yy_t45 = fMom[idxMom(threadIdx.x, threadIdx.y, M_MYY_INDEX, blockIdx.x, blockIdx.y)];

	dfloat pop[Q];

	pop_reconstruction(rhoVar, ux_t30, uy_t30, m_xx_t45, m_yy_t45, m_xy_t90, pop);

	// thread xyz
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// block xyz
	int bx = blockIdx.x;
	int by = blockIdx.y;

	if (threadIdx.x == 0)
	{ // w
		ghostInterface.fGhost.X_0[idxPopX(ty, 0, bx, by)] = pop[3];
		ghostInterface.fGhost.X_0[idxPopX(ty, 1, bx, by)] = pop[6];
		ghostInterface.fGhost.X_0[idxPopX(ty, 2, bx, by)] = pop[7];
	}
	else if (threadIdx.x == (BLOCK_NX - 1))
	{
		ghostInterface.fGhost.X_1[idxPopX(ty, 0, bx, by)] = pop[1];
		ghostInterface.fGhost.X_1[idxPopX(ty, 1, bx, by)] = pop[5];
		ghostInterface.fGhost.X_1[idxPopX(ty, 2, bx, by)] = pop[8];
	}

	if (threadIdx.y == 0)
	{ // s
		ghostInterface.fGhost.Y_0[idxPopY(tx, 0, bx, by)] = pop[4];
		ghostInterface.fGhost.Y_0[idxPopY(tx, 1, bx, by)] = pop[7];
		ghostInterface.fGhost.Y_0[idxPopY(tx, 2, bx, by)] = pop[8];
	}
	else if (threadIdx.y == (BLOCK_NY - 1))
	{
		ghostInterface.fGhost.Y_1[idxPopY(tx, 0, bx, by)] = pop[2];
		ghostInterface.fGhost.Y_1[idxPopY(tx, 1, bx, by)] = pop[5];
		ghostInterface.fGhost.Y_1[idxPopY(tx, 2, bx, by)] = pop[6];
	}
}

__global__ void gpuInitialization_nodeType(
	unsigned int *dNodeType)
{
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y;

	if (x >= NX || y >= NY)
		return;

	unsigned int nodeType;

	boundary_definition(&nodeType, x, y);

	dNodeType[idxScalarBlock(threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y)] = nodeType;
}

__host__ void hostInitialization_nodeType_bulk(
	unsigned int *hNodeType)
{
	int x, y;

	for (y = 0; y < NY; y++)
	{
		for (x = 0; x < NX; x++)
		{
			hNodeType[idxScalarBlock(x % BLOCK_NX, y % BLOCK_NY, x / BLOCK_NX, y / BLOCK_NY)] = BULK;
		}
	}
}

__host__ void hostInitialization_nodeType(
	unsigned int *hNodeType)
{
	int x, y;
	unsigned int nodeType;

	for (y = 0; y < NY; y++)
	{
		for (x = 0; x < NX; x++)
		{

			boundary_definition(&nodeType, x, y);

			if (nodeType != BULK)
				hNodeType[idxScalarBlock(x % BLOCK_NX, y % BLOCK_NY, x / BLOCK_NX, y / BLOCK_NY)] = (unsigned int)nodeType;
		}
	}
}

__host__ void hostInitialization_innerNodes(
	unsigned int *hNodeType,
	dfloat *D_MAX,
	cylinderProperties **cylinder_properties,
	size_t *contour_counter)
{
	int x, y;
	unsigned nodeType;

	float r = (float)(D * 0.5);

	for (x = 1; x < NX - 1; x++)
	{
		for (y = 1; y < NY - 1; y++)
		{
			float x_diff = x - xc;
			float y_diff = y - yc;

			float dist = sqrt(x_diff * x_diff + y_diff * y_diff);
			float value = dist - r;

			if ((value >= -0.5 && value <= 0.5) || dist <= r)
			{
				hNodeType[idxScalarBlock(x % BLOCK_NX, y % BLOCK_NY, x / BLOCK_NX, y / BLOCK_NY)] = SOLID_NODE;
			}
		}
	}

	int count = 0;
	int max_count = 1; // estimativa inicial

	float max_radius = 0.0;

	cudaMallocHost((void **)cylinder_properties, sizeof(cylinderProperties) * max_count);

	for (int y = L_bot - 2; y < L_bot + D + 2; y++)
	{
		for (int x = L_front - 2; x < L_front + D + 2; x++)
		{

			const unsigned short int xp1 = x + 1;
			const unsigned short int xm1 = x - 1;

			const unsigned short int yp1 = y + 1;
			const unsigned short int ym1 = y - 1;

			int node_0 = hNodeType[idxScalarBlock(x % BLOCK_NX, y % BLOCK_NY, x / BLOCK_NX, y / BLOCK_NY)];
			int node_1 = hNodeType[idxScalarBlock(xp1 % BLOCK_NX, y % BLOCK_NY, xp1 / BLOCK_NX, y / BLOCK_NY)];
			int node_2 = hNodeType[idxScalarBlock(x % BLOCK_NX, yp1 % BLOCK_NY, x / BLOCK_NX, yp1 / BLOCK_NY)];
			int node_3 = hNodeType[idxScalarBlock(xm1 % BLOCK_NX, y % BLOCK_NY, xm1 / BLOCK_NX, y / BLOCK_NY)];
			int node_4 = hNodeType[idxScalarBlock(x % BLOCK_NX, ym1 % BLOCK_NY, x / BLOCK_NX, ym1 / BLOCK_NY)];
			int node_5 = hNodeType[idxScalarBlock(xp1 % BLOCK_NX, yp1 % BLOCK_NY, xp1 / BLOCK_NX, yp1 / BLOCK_NY)];
			int node_6 = hNodeType[idxScalarBlock(xm1 % BLOCK_NX, yp1 % BLOCK_NY, xm1 / BLOCK_NX, yp1 / BLOCK_NY)];
			int node_7 = hNodeType[idxScalarBlock(xm1 % BLOCK_NX, ym1 % BLOCK_NY, xm1 / BLOCK_NX, ym1 / BLOCK_NY)];
			int node_8 = hNodeType[idxScalarBlock(xp1 % BLOCK_NX, ym1 % BLOCK_NY, xp1 / BLOCK_NX, ym1 / BLOCK_NY)];

			bool anyBulk =
				node_1 == BULK ||
				node_2 == BULK ||
				node_3 == BULK ||
				node_4 == BULK ||
				node_5 == BULK ||
				node_6 == BULK ||
				node_7 == BULK ||
				node_8 == BULK;

			if (node_0 == SOLID_NODE && anyBulk)
			{
				if (count >= max_count)
				{
					max_count++;
					cylinderProperties *new_cylinder_properties;
					cudaMallocHost((void **)&new_cylinder_properties, sizeof(cylinderProperties) * max_count);
					memcpy(new_cylinder_properties, *cylinder_properties, sizeof(cylinderProperties) * count);
					cudaFreeHost(*cylinder_properties);
					*cylinder_properties = new_cylinder_properties;
				}

				(*cylinder_properties)[count].xb = x;
				(*cylinder_properties)[count].yb = y;

				// Calculate radius
				float x_diff = x - xc;
				float y_diff = y - yc;

				float r2 = x_diff * x_diff + y_diff * y_diff;

				float lattice_radius = sqrt(r2);

				if (lattice_radius > max_radius)
				{
					max_radius = lattice_radius;
				}

				count++;

				int bit_1 = node_3 != BULK && node_4 != BULK && node_7 != BULK ? 0 : 1;
				int bit_2 = node_1 != BULK && node_4 != BULK && node_8 != BULK ? 0 : 1;
				int bit_4 = node_2 != BULK && node_3 != BULK && node_6 != BULK ? 0 : 1;
				int bit_8 = node_1 != BULK && node_2 != BULK && node_5 != BULK ? 0 : 1;

				int bc_number = bit_1 * 1 + bit_2 * 2 + bit_4 * 4 + bit_8 * 8;

				hNodeType[idxScalarBlock(x % BLOCK_NX, y % BLOCK_NY, x / BLOCK_NX, y / BLOCK_NY)] = bc_number + 100;
			}
		}
	}

	for (int i = 0; i < count; i++)
	{
		int xb = (*cylinder_properties)[i].xb;
		int yb = (*cylinder_properties)[i].yb;

		// Calculate radius
		dfloat xb_diff = xb - xc;
		dfloat yb_diff = yb - yc;

		dfloat r2 = xb_diff * xb_diff + yb_diff * yb_diff;

		dfloat lattice_radius = sqrt(r2);

		dfloat inv_radius = 1.0 / lattice_radius;
		// ----------------------------

		// location at wall
		dfloat unit_nx = xb_diff * inv_radius;
		dfloat unit_ny = yb_diff * inv_radius;

		dfloat xw = xc + max_radius * unit_nx;
		dfloat yw = yc + max_radius * unit_ny;

		(*cylinder_properties)[i].xw = xw;
		(*cylinder_properties)[i].yw = yw;

		dfloat xwb_diff = xw - xb;
		dfloat ywb_diff = yw - yb;

		dfloat dr2 = xwb_diff * xwb_diff + ywb_diff * ywb_diff;
		(*cylinder_properties)[i].dr = sqrt(dr2);
		// ---------------------------

		// fluid point 1
		dfloat delx = sqrt(2.0);
		dfloat x1 = xw + delx * unit_nx;
		dfloat y1 = yw + delx * unit_ny;

		(*cylinder_properties)[i].x1 = x1;
		(*cylinder_properties)[i].y1 = y1;
		// ---------------------------

		// fluid point 2
		dfloat x2 = xw + 2.0 * delx * unit_nx;
		dfloat y2 = yw + 2.0 * delx * unit_ny;

		(*cylinder_properties)[i].x2 = x2;
		(*cylinder_properties)[i].y2 = y2;
		// ---------------------------

		// fluid point 3
		dfloat x3 = xw + 3.0 * delx * unit_nx;
		dfloat y3 = yw + 3.0 * delx * unit_ny;

		(*cylinder_properties)[i].x3 = x3;
		(*cylinder_properties)[i].y3 = y3;
		// ---------------------------

		dfloat theta = atan2f(yb_diff, xb_diff);
		if (theta < 0)
		{
			theta += 2.0 * M_PI;
		}

		(*cylinder_properties)[i].theta = theta;

		nodeType = hNodeType[idxScalarBlock(xb % BLOCK_NX, yb % BLOCK_NY, xb / BLOCK_NX, yb / BLOCK_NY)];
		nodeType -= 100;

		(*cylinder_properties)[i].is[0] = 0;
		(*cylinder_properties)[i].is[1] = 0;
		(*cylinder_properties)[i].is[2] = 0;
		(*cylinder_properties)[i].is[3] = 0;
		(*cylinder_properties)[i].is[4] = 0;
		(*cylinder_properties)[i].is[5] = 0;
		(*cylinder_properties)[i].is[6] = 0;
		(*cylinder_properties)[i].is[7] = 0;
		(*cylinder_properties)[i].is[8] = 0;

		(*cylinder_properties)[i].os[0] = 0;
		(*cylinder_properties)[i].os[1] = 0;
		(*cylinder_properties)[i].os[2] = 0;
		(*cylinder_properties)[i].os[3] = 0;
		(*cylinder_properties)[i].os[4] = 0;
		(*cylinder_properties)[i].os[5] = 0;
		(*cylinder_properties)[i].os[6] = 0;
		(*cylinder_properties)[i].os[7] = 0;
		(*cylinder_properties)[i].os[8] = 0;

		if (nodeType & 1)
		{
			(*cylinder_properties)[i].is[0] = 1;
			(*cylinder_properties)[i].is[1] = 1;
			(*cylinder_properties)[i].is[2] = 1;
			(*cylinder_properties)[i].is[5] = 1;

			(*cylinder_properties)[i].os[0] = 1;
			(*cylinder_properties)[i].os[3] = 1;
			(*cylinder_properties)[i].os[4] = 1;
			(*cylinder_properties)[i].os[7] = 1;
		}

		if (nodeType & 2)
		{
			(*cylinder_properties)[i].is[0] = 1;
			(*cylinder_properties)[i].is[2] = 1;
			(*cylinder_properties)[i].is[3] = 1;
			(*cylinder_properties)[i].is[6] = 1;

			(*cylinder_properties)[i].os[0] = 1;
			(*cylinder_properties)[i].os[1] = 1;
			(*cylinder_properties)[i].os[4] = 1;
			(*cylinder_properties)[i].os[8] = 1;
		}

		if (nodeType & 4)
		{
			(*cylinder_properties)[i].is[0] = 1;
			(*cylinder_properties)[i].is[1] = 1;
			(*cylinder_properties)[i].is[4] = 1;
			(*cylinder_properties)[i].is[8] = 1;

			(*cylinder_properties)[i].os[0] = 1;
			(*cylinder_properties)[i].os[2] = 1;
			(*cylinder_properties)[i].os[3] = 1;
			(*cylinder_properties)[i].os[6] = 1;
		}

		if (nodeType & 8)
		{
			(*cylinder_properties)[i].is[0] = 1;
			(*cylinder_properties)[i].is[3] = 1;
			(*cylinder_properties)[i].is[4] = 1;
			(*cylinder_properties)[i].is[7] = 1;

			(*cylinder_properties)[i].os[0] = 1;
			(*cylinder_properties)[i].os[1] = 1;
			(*cylinder_properties)[i].os[2] = 1;
			(*cylinder_properties)[i].os[5] = 1;
		}
	}

	*D_MAX = (float)2 * max_radius;
	*contour_counter = count;
}