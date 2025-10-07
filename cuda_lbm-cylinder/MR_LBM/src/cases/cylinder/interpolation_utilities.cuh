#ifndef INTERPOLATION_UTILITIES_CUH
#define INTERPOLATION_UTILITIES_CUH

#include "../../var.h"
#include "../../globalFunctions.h"

__device__
inline void bilinear_density_interpolation(dfloat x, dfloat y, int x0, int y0, int x1, int y1, dfloat* moms, int mom_index, dfloat* r_f) {
	const dfloat xd = x - x0;
	const dfloat yd = y - y0;

	const dfloat r1 = RHO_0 + moms[idxMom(x0 % BLOCK_NX, y0 % BLOCK_NY, mom_index, x0 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat r2 = RHO_0 + moms[idxMom(x1 % BLOCK_NX, y0 % BLOCK_NY, mom_index, x1 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat r3 = RHO_0 + moms[idxMom(x0 % BLOCK_NX, y1 % BLOCK_NY, mom_index, x0 / BLOCK_NX, y1 / BLOCK_NY)];
	const dfloat r4 = RHO_0 + moms[idxMom(x1 % BLOCK_NX, y1 % BLOCK_NY, mom_index, x1 / BLOCK_NX, y1 / BLOCK_NY)];

	const dfloat r_y0 = (1.0 - xd) * r1 + xd * r2;
	const dfloat r_y1 = (1.0 - xd) * r3 + xd * r4;

	*r_f = (1.0 - yd) * r_y0 + yd * r_y1;
}

__device__
inline void bilinear_velocity_interpolation(dfloat x, dfloat y, int x0, int y0, int x1, int y1, dfloat* moms, int mom_index, dfloat* u_f) {
	const dfloat xd = x - x0;
	const dfloat yd = y - y0;

	const dfloat ux1 = moms[idxMom(x0 % BLOCK_NX, y0 % BLOCK_NY, mom_index, x0 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat ux2 = moms[idxMom(x1 % BLOCK_NX, y0 % BLOCK_NY, mom_index, x1 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat ux3 = moms[idxMom(x0 % BLOCK_NX, y1 % BLOCK_NY, mom_index, x0 / BLOCK_NX, y1 / BLOCK_NY)];
	const dfloat ux4 = moms[idxMom(x1 % BLOCK_NX, y1 % BLOCK_NY, mom_index, x1 / BLOCK_NX, y1 / BLOCK_NY)];

	const dfloat ux_y0 = (1.0 - xd) * ux1 + xd * ux2;
	const dfloat ux_y1 = (1.0 - xd) * ux3 + xd * ux4;

	*u_f = (1.0 - yd) * ux_y0 + yd * ux_y1;
}

__device__
inline void bilinear_moment_interpolation(dfloat x, dfloat y, int x0, int y0, int x1, int y1, dfloat* moms, dfloat* mxx_f, dfloat* myy_f) {

	const dfloat xd = x - dfloat(x0);
	const dfloat yd = y - dfloat(y0);

	const dfloat radius_1 = sqrt((dfloat(x0) - xc) * (dfloat(x0) - xc) + (dfloat(y0) - yc) * (dfloat(y0) - yc));
	const dfloat radius_2 = sqrt((dfloat(x1) - xc) * (dfloat(x1) - xc) + (dfloat(y0) - yc) * (dfloat(y0) - yc));
	const dfloat radius_3 = sqrt((dfloat(x0) - xc) * (dfloat(x0) - xc) + (dfloat(y1) - yc) * (dfloat(y1) - yc));
	const dfloat radius_4 = sqrt((dfloat(x1) - xc) * (dfloat(x1) - xc) + (dfloat(y1) - yc) * (dfloat(y1) - yc));

	const dfloat cos_theta_1 = (dfloat(x0) - xc) / radius_1;
	const dfloat cos_theta_2 = (dfloat(x1) - xc) / radius_2;
	const dfloat cos_theta_3 = (dfloat(x0) - xc) / radius_3;
	const dfloat cos_theta_4 = (dfloat(x1) - xc) / radius_4;

	const dfloat sen_theta_1 = (dfloat(y0) - yc) / radius_1;
	const dfloat sen_theta_2 = (dfloat(y0) - yc) / radius_2;
	const dfloat sen_theta_3 = (dfloat(y1) - yc) / radius_3;
	const dfloat sen_theta_4 = (dfloat(y1) - yc) / radius_4;

	const dfloat sen_two_theta_1 = 2.0 * cos_theta_1 * sen_theta_1;
	const dfloat sen_two_theta_2 = 2.0 * cos_theta_2 * sen_theta_2;
	const dfloat sen_two_theta_3 = 2.0 * cos_theta_3 * sen_theta_3;
	const dfloat sen_two_theta_4 = 2.0 * cos_theta_4 * sen_theta_4;

	const dfloat mxx1 = moms[idxMom(x0 % BLOCK_NX, y0 % BLOCK_NY, M_MXX_INDEX, x0 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat mxx2 = moms[idxMom(x1 % BLOCK_NX, y0 % BLOCK_NY, M_MXX_INDEX, x1 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat mxx3 = moms[idxMom(x0 % BLOCK_NX, y1 % BLOCK_NY, M_MXX_INDEX, x0 / BLOCK_NX, y1 / BLOCK_NY)];
	const dfloat mxx4 = moms[idxMom(x1 % BLOCK_NX, y1 % BLOCK_NY, M_MXX_INDEX, x1 / BLOCK_NX, y1 / BLOCK_NY)];

	const dfloat myy1 = moms[idxMom(x0 % BLOCK_NX, y0 % BLOCK_NY, M_MYY_INDEX, x0 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat myy2 = moms[idxMom(x1 % BLOCK_NX, y0 % BLOCK_NY, M_MYY_INDEX, x1 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat myy3 = moms[idxMom(x0 % BLOCK_NX, y1 % BLOCK_NY, M_MYY_INDEX, x0 / BLOCK_NX, y1 / BLOCK_NY)];
	const dfloat myy4 = moms[idxMom(x1 % BLOCK_NX, y1 % BLOCK_NY, M_MYY_INDEX, x1 / BLOCK_NX, y1 / BLOCK_NY)];

	const dfloat mxy1 = moms[idxMom(x0 % BLOCK_NX, y0 % BLOCK_NY, M_MXY_INDEX, x0 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat mxy2 = moms[idxMom(x1 % BLOCK_NX, y0 % BLOCK_NY, M_MXY_INDEX, x1 / BLOCK_NX, y0 / BLOCK_NY)];
	const dfloat mxy3 = moms[idxMom(x0 % BLOCK_NX, y1 % BLOCK_NY, M_MXY_INDEX, x0 / BLOCK_NX, y1 / BLOCK_NY)];
	const dfloat mxy4 = moms[idxMom(x1 % BLOCK_NX, y1 % BLOCK_NY, M_MXY_INDEX, x1 / BLOCK_NX, y1 / BLOCK_NY)];

	const dfloat mxx_p1 = mxx1 * cos_theta_1 * cos_theta_1 + myy1 * sen_theta_1 * sen_theta_1 + mxy1 * sen_two_theta_1;
	const dfloat myy_p1 = mxx1 * sen_theta_1 * sen_theta_1 + myy1 * cos_theta_1 * cos_theta_1 - mxy1 * sen_two_theta_1;

	const dfloat mxx_p2 = mxx2 * cos_theta_2 * cos_theta_2 + myy2 * sen_theta_2 * sen_theta_2 + mxy2 * sen_two_theta_2;
	const dfloat myy_p2 = mxx2 * sen_theta_2 * sen_theta_2 + myy2 * cos_theta_2 * cos_theta_2 - mxy2 * sen_two_theta_2;

	const dfloat mxx_p3 = mxx3 * cos_theta_3 * cos_theta_3 + myy3 * sen_theta_3 * sen_theta_3 + mxy3 * sen_two_theta_3;
	const dfloat myy_p3 = mxx3 * sen_theta_3 * sen_theta_3 + myy3 * cos_theta_3 * cos_theta_3 - mxy3 * sen_two_theta_3;

	const dfloat mxx_p4 = mxx4 * cos_theta_4 * cos_theta_4 + myy4 * sen_theta_4 * sen_theta_4 + mxy4 * sen_two_theta_4;
	const dfloat myy_p4 = mxx4 * sen_theta_4 * sen_theta_4 + myy4 * cos_theta_4 * cos_theta_4 - mxy4 * sen_two_theta_4;

	const dfloat mxx_temp = (1.0 - xd) * mxx_p1 + xd * mxx_p2;
	const dfloat mxx_temp2 = (1.0 - xd) * mxx_p3 + xd * mxx_p4;

	const dfloat myy_temp = (1.0 - xd) * myy_p1 + xd * myy_p2;
	const dfloat myy_temp2 = (1.0 - xd) * myy_p3 + xd * myy_p4;

	*mxx_f = (1.0 - yd) * mxx_temp + yd * mxx_temp2;
	*myy_f = (1.0 - yd) * myy_temp + yd * myy_temp2;
}

__device__
inline void pressure_extrapolation(dfloat xw, dfloat yw, dfloat x1, dfloat y1, dfloat x2, dfloat y2, dfloat x3, dfloat y3, dfloat rho1, dfloat rho2, dfloat rho3, dfloat* pressure) {

	// pressure interpolation
	dfloat xw_diff = xw - xc;
	dfloat yw_diff = yw - yc;

	dfloat x1_diff = x1 - xc;
	dfloat y1_diff = y1 - yc;

	dfloat x2_diff = x2 - xc;
	dfloat y2_diff = y2 - yc;

	dfloat x3_diff = x3 - xc;
	dfloat y3_diff = y3 - yc;

	dfloat rw2 = xw_diff * xw_diff + yw_diff * yw_diff;
	dfloat r12 = x1_diff * x1_diff + y1_diff * y1_diff;
	dfloat r22 = x2_diff * x2_diff + y2_diff * y2_diff;
	dfloat r32 = x3_diff * x3_diff + y3_diff * y3_diff;

	dfloat rw = sqrt(rw2);
	dfloat r1 = sqrt(r12);
	dfloat r2 = sqrt(r22);
	dfloat r3 = sqrt(r32);

	dfloat denom = (r1 - r2) * (r1 - r3) * (r2 - r3);

	dfloat p1 = rho1 * cs2;
	dfloat p2 = rho2 * cs2;
	dfloat p3 = rho3 * cs2;

	dfloat a0 = (r1 * r3 * p2 * (r3 - r1) + (r2 * r2) * (r3 * p1 - r1 * p3) + r2 * ((r1 * r1) * p3 - (r3 * r3) * p1)) / denom;
	dfloat a1 = ((r3 * r3) * (p1 - p2) + (r1 * r1) * (p2 - p3) + (r2 * r2) * (p3 - p1)) / denom;
	dfloat a2 = (r3 * (p2 - p1) + r2 * (p1 - p3) + r1 * (p3 - p2)) / denom;

	*pressure = a0 + a1 * rw + a2 * (rw * rw);
}

__device__
inline dfloat extrapolation(dfloat delta, dfloat value1, dfloat value2) {
	const dfloat deltax = sqrt(2.0);

	return (delta * (delta - 2.0 * deltax) / (deltax * deltax)) * value1 - (delta * (delta - deltax) / (2.0 * deltax * deltax)) * value2;
}

#endif