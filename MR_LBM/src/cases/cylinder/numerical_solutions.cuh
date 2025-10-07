#ifndef NUMERICAL_SOLUTION_CUH
#define NUMERICAL_SOLUTION_CUH

#include "../../var.h"

__device__
inline void numericalSolution(
	dfloat* rhoVar, const dfloat ux, const dfloat uy, dfloat* mxx_I, dfloat* mxy_I, dfloat* myy_I,
	const int(&incomings)[9], const int(&outgoings)[9], const dfloat OMEGA)
{
	dfloat rho;

	dfloat A_prime = 0;
	dfloat E_prime = 0;

	dfloat B11_prime = 0;
	dfloat B22_prime = 0;
	dfloat B12_prime = 0;

	dfloat D_xx_prime = 0;
	dfloat D_yy_prime = 0;
	dfloat D_xy_prime = 0;

	dfloat F11_xx_prime = 0;
	dfloat F22_xx_prime = 0;
	dfloat F12_xx_prime = 0;

	dfloat F11_yy_prime = 0;
	dfloat F22_yy_prime = 0;
	dfloat F12_yy_prime = 0;

	dfloat F11_xy_prime = 0;
	dfloat F22_xy_prime = 0;
	dfloat F12_xy_prime = 0;

	dfloat G_prime;

	dfloat J11_prime;
	dfloat J22_prime;
	dfloat J12_prime;

	dfloat J11_xx_star;
	dfloat J22_xx_star;
	dfloat J12_xx_star;

	dfloat J11_yy_star;
	dfloat J22_yy_star;
	dfloat J12_yy_star;

	dfloat J11_xy_star;
	dfloat J22_xy_star;
	dfloat J12_xy_star;

	dfloat a1, a2, a3;
	dfloat b1, b2, b3;
	dfloat c1, c2, c3;
	dfloat d1, d2, d3;

	dfloat denominator;
	dfloat inv_denominator;

	for (int i = 0; i < 9; i++) {
		dfloat Hxx = cx[i] * cx[i] - cs2;
		dfloat Hyy = cy[i] * cy[i] - cs2;
		dfloat Hxy = cx[i] * cy[i];

		dfloat A_i = w[i] * (1.0 + as2 * ux * cx[i] + as2 * uy * cy[i]);

		dfloat B11_i = 4.5 * w[i] * Hxx;
		dfloat B22_i = 4.5 * w[i] * Hyy;
		dfloat B12_i = 4.5 * w[i] * Hxy;

		if (outgoings[i] == 1)
		{
			A_prime += A_i;
			E_prime += w[i] * (1.0 + as2 * ux * cx[i] + as2 * uy * cy[i] + 4.5 * ux * ux * Hxx + 4.5 * uy * uy * Hyy + 9.0 * ux * uy * Hxy);

			B11_prime += B11_i;
			B22_prime += B22_i;
			B12_prime += B12_i;
		}

		if (incomings[i] == 1)
		{
			D_xx_prime += A_i * Hxx;
			D_yy_prime += A_i * Hyy;
			D_xy_prime += A_i * Hxy;

			F11_xx_prime += B11_i * Hxx;
			F22_xx_prime += B22_i * Hxx;
			F12_xx_prime += B12_i * Hxx;

			F11_yy_prime += B11_i * Hyy;
			F22_yy_prime += B22_i * Hyy;
			F12_yy_prime += B12_i * Hyy;

			F11_xy_prime += B11_i * Hxy;
			F22_xy_prime += B22_i * Hxy;
			F12_xy_prime += B12_i * Hxy;
		}
	}

	G_prime = (1.0 - OMEGA) * A_prime + OMEGA * E_prime;

	J11_prime = (1.0 - OMEGA) * B11_prime;
	J22_prime = (1.0 - OMEGA) * B22_prime;
	J12_prime = (1.0 - OMEGA) * B12_prime;

	J11_xx_star = J11_prime * (*mxx_I);
	J22_xx_star = J22_prime * (*mxx_I);
	J12_xx_star = J12_prime * (*mxx_I);

	J11_yy_star = J11_prime * (*myy_I);
	J22_yy_star = J22_prime * (*myy_I);
	J12_yy_star = J12_prime * (*myy_I);

	J11_xy_star = J11_prime * (*mxy_I);
	J22_xy_star = J22_prime * (*mxy_I);
	J12_xy_star = J12_prime * (*mxy_I);

	a1 = J11_xx_star - F11_xx_prime;
	a2 = J11_yy_star - F11_yy_prime;
	a3 = J11_xy_star - F11_xy_prime;

	b1 = J22_xx_star - F22_xx_prime;
	b2 = J22_yy_star - F22_yy_prime;
	b3 = J22_xy_star - F22_xy_prime;

	c1 = 2.0 * (J12_xx_star - F12_xx_prime);
	c2 = 2.0 * (J12_yy_star - F12_yy_prime);
	c3 = 2.0 * (J12_xy_star - F12_xy_prime);

	d1 = D_xx_prime - G_prime * (*mxx_I);
	d2 = D_yy_prime - G_prime * (*myy_I);

	d3 = D_xy_prime - G_prime * (*mxy_I);
	denominator = a3 * b2 * c1 - a2 * b3 * c1 - a3 * b1 * c2 + a1 * b3 * c2 + a2 * b1 * c3 - a1 * b2 * c3;
	inv_denominator = 1.0 / denominator;

	*mxx_I = -(-b3 * c2 * d1 + b2 * c3 * d1 + b3 * c1 * d2 - b1 * c3 * d2 - b2 * c1 * d3 + b1 * c2 * d3) * inv_denominator;
	*myy_I = -(a3 * c2 * d1 - a2 * c3 * d1 - a3 * c1 * d2 + a1 * c3 * d2 + a2 * c1 * d3 - a1 * c2 * d3) * inv_denominator;
	*mxy_I = -(-a3 * b2 * d1 + a2 * b3 * d1 + a3 * b1 * d2 - a1 * b3 * d2 - a2 * b1 * d3 + a1 * b2 * d3) * inv_denominator;

	denominator = (1.0 - OMEGA) * A_prime + (1.0 - OMEGA) * ((*mxx_I) * B11_prime + (*myy_I) * B22_prime + 2.0 * (*mxy_I) * B12_prime) + OMEGA * E_prime;
	inv_denominator = 1.0 / denominator;

	rho = (*rhoVar) * inv_denominator;

	*rhoVar = rho;
}

__device__
inline void numericalSolution_rotation(
	dfloat* rhoVar, const dfloat ux, const dfloat uy, dfloat* mxx_I, dfloat* mxy_I, dfloat* myy_I,
	dfloat mxx, dfloat myy,
	const int(&incomings)[9], const int(&outgoings)[9], const dfloat OMEGA, int x, int y)
{
	dfloat A_prime = 0;
	dfloat E_prime = 0;

	dfloat B11_prime = 0;
	dfloat B22_prime = 0;
	dfloat B12_prime = 0;

	dfloat D_xy_prime = 0;

	dfloat F11_xy_prime = 0;
	dfloat F22_xy_prime = 0;
	dfloat F12_xy_prime = 0;

	float xc = (float)(L_front + D * 0.5);
	float yc = (float)(L_bot + D * 0.5);

	dfloat radius = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));

	dfloat cos_theta = (x - xc) / radius;
	dfloat sen_theta = (y - yc) / radius;
	dfloat sen_two_theta = 2.0 * cos_theta * sen_theta;
	dfloat cos_two_theta = cos_theta * cos_theta - sen_theta * sen_theta;

	for (int i = 0; i < 9; i++) {
		dfloat cx_p = cx[i] * cos_theta + cy[i] * sen_theta;
		dfloat cy_p = cy[i] * cos_theta - cx[i] * sen_theta;

		const dfloat Hxx_p = cx_p * cx_p - cs2;
		const dfloat Hyy_p = cy_p * cy_p - cs2;
		const dfloat Hxy_p = cx_p * cy_p;

		dfloat ux_p = ux * cos_theta + uy * sen_theta;
		dfloat uy_p = uy * cos_theta - ux * sen_theta;

		dfloat A_i = w[i] * (1.0 + as2 * ux_p * cx_p + as2 * uy_p * cy_p);

		dfloat B11_i = 4.5 * w[i] * Hxx_p;
		dfloat B22_i = 4.5 * w[i] * Hyy_p;
		dfloat B12_i = 4.5 * w[i] * Hxy_p;

		if (outgoings[i] == 1)
		{
			A_prime += A_i;
			E_prime += w[i] * (1.0 + as2 * ux_p * cx_p + as2 * uy_p * cy_p + 4.5 * ux_p * ux_p * Hxx_p + 4.5 * uy_p * uy_p * Hyy_p + 9.0 * ux_p * uy_p * Hxy_p);

			B11_prime += B11_i;
			B22_prime += B22_i;
			B12_prime += B12_i;
		}

		if (incomings[i] == 1)
		{
			D_xy_prime += A_i * Hxy_p;

			F11_xy_prime += B11_i * Hxy_p;
			F22_xy_prime += B22_i * Hxy_p;
			F12_xy_prime += B12_i * Hxy_p;
		}
	}

	const dfloat G_prime = (1.0 - OMEGA) * A_prime + OMEGA * E_prime;

	const dfloat J11_prime = (1.0 - OMEGA) * B11_prime;
	const dfloat J22_prime = (1.0 - OMEGA) * B22_prime;
	const dfloat J12_prime = (1.0 - OMEGA) * B12_prime;

	const dfloat J11_xy_star = J11_prime * (*mxy_I);
	const dfloat J22_xy_star = J22_prime * (*mxy_I);
	const dfloat J12_xy_star = J12_prime * (*mxy_I);

	const dfloat L_11_xy = J11_xy_star - F11_xy_prime;

	const dfloat L_22_xy = J22_xy_star - F22_xy_prime;

	const dfloat L_12_xy = 2.0 * (J12_xy_star - F12_xy_prime);

	const dfloat R_xy = D_xy_prime - G_prime * (*mxy_I);

	dfloat mxx_p = mxx;
	dfloat myy_p = myy;
	dfloat mxy_p = (R_xy - mxx_p * L_11_xy - myy_p * L_22_xy) / L_12_xy;

	const dfloat denominator = (1.0 - OMEGA) * A_prime + (1.0 - OMEGA) * (mxx_p * B11_prime + myy_p * B22_prime + 2.0 * mxy_p * B12_prime) + OMEGA * E_prime;
	const dfloat inv_denominator = 1.0 / denominator;

	const dfloat rho = (*rhoVar) * inv_denominator;

	*rhoVar = rho;
	*mxx_I = mxx_p * cos_theta * cos_theta + myy_p * sen_theta * sen_theta - mxy_p * sen_two_theta;
	*myy_I = mxx_p * sen_theta * sen_theta + myy_p * cos_theta * cos_theta + mxy_p * sen_two_theta;
	*mxy_I = (mxx_p - myy_p) * 0.5 * sen_two_theta + mxy_p * cos_two_theta;
}

__device__
inline void numericalSolution_rotation_rhoeq(
	dfloat* rhoVar, const dfloat ux, const dfloat uy, dfloat* mxx_I, dfloat* mxy_I, dfloat* myy_I,
	dfloat mxx, dfloat myy,
	const int(&incomings)[9], const int(&outgoings)[9], const dfloat OMEGA, int x, int y)
{
	dfloat A_prime = 0;
	dfloat E_prime = 0;

	dfloat B11_prime = 0;
	dfloat B22_prime = 0;
	dfloat B12_prime = 0;

	dfloat D_xy_prime = 0;

	dfloat F11_xy_prime = 0;
	dfloat F22_xy_prime = 0;
	dfloat F12_xy_prime = 0;

	float xc = (float)(L_front + D * 0.5);
	float yc = (float)(L_bot + D * 0.5);

	dfloat radius = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));

	dfloat cos_theta = (x - xc) / radius;
	dfloat sen_theta = (y - yc) / radius;
	dfloat sen_two_theta = 2.0 * cos_theta * sen_theta;
	dfloat cos_two_theta = cos_theta * cos_theta - sen_theta * sen_theta;

	for (int i = 0; i < 9; i++) {
		dfloat cx_p = cx[i] * cos_theta + cy[i] * sen_theta;
		dfloat cy_p = cy[i] * cos_theta - cx[i] * sen_theta;

		const dfloat Hxx_p = cx_p * cx_p - cs2;
		const dfloat Hyy_p = cy_p * cy_p - cs2;
		const dfloat Hxy_p = cx_p * cy_p;

		dfloat ux_p = ux * cos_theta + uy * sen_theta;
		dfloat uy_p = uy * cos_theta - ux * sen_theta;

		dfloat A_i = w[i] * (1.0 + as2 * ux_p * cx_p + as2 * uy_p * cy_p);

		dfloat B11_i = 4.5 * w[i] * Hxx_p;
		dfloat B22_i = 4.5 * w[i] * Hyy_p;
		dfloat B12_i = 4.5 * w[i] * Hxy_p;

		if (outgoings[i] == 1)
		{
			A_prime += A_i;
			E_prime += w[i] * (1.0 + as2 * ux_p * cx_p + as2 * uy_p * cy_p + 4.5 * ux_p * ux_p * Hxx_p + 4.5 * uy_p * uy_p * Hyy_p + 9.0 * ux_p * uy_p * Hxy_p);

			B11_prime += B11_i;
			B22_prime += B22_i;
			B12_prime += B12_i;
		}

		if (incomings[i] == 1)
		{
			D_xy_prime += A_i * Hxy_p;

			F11_xy_prime += B11_i * Hxy_p;
			F22_xy_prime += B22_i * Hxy_p;
			F12_xy_prime += B12_i * Hxy_p;
		}
	}

	const dfloat L_11_xy = F11_xy_prime;
	const dfloat L_22_xy = F22_xy_prime;
	const dfloat L_12_xy = 2.0 * F12_xy_prime;

	const dfloat R_xy = E_prime * (*mxy_I) - D_xy_prime ;

	dfloat mxx_p = mxx;
	dfloat myy_p = myy;
	dfloat mxy_p = (R_xy - mxx_p * L_11_xy - myy_p * L_22_xy) / L_12_xy;

	const dfloat denominator = E_prime;
	const dfloat inv_denominator = 1.0 / denominator;

	const dfloat rho = (*rhoVar) * inv_denominator;

	*rhoVar = rho;
	*mxx_I = mxx_p * cos_theta * cos_theta + myy_p * sen_theta * sen_theta - mxy_p * sen_two_theta;
	*myy_I = mxx_p * sen_theta * sen_theta + myy_p * cos_theta * cos_theta + mxy_p * sen_two_theta;
	*mxy_I = (mxx_p - myy_p) * 0.50 * sen_two_theta + mxy_p * cos_two_theta;
}


#endif