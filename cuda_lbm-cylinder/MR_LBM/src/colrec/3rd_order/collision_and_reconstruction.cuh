#ifndef COLLISION_AND_RECONSTRUCTION
#define COLLISION_AND_RECONSTRUCTION

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "../../var.h"

__device__
inline void pop_reconstruction(dfloat rhoVar, dfloat ux, dfloat uy, dfloat mxx, dfloat myy, dfloat mxy, dfloat* pop) {
	dfloat mxxy = ux * mxy + uy * mxx;
	dfloat mxyy = uy * mxy + ux * myy;

	dfloat one_cs2 = 1.0 - cs2;
	dfloat cs2_mxxy = cs2 * mxxy;
	dfloat cs2_mxyy = cs2 * mxyy;

	dfloat pics2 = 1 - cs2 * (mxx + myy);

	dfloat multiplyTerm = W0 * rhoVar;
	pop[0] = multiplyTerm * (pics2);

	multiplyTerm = W1 * rhoVar;
	pop[1] = multiplyTerm * (pics2 + ux + mxx - cs2_mxyy);
	pop[2] = multiplyTerm * (pics2 + uy + myy - cs2_mxxy);
	pop[3] = multiplyTerm * (pics2 - ux + mxx + cs2_mxyy);
	pop[4] = multiplyTerm * (pics2 - uy + myy + cs2_mxxy);

	multiplyTerm = W2 * rhoVar;
	pop[5] = multiplyTerm * (pics2 + ux + uy + mxx + myy + mxy + one_cs2 * mxxy + one_cs2 * mxyy);
	pop[6] = multiplyTerm * (pics2 - ux + uy + mxx + myy - mxy + one_cs2 * mxxy - one_cs2 * mxyy);
	pop[7] = multiplyTerm * (pics2 - ux - uy + mxx + myy + mxy - one_cs2 * mxxy - one_cs2 * mxyy);
	pop[8] = multiplyTerm * (pics2 + ux - uy + mxx + myy - mxy - one_cs2 * mxxy + one_cs2 * mxyy);
}

__device__
inline void moment_collision(dfloat ux, dfloat uy, dfloat* mxx, dfloat* myy, dfloat* mxy,
#ifdef CYLINDER
dfloat OMEGA
#endif
) {
	const dfloat omegaVar = OMEGA;
	const dfloat t_omegaVar = 1 - omegaVar;
	const dfloat omegaVar_d2 = omegaVar / 2;

	*mxx = (t_omegaVar * (*mxx) + omegaVar_d2 * ux * ux);
	*myy = (t_omegaVar * (*myy) + omegaVar_d2 * uy * uy);

	*mxy = (t_omegaVar * (*mxy) + omegaVar * ux * uy);
}

#endif // COLLISION_AND_RECONSTRUCTION
