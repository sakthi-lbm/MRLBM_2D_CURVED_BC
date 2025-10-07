#ifndef COLLISION_AND_RECONSTRUCTION
#define COLLISION_AND_RECONSTRUCTION

// CUDA INCLUDE
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include "../../var.h"

__device__
inline void pop_reconstruction(dfloat rhoVar, dfloat ux, dfloat uy, dfloat mxx, dfloat myy, dfloat mxy, dfloat* pop) {
	dfloat pics2 = 1 - cs2 * (mxx + myy);

	dfloat multiplyTerm = W0 * rhoVar;
	pop[0] = multiplyTerm * (pics2) - W0;

	multiplyTerm = W1 * rhoVar;
	pop[1] = multiplyTerm * (pics2 + ux + mxx) - W1;
	pop[2] = multiplyTerm * (pics2 + uy + myy) - W1;
	pop[3] = multiplyTerm * (pics2 - ux + mxx) - W1;
	pop[4] = multiplyTerm * (pics2 - uy + myy) - W1;

	multiplyTerm = W2 * rhoVar;
	pop[5] = multiplyTerm * (pics2 + ux + uy + mxx + myy + mxy) - W2;
	pop[6] = multiplyTerm * (pics2 - ux + uy + mxx + myy - mxy) - W2;
	pop[7] = multiplyTerm * (pics2 - ux - uy + mxx + myy + mxy) - W2;
	pop[8] = multiplyTerm * (pics2 + ux - uy + mxx + myy - mxy) - W2;
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
