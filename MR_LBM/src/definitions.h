#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

#include "var.h"

/* --------------------------- CONSTANTS --------------------------- */

#define SQRT_2 (1.41421356237309504880168872420969807856967187537)
#define SQRT_10 (3.162277660168379331998893544432718533719555139325)

constexpr dfloat ONESIXTH = 1.0 / 6.0;
constexpr dfloat ONETHIRD = 1.0 / 3.0;

/* --------------------------- AUXILIARY DEFINES --------------------------- */
#define IN_HOST 1    // variable accessible only for host
#define IN_VIRTUAL 2 // variable accessible for device and host

constexpr size_t BYTES_PER_GB = (1 << 30);
constexpr size_t BYTES_PER_MB = (1 << 20);

/* ------------------------------ VELOCITY SET ------------------------------ */
constexpr unsigned char Q = 9;  // number of velocities
constexpr unsigned char QF = 3; // number of velocities on each face
constexpr dfloat W0 = 4.0 / 9;  // population 0 weight (0, 0, 0)
constexpr dfloat W1 = 1.0 / 9;  // adjacent populations (1, 0, 0)
constexpr dfloat W2 = 1.0 / 36; // diagonal populations (1, 1, 0)

// velocities weight vector
__device__ const dfloat w[Q] = { W0,
								W1, W1, W1, W1,
								W2, W2, W2, W2 };

constexpr dfloat as2 = 3.0;
constexpr dfloat cs2 = 1.0 / as2;

// populations velocities      0  1  2  3  4  5  6  7  8
__device__ constexpr dfloat cx[Q] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
__device__ constexpr dfloat cy[Q] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };

constexpr dfloat F_M_0_SCALE = 1.0;
constexpr dfloat F_M_I_SCALE = as2;
constexpr dfloat F_M_II_SCALE = as2 * as2 / 2;
constexpr dfloat F_M_IJ_SCALE = as2 * as2;

/* ------------------------------ MEMORY SIZE ------------------------------ */
#include "arrayIndex.h"

constexpr int SHARED_MEMORY_ELEMENT_SIZE = sizeof(dfloat) * (Q - 1);
 constexpr int MAX_ELEMENTS_IN_BLOCK = 48128 / SHARED_MEMORY_ELEMENT_SIZE;

constexpr BlockDim optimalBlockDimArray = findOptimalBlockDimensions(MAX_ELEMENTS_IN_BLOCK);

const int BLOCK_NX = optimalBlockDimArray.x; // number of threads in x
const int BLOCK_NY = optimalBlockDimArray.y; // number of threads in y

// const int BLOCK_NX = 16; // number of threads in x
// const int BLOCK_NY = 16; // number of threads in y

#define BLOCK_LBM_SIZE (BLOCK_NX * BLOCK_NY) // size of a block

const size_t BLOCK_GHOST_SIZE = BLOCK_NX + BLOCK_NY;

const size_t BLOCK_SIZE = BLOCK_LBM_SIZE + BLOCK_GHOST_SIZE;

const size_t X_CORRECTION = NX % BLOCK_NX > 0.0 ? 1 : 0;
const size_t Y_CORRECTION = NY % BLOCK_NY > 0.0 ? 1 : 0;

const size_t NUM_BLOCK_X = NX / BLOCK_NX + X_CORRECTION;
const size_t NUM_BLOCK_Y = NY / BLOCK_NY + Y_CORRECTION;


const size_t NUM_BLOCK = NUM_BLOCK_X * NUM_BLOCK_Y;

const size_t NUMBER_LBM_NODES = NUM_BLOCK * BLOCK_LBM_SIZE;
const size_t NUMBER_GHOST_FACE_X = BLOCK_NY * NUM_BLOCK_X * NUM_BLOCK_Y;
const size_t NUMBER_GHOST_FACE_Y = BLOCK_NX * NUM_BLOCK_X * NUM_BLOCK_Y;

const size_t MEM_SIZE_BLOCK_LBM = sizeof(dfloat) * BLOCK_LBM_SIZE * NUMBER_MOMENTS;
const size_t MEM_SIZE_BLOCK_GHOST = sizeof(dfloat) * BLOCK_GHOST_SIZE * Q;
const size_t MEM_SIZE_BLOCK_TOTAL = MEM_SIZE_BLOCK_GHOST + MEM_SIZE_BLOCK_LBM;

const size_t NUMBER_LBM_POP_NODES = NX * NY;

// memory size
const size_t MEM_SIZE_SCALAR = sizeof(dfloat) * NUMBER_LBM_POP_NODES;
const size_t MEM_SIZE_POP = sizeof(dfloat) * NUMBER_LBM_POP_NODES * Q;
const size_t MEM_SIZE_MOM = sizeof(dfloat) * NUMBER_LBM_NODES * NUMBER_MOMENTS;

#endif //!__DEFINITIONS_H