#ifndef INTERFACE_HANDLING_CUH
#define INTERFACE_HANDLING_CUH

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "../var.h"
#include "../globalStructs.h"
#include "../globalFunctions.h"
#include "interface.h"

__device__
inline void pop_load(ghostInterfaceData ghostInterface, size_t tx, size_t ty, size_t bx, size_t by, dfloat* pop) {
    const int txm1 = (tx - 1 + BLOCK_NX) % BLOCK_NX;
    const int txp1 = (tx + 1 + BLOCK_NX) % BLOCK_NX;

    const int tym1 = (ty - 1 + BLOCK_NY) % BLOCK_NY;
    const int typ1 = (ty + 1 + BLOCK_NY) % BLOCK_NY;

    const int bxm1 = (bx - 1 + NUM_BLOCK_X) % NUM_BLOCK_X;
    const int bxp1 = (bx + 1 + NUM_BLOCK_X) % NUM_BLOCK_X;

    const int bym1 = (by - 1 + NUM_BLOCK_Y) % NUM_BLOCK_Y;
    const int byp1 = (by + 1 + NUM_BLOCK_Y) % NUM_BLOCK_Y;

    if (tx == 0)
    { // w
        pop[1] = ghostInterface.fGhost.X_1[idxPopX(ty, 0, bxm1, by)];
        pop[5] = ghostInterface.fGhost.X_1[idxPopX(tym1, 1, bxm1, ((ty == 0) ? bym1 : by))];
        pop[8] = ghostInterface.fGhost.X_1[idxPopX(typ1, 2, bxm1, ((ty == BLOCK_NY - 1) ? byp1 : by))];
    }

    else if (tx == (BLOCK_NX - 1))
    { // e
        pop[3] = ghostInterface.fGhost.X_0[idxPopX(ty, 0, bxp1, by)];
        pop[6] = ghostInterface.fGhost.X_0[idxPopX(tym1, 1, bxp1, ((ty == 0) ? bym1 : by))];
        pop[7] = ghostInterface.fGhost.X_0[idxPopX(typ1, 2, bxp1, ((ty == BLOCK_NY - 1) ? byp1 : by))];
    }

    if (ty == 0)
    { // s
        pop[2] = ghostInterface.fGhost.Y_1[idxPopY(tx, 0, bx, bym1)];
        pop[5] = ghostInterface.fGhost.Y_1[idxPopY(txm1, 1, ((tx == 0) ? bxm1 : bx), bym1)];
        pop[6] = ghostInterface.fGhost.Y_1[idxPopY(txp1, 2, ((tx == (BLOCK_NX - 1)) ? bxp1 : bx), bym1)];
    }
    else if (ty == (BLOCK_NY - 1))
    { // n
        pop[4] = ghostInterface.fGhost.Y_0[idxPopY(tx, 0, bx, byp1)];
        pop[7] = ghostInterface.fGhost.Y_0[idxPopY(txp1, 1, ((tx == (BLOCK_NX - 1)) ? bxp1 : bx), byp1)];
        pop[8] = ghostInterface.fGhost.Y_0[idxPopY(txm1, 2, ((tx == 0) ? bxm1 : bx), byp1)];
    }
}

__device__
inline void pop_save(ghostInterfaceData ghostInterface, size_t tx, size_t ty, size_t bx, size_t by, unsigned int x, unsigned int y, dfloat* pop) {
    /* write to global pop */
    if (INTERFACE_BC_WEST)
    { // w
        ghostInterface.gGhost.X_0[idxPopX(ty, 0, bx, by)] = pop[3];
        ghostInterface.gGhost.X_0[idxPopX(ty, 1, bx, by)] = pop[6];
        ghostInterface.gGhost.X_0[idxPopX(ty, 2, bx, by)] = pop[7];
    }

    if (INTERFACE_BC_EAST)
    { // e
        ghostInterface.gGhost.X_1[idxPopX(ty, 0, bx, by)] = pop[1];
        ghostInterface.gGhost.X_1[idxPopX(ty, 1, bx, by)] = pop[5];
        ghostInterface.gGhost.X_1[idxPopX(ty, 2, bx, by)] = pop[8];
    }

    if (INTERFACE_BC_SOUTH)
    { // s
        ghostInterface.gGhost.Y_0[idxPopY(tx, 0, bx, by)] = pop[4];
        ghostInterface.gGhost.Y_0[idxPopY(tx, 1, bx, by)] = pop[7];
        ghostInterface.gGhost.Y_0[idxPopY(tx, 2, bx, by)] = pop[8];
    }

    if (INTERFACE_BC_NORTH)
    { // n
        ghostInterface.gGhost.Y_1[idxPopY(tx, 0, bx, by)] = pop[2];
        ghostInterface.gGhost.Y_1[idxPopY(tx, 1, bx, by)] = pop[5];
        ghostInterface.gGhost.Y_1[idxPopY(tx, 2, bx, by)] = pop[6];
    }
}

#endif // !POP_LOAD_CUH