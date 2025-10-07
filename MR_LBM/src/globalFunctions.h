
#ifndef __GLOBAL_FUNCTIONS_H
#define __GLOBAL_FUNCTIONS_H

#include <builtin_types.h> // for device variables
#include "var.h"
#include "globalStructs.h"

__host__ __device__
    size_t __forceinline__
    idxMom(
        const int tx,
        const int ty,
        const int mom,
        const int bx,
        const int by)
{
    return tx + BLOCK_NX * (ty + BLOCK_NY * (mom + NUMBER_MOMENTS * (bx + NUM_BLOCK_X * by)));
}

__device__ int __forceinline__
idxPopX(
    const int ty,
    const int pop,
    const int bx,
    const int by)
{
    return ty + BLOCK_NY * (pop + QF * (bx + NUM_BLOCK_X * by));
}

__device__ int __forceinline__
idxPopY(
    const int tx,
    const int pop,
    const int bx,
    const int by)
{
    return tx + BLOCK_NX * (pop + QF * (bx + NUM_BLOCK_X * by));
}

__host__ __device__
    size_t __forceinline__
    idxScalarBlock(
        const int tx,
        const int ty,
        const int bx,
        const int by)
{
    return tx + BLOCK_NX * (ty + BLOCK_NY * (bx + NUM_BLOCK_X * by));
}

__host__ __device__
    size_t __forceinline__
    idxPopBlock(const unsigned int tx, const unsigned int ty, const unsigned int pop)
{
    return tx + BLOCK_NX * (ty + BLOCK_NY * (pop));
}

__host__ __device__
size_t __forceinline__
idxScalarGlobal(unsigned int x, unsigned int y)
{
    return x + NX * y;
}

__host__ __device__
size_t __forceinline__
idxCylinder(unsigned int x, unsigned int y)
{
    return x + NX * y;
}

#endif // !__GLOBAL_FUNCTIONS_H
