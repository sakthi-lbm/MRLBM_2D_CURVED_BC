#ifndef __ARRAYINDEX_H
#define __ARRAYINDEX_H

#include "var.h"
// #include "definitions.h"

constexpr int M_RHO_INDEX = 0;
constexpr int M_UX_INDEX = 1;
constexpr int M_UY_INDEX = 2;
constexpr int M_MXX_INDEX = 3;
constexpr int M_MXY_INDEX = 4;
constexpr int M_MYY_INDEX = 5;


#ifdef M_OFFSET
#undef M_OFFSET
#endif

#define M_OFFSET M_MYY_INDEX

const size_t NUMBER_MOMENTS = M_OFFSET + 1;

#endif