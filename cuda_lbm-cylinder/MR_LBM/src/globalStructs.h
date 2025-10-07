#ifndef __GLOBAL_STRUCTS_H
#define __GLOBAL_STRUCTS_H

#include "var.h"
#include "errorDef.h"

typedef struct ghostData
{
    dfloat *X_0;
    dfloat *X_1;
    dfloat *Y_0;
    dfloat *Y_1;
} GhostData;

typedef struct ghostInterfaceData
{
    ghostData fGhost;
    ghostData gGhost;
    ghostData h_fGhost;
} GhostInterfaceData;

typedef struct cylinderProperties
{
    // x coordinates
    dfloat xb;
    dfloat xw;
    dfloat x1;
    dfloat x2;
    dfloat x3;

    // y coordinates
    dfloat yb;
    dfloat yw;
    dfloat y1;
    dfloat y2;
    dfloat y3;

    // incomings and outgoings
    int is[9];
    int os[9];

    // other properties
    double Fx;
    double Fy;
    dfloat ps;
    dfloat dr;
    dfloat theta;
} CylinderProperties;

#endif //__GLOBAL_STRUCTS
