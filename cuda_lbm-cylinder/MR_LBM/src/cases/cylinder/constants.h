#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../../var.h"

constexpr dfloat RE = 100;

constexpr int SCALE = 1;

constexpr int D = 32;       //Diameter of the cylinder

constexpr int HD = 10 * D;    //Height of the Domain !!!!should be a odd number!!!!
constexpr int LD = 60 * D;    //Length of the Domain
constexpr int L_front = 15 * D;

constexpr int L_bot = (HD - D) / 2;
constexpr int L_top = L_bot;

constexpr int L_back = LD - L_front - D;

constexpr int N = 1 * SCALE;

constexpr int NX = LD; // size x of the grid
constexpr int NY = HD;      // size y of the grid

constexpr dfloat xc = (dfloat)(L_front + D * 0.5);
constexpr dfloat yc = (dfloat)(L_bot + D * 0.5);

constexpr dfloat U_MAX = 0.1;
constexpr dfloat L = N;

constexpr int MACR_SAVE = 10 * D / U_MAX;
constexpr int tstar = 2000;         //non-dimensional time to start the statistics
constexpr int stat_period = 100;    //period of statistics * tstar 

constexpr int N_STEPS = (tstar + stat_period) * D / U_MAX;

// value for the velocity initial condition in the domain
constexpr dfloat U_0_X = 0.0;
constexpr dfloat U_0_Y = 0.0;
constexpr dfloat U_0_Z = 0.0;
constexpr dfloat RHO_0 = 1.0;

constexpr dfloat MACH_NUMBER = U_MAX / 0.57735026918962;

/* --------------------- INITIALIZATION LOADING DEFINES -------------------- */
constexpr int INI_STEP = 0; // initial simulation step (0 default)

#define BC_X_WALL
#define BC_Y_WALL
// #define BC_Y_PERIODIC

constexpr bool IRBC = false;
constexpr bool ROTATIONAL_COORDINATES = true;
constexpr bool RHO_STRONG = true;
constexpr bool RHO_EQ = false;

#define CYLINDER

#endif // !CONSTANTS_H