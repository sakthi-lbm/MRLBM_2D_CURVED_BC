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