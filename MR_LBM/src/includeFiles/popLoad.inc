
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
    pop[8] = ghostInterface.fGhost.Y_0[idxPopY(txm1, 2,  ((tx == 0) ? bxm1 : bx), byp1)];
}