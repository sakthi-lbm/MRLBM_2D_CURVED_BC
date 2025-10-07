#ifndef __NODE_TYPE_MAP_H
#define __NODE_TYPE_MAP_H

#include <builtin_types.h>
#include <stdint.h>

// DIRECTION DEFINES 00000000

#define BULK (15)
// FACE
#define NORTH (3)
#define SOUTH (12)
#define WEST (10)
#define EAST (5) 
// CORNER
#define NORTH_WEST (2)
#define NORTH_EAST (1)
#define SOUTH_WEST (8)
#define SOUTH_EAST (4)
// IMMERSED
#define IMMERSED_01 (101)
#define IMMERSED_02 (102)
#define IMMERSED_03 (103)
#define IMMERSED_04 (104)
#define IMMERSED_05 (105)
#define IMMERSED_06 (106)
#define IMMERSED_07 (107)
#define IMMERSED_08 (108)
#define IMMERSED_09 (109)
#define IMMERSED_10 (110)
#define IMMERSED_11 (111)
#define IMMERSED_12 (112)
#define IMMERSED_13 (113)
#define IMMERSED_14 (114)
#define IMMERSED_15 (115)

#define SOLID_NODE (0)

#define MISSING_DEFINITION (0b11111111111111111111111111111111)

#define DIRECTION_BITS (0b11111 << DIRECTION_OFFSET)

#endif // !__NODE_TYPE_MAP
