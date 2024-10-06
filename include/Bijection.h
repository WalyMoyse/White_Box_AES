#ifndef BIJECTION_H
#define BIJECTION_H
#include <stdio.h>

extern unsigned int Mat32[9][4][32][32]; // 32*32 invertable matrix
extern unsigned int InvMat32[9][4][32][32]; // Their inverse
extern unsigned char Mat8[9][16][8][8]; // Same for 8*8 matrix
extern unsigned char InvMat8[9][16][8][8];


extern unsigned char Bij[10000][16]; // Will be used for the encoding part when functional
extern unsigned char InvBij[10000][16]; // Same

#endif /* BIJECTION_H */