#ifndef TOOLS_H
#define TOOLS_H
#include <stdio.h>
#define size_block 16

// Table-based AES construction (1st step)
void T_box_calculation(unsigned char* key, unsigned char T_box[10][16][256]);
int Mult_Vect(unsigned char c, unsigned char vect[4]);
void Tyi_tables(unsigned int tyi_tables[4][256]) ;
void construc_table(unsigned int TBoxesTyiTables[9][16][256], unsigned char tbox[16][256], unsigned char XOR[9][192][16][16], unsigned char* key);

//MixingTable Part (2nd step)
unsigned char Matrix8Mult(unsigned char mat[8][8], unsigned char X);
unsigned int Matrix32Mult(unsigned int mat[32][32], unsigned int X);
void construc_mixing_bijection(unsigned int mixin[10][16][256], unsigned int mixout[9][16][256], unsigned int TBoxesTyiTables[9][16][256],  unsigned char tbox[16][256]);

//Encoding part (Last step)
void concatenate_bij(unsigned char bij1[16], unsigned char bij2[16], unsigned char bijres[256]);
int apply_concatenate_bij(int currentbij, int t);
int apply_concatenate_invbij(int currentbij, int t);
void encode_tables(unsigned int mixin[10][16][256], unsigned int mixout[9][16][256], unsigned char XOR[9][192][16][16]);



#endif /* TOOLS_H */