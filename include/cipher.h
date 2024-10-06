#ifndef CIPHER_H
#define CIPHER_H
#include <stdio.h>

// Classic AES
int GetPowerGenerator(unsigned char c);
unsigned char FastMult(unsigned char a, unsigned char b);
void State_mult(unsigned char* state, const unsigned char* mixcol);
unsigned char GetSboxValue(unsigned char c);
int get_nround(int key_size);
void SubBytes(unsigned char* state);
void ShiftRows(unsigned char* state);
void MixColumns(unsigned char* state);
void AddRoundKey(unsigned char* state, unsigned char* expandedkey, int step);
void cipher_block(unsigned char* state, unsigned char* expkey, int Nround);

// White-Box AES
int applyXOR(unsigned char XOR[16][16], int t1, int t2);
void WBC_AES_No_MixBij(unsigned char block[16], unsigned char tbox[16][256], unsigned int TBoxesTyiTables[9][16][256], unsigned char XOR[9][192][16][16]);

int applyXOR2(unsigned char XOR[9][192][16][16], int t1, int t2,int r, int i);
void WBC_AES_MixBij(unsigned char block[16], unsigned char XOR[9][192][16][16], unsigned int mixin[10][16][256], unsigned int mixout[9][16][256]);


#endif /* CIPHER_H */