#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include "../include/keyexp.h"
#include "../include/cipher.h"
#include "../include/Bijection.h"
#include "../include/tools.h"

void T_box_calculation(unsigned char* key, unsigned char T_box[10][16][256]) {
    // Calculation of all the Tbox (160 in total : 16 for each round of AES)
    unsigned char* expkey = ExpandKey(key, 16);  // Only time we will use the expand key
    unsigned char temp[11][16];
    int i; int j; int r; unsigned char c;
    for(i = 0; i < 11; i++) {
        for( j = 0; j < 16; j++) {
            temp[i][j] = expkey[16*i + j]; // To make it easier whe split the key in 11
        }
    }
    for(i = 0; i < 10; i++) {
        ShiftRows(temp[i]); // applying ShiftRows for each 128 bits block of the expkey except for the last one
    }
    for(r = 0; r < 9; r++) {
        for(i = 0; i < 16; i++) {
            c = 0;
            for(j = 0; j < 256; j++) {
                T_box[r][i][j] = GetSboxValue(c^temp[r][i]); //Calulation of the 144 first Tbox
                c++;
            }
        }
    }
    for(i = 0; i < 16; i++) {
        c = 0;
        for(j = 0; j < 256; j++) {
            T_box[9][i][j] = GetSboxValue(c^temp[9][i])^temp[10][i]; // different calculation for the 16 last
            c++;
        }
    }
    free(expkey);
}

int Mult_Vect(unsigned char c, unsigned char vect[4]) {
    // This function is used to help the calculation of Tyi_tables 
    // It return the value of c multiplied by a 1*4 Vector. The result is store in a 32 bit int
    int ab = (int)FastMult(c,vect[0]);  // Calculate each byte of the result
    int cd = (int)FastMult(c,vect[1]);
    int ef = (int)FastMult(c,vect[2]);
    int gh = (int)FastMult(c,vect[3]);
    int res = 0;

    res |= ab << 24; // reconstruction
    res |= cd << 16; 
    res |= ef << 8;
    res |= gh; 
    return res;
}

void Tyi_tables(unsigned int tyi_tables[4][256]) { 
    // This function calculate the 4 Tyi_tables
    unsigned char Y[4][4] = {{2,1,1,3},{3,2,1,1},{1,3,2,1},{1,1,3,2}}; // Each vector is a column of the matrix used for MixColumns
    unsigned char c;
    for(int i = 0; i < 4; i++) {
        c = 0;
        for(int j = 0; j < 256; j++) {
            tyi_tables[i][j] = Mult_Vect(c, Y[i]); // using our previous function
            c++;
        }
    }
}

void construc_table(unsigned int TBoxesTyiTables[9][16][256], unsigned char tbox[16][256], unsigned char XOR[9][192][16][16], unsigned char* key) {
    // This function calculate all  necessaries tables for a table based AES implementation. Parameters are just here to store the results (except the key)
    int i, j,r;
    unsigned char k, l;

    //creation  of XOR tables
    for(i = 0; i < 9; i++) {
        for(r = 0; r < 192; r++){
                for(k = 0; k < 16; k++) {
                    for(l = 0; l < 16; l++) {
                        XOR[i][r][k][l] = k ^ l;
                    }
                }
            }
        }

    //creation of TBoxesTyiTables and Tboxs
    unsigned int tyi_tables[4][256];
    Tyi_tables(tyi_tables);
    unsigned char temp_tbox[10][16][256];
    T_box_calculation(key, temp_tbox);
    for(r = 0; r < 9; r++) {
        for(i = 0; i < 16; i++) {
            for(j = 0; j < 256; j++) {
                TBoxesTyiTables[r][i][j] = tyi_tables[i%4][temp_tbox[r][i][j]]; // composition of Tbox with TBoxesTyiTables
            }
        }
    }

    //creation of the 16 tbox left
    for(i = 0; i < 16; i++) {
        for(j = 0; j < 256; j++) {
            tbox[i][j] =  temp_tbox[9][i][j];
        }
    }
}


unsigned char Matrix8Mult(unsigned char mat[8][8], unsigned char X) { // calculate the product mat*X where X is seen as an 8 bits vector. 
    unsigned char res = 0;
    unsigned char temp;
    for(int i = 0; i < 8; i++) {
        temp = 0;
        for(int j = 0; j < 8; j++) {
            temp ^= mat[i][j] & (X >> (7 - j) & 1); 
        }
        res = res | (temp << (7 - i));
    }
    return res;
}

unsigned int Matrix32Mult(unsigned int mat[32][32], unsigned int X) { // calculate the product mat*X where X is seen as a 32 bits vector.
    unsigned int res = 0;
    unsigned int temp;
    for(int i = 0; i < 32; i++) {
        temp = 0;
        for(int j = 0; j < 32; j++) {
            temp ^= mat[i][j] & (X >> (31 - j) & 1);
        }
        res = res | (temp << (31 - i));
    }
    return res;
}


void construc_mixing_bijection(unsigned int mixin[10][16][256], unsigned int mixout[9][16][256], unsigned int TBoxesTyiTables[9][16][256],  unsigned char tbox[16][256]) {
    // Calculate all the new table used for the table based AES implementation with mixing bij included
    // Result will be stored in mixin and mixout. mixin will replace TBoxesTyiTables and tbox
    int r,i,j;

    // calculation of mixin
    for(r = 0; r < 10; r++) {
        for(i = 0; i < 16; i++) {
            for(j = 0; j < 256; j++) {
                if(r) {
                    mixin[r][i][j] = (int)Matrix8Mult(InvMat8[r - 1][i],(unsigned char) j); // We start by applying the inverse of an 8*8 matrix
                }
                else {
                    mixin[r][i][j] = j; // TBoxesTyiTables of the first round are not composed with those mixing bijections
                }
            }
        }
    }
    for(r = 0; r < 9; r++) {
        for(i = 0; i < 16; i++) {
            for(j = 0; j < 256; j++) {
                mixin[r][i][j] = TBoxesTyiTables[r][i][mixin[r][i][j]]; // Applying TBoxesTyiTables to the result
            }
        }
    }
    for(i = 0; i < 16; i++) {
        for(j = 0; j < 256; j++) {
            mixin[9][i][j] = (int)tbox[i][mixin[9][i][j]];  // tbox for the last ones
        }
    }

    for(r = 0; r < 9; r++) {
        for(i = 0; i < 4; i++) {
            for(j = 0; j < 256; j++) {
                mixin[r][4*i][j] = Matrix32Mult(Mat32[r][i], mixin[r][4*i][j]); // Applying the 32*32 matrix tranformation
                mixin[r][4*i + 1][j] = Matrix32Mult(Mat32[r][i], mixin[r][4*i + 1][j]);
                mixin[r][4*i + 2][j] = Matrix32Mult(Mat32[r][i], mixin[r][4*i + 2][j]);
                mixin[r][4*i + 3][j] = Matrix32Mult(Mat32[r][i], mixin[r][4*i + 3][j]);
            }
        }
    }

    // Calculation of mixout
    for(r = 0; r < 9; r++) {
        for(i = 0; i < 4; i++) {
            for(j = 0; j < 256; j++) {
                mixout[r][4*i][j] = Matrix32Mult(InvMat32[r][i], j<<24); // Aplying the inverse of the matrix we just applyied to mixin
                mixout[r][4*i + 1][j] = Matrix32Mult(InvMat32[r][i], j<<16);
                mixout[r][4*i + 2][j] = Matrix32Mult(InvMat32[r][i], j<<8);
                mixout[r][4*i + 3][j] = Matrix32Mult(InvMat32[r][i], j);
            }
        }
    }
    
    /*
    I'm deeply sorry for the way I implemented the next step but it is the only way I managed to make it work
    */

    int a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
    for(r = 0; r < 9; r++) {
        for(j = 0; j < 256; j++) {
            a0 = mixout[r][0][j];
            a1 = mixout[r][1][j];
            a2 = mixout[r][2][j];
            a3 = mixout[r][3][j];
            a4 = mixout[r][4][j];
            a5 = mixout[r][5][j];
            a6 = mixout[r][6][j];
            a7 = mixout[r][7][j];
            a8 = mixout[r][8][j];
            a9 = mixout[r][9][j];
            a10 = mixout[r][10][j];
            a11 = mixout[r][11][j];
            a12 = mixout[r][12][j];
            a13 = mixout[r][13][j];
            a14 = mixout[r][14][j];
            a15 = mixout[r][15][j];

            // The idea is basically to apply a concatenated matrix transformation that will compensate next-step's ShiftRows

            mixout[r][0][j] = ((int)Matrix8Mult(Mat8[r][0],(unsigned char)((a0  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][13],(unsigned char)((a0  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][10],(unsigned char)((a0  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][7],(unsigned char)(a0 & 0xFF)));
            mixout[r][1][j] = ((int)Matrix8Mult(Mat8[r][0],(unsigned char)((a1  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][13],(unsigned char)((a1  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][10],(unsigned char)((a1  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][7],(unsigned char)(a1 & 0xFF)));
            mixout[r][2][j] = ((int)Matrix8Mult(Mat8[r][0],(unsigned char)((a2  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][13],(unsigned char)((a2  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][10],(unsigned char)((a2  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][7],(unsigned char)(a2 & 0xFF)));
            mixout[r][3][j] = ((int)Matrix8Mult(Mat8[r][0],(unsigned char)((a3  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][13],(unsigned char)((a3  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][10],(unsigned char)((a3  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][7],(unsigned char)(a3 & 0xFF)));

            mixout[r][4][j] = ((int)Matrix8Mult(Mat8[r][4],(unsigned char)((a4  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][1],(unsigned char)((a4  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][14],(unsigned char)((a4  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][11],(unsigned char)(a4 & 0xFF)));
            mixout[r][5][j] = ((int)Matrix8Mult(Mat8[r][4],(unsigned char)((a5  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][1],(unsigned char)((a5  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][14],(unsigned char)((a5  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][11],(unsigned char)(a5 & 0xFF)));
            mixout[r][6][j] = ((int)Matrix8Mult(Mat8[r][4],(unsigned char)((a6  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][1],(unsigned char)((a6  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][14],(unsigned char)((a6  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][11],(unsigned char)(a6 & 0xFF)));
            mixout[r][7][j] = ((int)Matrix8Mult(Mat8[r][4],(unsigned char)((a7  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][1],(unsigned char)((a7  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][14],(unsigned char)((a7  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][11],(unsigned char)(a7 & 0xFF)));

            mixout[r][8][j] = ((int)Matrix8Mult(Mat8[r][8],(unsigned char)((a8  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][5],(unsigned char)((a8  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][2],(unsigned char)((a8  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][15],(unsigned char)(a8 & 0xFF)));
            mixout[r][9][j] = ((int)Matrix8Mult(Mat8[r][8],(unsigned char)((a9  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][5],(unsigned char)((a9  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][2],(unsigned char)((a9  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][15],(unsigned char)(a9 & 0xFF)));
            mixout[r][10][j] = ((int)Matrix8Mult(Mat8[r][8],(unsigned char)((a10  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][5],(unsigned char)((a10  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][2],(unsigned char)((a10  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][15],(unsigned char)(a10 & 0xFF)));
            mixout[r][11][j] = ((int)Matrix8Mult(Mat8[r][8],(unsigned char)((a11  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][5],(unsigned char)((a11  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][2],(unsigned char)((a11  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][15],(unsigned char)(a11 & 0xFF)));

            mixout[r][12][j] = ((int)Matrix8Mult(Mat8[r][12],(unsigned char)((a12  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][9],(unsigned char)((a12  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][6],(unsigned char)((a12  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][3],(unsigned char)(a12 & 0xFF)));
            mixout[r][13][j] = ((int)Matrix8Mult(Mat8[r][12],(unsigned char)((a13  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][9],(unsigned char)((a13  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][6],(unsigned char)((a13  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][3],(unsigned char)(a13 & 0xFF)));
            mixout[r][14][j] = ((int)Matrix8Mult(Mat8[r][12],(unsigned char)((a14  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][9],(unsigned char)((a14  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][6],(unsigned char)((a14  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][3],(unsigned char)(a14 & 0xFF)));
            mixout[r][15][j] = ((int)Matrix8Mult(Mat8[r][12],(unsigned char)((a15  >> 24) & 0xFF)) << 24) | ((int)Matrix8Mult(Mat8[r][9],(unsigned char)((a15  >> 16) & 0xFF)) << 16) | ((int)Matrix8Mult(Mat8[r][6],(unsigned char)((a15  >> 8) & 0xFF)) << 8) | ((int)Matrix8Mult(Mat8[r][3],(unsigned char)(a15 & 0xFF)));
        }
    }
}
// Store the concatenation of bij1 and bij2 in bijres
void concatenate_bij(unsigned char bij1[16], unsigned char bij2[16], unsigned char bijres[256]) {
    for(int i = 0; i < 16; i++) {
        for(int j = 0; j < 16; j++) {
            bijres[16*i + j] = bij1[i] << 4 | bij2[j];
        }
    }
}

// This function apply a concatenation of 8 bijection on t. currentbij track for which bijection we must start and continue with the 7 nexts

int apply_concatenate_bij(int currentbij, int t) {
    int res = 0;
    int a = (unsigned int)Bij[currentbij][(t >> 28) & 0xF];
    int b = (unsigned int)Bij[currentbij + 1][(t >> 24) & 0xF];
    int c = (unsigned int)Bij[currentbij + 2][(t >> 20) & 0xF];
    int d = (unsigned int)Bij[currentbij + 3][(t >> 16) & 0xF];
    int e = (unsigned int)Bij[currentbij + 4][(t >> 12) & 0xF];
    int f = (unsigned int)Bij[currentbij + 5][(t >> 8) & 0xF];
    int g = (unsigned int)Bij[currentbij + 6][(t >> 4) & 0xF];
    int h = (unsigned int)Bij[currentbij + 7][(t) & 0xF];
    res = (a << 28) | (b << 24) | (c << 20) | (d << 16) | (e << 12) | (f << 8) | (g << 4) | h ; // reconstruction of the result
    return res;
}

int apply_concatenate_invbij(int currentbij, int t) {
    int res = 0;
    int a = (unsigned int)InvBij[currentbij][(t >> 28) & 0xF];
    int b = (unsigned int)InvBij[currentbij + 1][(t >> 24) & 0xF];
    int c = (unsigned int)InvBij[currentbij + 2][(t >> 20) & 0xF];
    int d = (unsigned int)InvBij[currentbij + 3][(t >> 16) & 0xF];
    int e = (unsigned int)InvBij[currentbij + 4][(t >> 12) & 0xF];
    int f = (unsigned int)InvBij[currentbij + 5][(t >> 8) & 0xF];
    int g = (unsigned int)InvBij[currentbij + 6][(t >> 4) & 0xF];
    int h = (unsigned int)InvBij[currentbij + 7][(t) & 0xF];
    res = (a << 28) | (b << 24) | (c << 20) | (d << 16) | (e << 12) | (f << 8) | (g << 4) | h ; // reconstruction of the result
    return res;
}


void encode_tables(unsigned int mixin[10][16][256], unsigned int mixout[9][16][256], unsigned char XOR[9][192][16][16]) {
    int currentbij = 0; // Track where we are in terms of bijection
    int currentinvbij = 0; // track where we are in terms of reciprocal bijection
    int r,i,j,m,k;
    unsigned char tempbij[256];
    unsigned int new_mixout[9][16][256];
    for(r = 0; r < 9; r++) {
        for(i = 0; i < 16; i++) {
            for(j = 0; j < 256; j++) {
                new_mixout[r][i][j] = mixout[r][i][j]; // Copy of mixout for later
            }
        }
    }

    for(r = 0; r < 9; r++) {
        for(i = 0; i < 16; i++) {
            for(j = 0; j < 256; j++) {
                    mixin[r][i][j] =  apply_concatenate_bij(currentbij, mixin[r][i][j]); // ouput encoding for TBoxesTyiTables
            }
            currentbij+=8;
        }
    }

    unsigned char original_xor[16][16]; // Will be used to get the original value of all XOR tables
    for(unsigned char c1 = 0; c1 < 16; c1++) {
        for(unsigned char c2 = 0; c2 < 16; c2++) {
            original_xor[c1][c2] = c1 ^ c2; 
        }
    }

    for(r = 0; r < 9; r++) {
        for(i = 0; i < 8; i++) {
            for(m = 0; m < 16; m++) {
                for(k = 0; k < 16; k++) {
                    XOR[r][8 * i][m][k] = original_xor[InvBij[currentinvbij][m]][InvBij[currentinvbij + 8][k]]; // input encoding of the first XOR tables (64 tables per round in total)
                    XOR[r][8 * i + 1][m][k] = original_xor[InvBij[currentinvbij + 1][m]][InvBij[currentinvbij + 9][k]];
                    XOR[r][8 * i + 2][m][k] = original_xor[InvBij[currentinvbij + 2][m]][InvBij[currentinvbij + 10][k]];
                    XOR[r][8 * i + 3][m][k] = original_xor[InvBij[currentinvbij + 3][m]][InvBij[currentinvbij + 11][k]];
                    XOR[r][8 * i + 4][m][k] = original_xor[InvBij[currentinvbij + 4][m]][InvBij[currentinvbij + 12][k]];
                    XOR[r][8 * i + 5][m][k] = original_xor[InvBij[currentinvbij + 5][m]][InvBij[currentinvbij + 13][k]];
                    XOR[r][8 * i + 6][m][k] = original_xor[InvBij[currentinvbij + 6][m]][InvBij[currentinvbij + 14][k]];
                    XOR[r][8 * i + 7][m][k] = original_xor[InvBij[currentinvbij + 7][m]][InvBij[currentinvbij + 15][k]];
                }
            }
            currentinvbij += 16;
        }
    }

    for(r = 0; r < 9; r++) {
        for(i = 0; i < 64; i++) {
            for(m = 0; m < 16; m++) {
                for(k = 0; k < 16; k++) {
                    XOR[r][i][m][k] = Bij[currentbij][XOR[r][i][m][k]]; // output encoding for those same XOR tables
                }
            }
            currentbij++;
        }
    }

    for(r = 0; r < 9; r++) {
        for(i = 0; i < 4; i++) {
            for(m = 0; m < 16; m++) {
                for(k = 0; k < 16; k++) {
                    XOR[r][64 + 8 * i][m][k] = original_xor[InvBij[currentinvbij][m]][InvBij[currentinvbij + 8][k]]; // input encoding the the 2nd XOR step tables (32 tables per round)
                    XOR[r][64 + 8 * i + 1][m][k] = original_xor[InvBij[currentinvbij + 1][m]][InvBij[currentinvbij + 9][k]];
                    XOR[r][64 + 8 * i + 2][m][k] = original_xor[InvBij[currentinvbij + 2][m]][InvBij[currentinvbij + 10][k]];
                    XOR[r][64 + 8 * i + 3][m][k] = original_xor[InvBij[currentinvbij + 3][m]][InvBij[currentinvbij + 11][k]];
                    XOR[r][64 + 8 * i + 4][m][k] = original_xor[InvBij[currentinvbij + 4][m]][InvBij[currentinvbij + 12][k]];
                    XOR[r][64 + 8 * i + 5][m][k] = original_xor[InvBij[currentinvbij + 5][m]][InvBij[currentinvbij + 13][k]];
                    XOR[r][64 + 8 * i + 6][m][k] = original_xor[InvBij[currentinvbij + 6][m]][InvBij[currentinvbij + 14][k]];
                    XOR[r][64 + 8 * i + 7][m][k] = original_xor[InvBij[currentinvbij + 7][m]][InvBij[currentinvbij + 15][k]];
                }
            }
            currentinvbij += 16;
        }
    }

    for(r = 0; r < 9; r++) {
        for(i = 0; i < 32; i++) {
            for(m = 0; m < 16; m++) {
                for(k = 0; k < 16; k++) {
                    XOR[r][64 + i][m][k] = Bij[currentbij][XOR[r][64 + i][m][k]]; //ouput encoding for those table
                }
            }
            currentbij++;
        }
    }
    
    
    for(r = 0; r < 9; r++) {
        for(i = 0; i < 16; i++) {
            concatenate_bij(InvBij[currentinvbij], InvBij[currentinvbij + 1], tempbij);
            for(j = 0; j < 256; j++) {
                mixout[r][i][j] = new_mixout[r][i][tempbij[j]]; //input encoding for mixing tables  
            }
            currentinvbij+=2;
        }
    }

    for(r = 0; r < 9; r++) {
        for(i = 0; i < 16; i++) {
            for(j = 0; j < 256; j++) {
                mixout[r][i][j] = apply_concatenate_bij(currentbij, mixout[r][i][j]); // output encoding for mixing tables
            }
            currentbij+=8;
        }
    }

    for(r = 0; r < 9; r++) {
        for(i = 0; i < 8; i++) {
            for(m = 0; m < 16; m++) {
                for(k = 0; k < 16; k++) {
                    XOR[r][96 + 8 * i][m][k] = original_xor[InvBij[currentinvbij][m]][InvBij[currentinvbij + 8][k]]; // input encoding for the 3rd step XOR tables (64 tables per round)
                    XOR[r][96 + 8 * i + 1][m][k] = original_xor[InvBij[currentinvbij + 1][m]][InvBij[currentinvbij + 9][k]];
                    XOR[r][96 + 8 * i + 2][m][k] = original_xor[InvBij[currentinvbij + 2][m]][InvBij[currentinvbij + 10][k]];
                    XOR[r][96 + 8 * i + 3][m][k] = original_xor[InvBij[currentinvbij + 3][m]][InvBij[currentinvbij + 11][k]];
                    XOR[r][96 + 8 * i + 4][m][k] = original_xor[InvBij[currentinvbij + 4][m]][InvBij[currentinvbij + 12][k]];
                    XOR[r][96 + 8 * i + 5][m][k] = original_xor[InvBij[currentinvbij + 5][m]][InvBij[currentinvbij + 13][k]];
                    XOR[r][96 + 8 * i + 6][m][k] = original_xor[InvBij[currentinvbij + 6][m]][InvBij[currentinvbij + 14][k]];
                    XOR[r][96 + 8 * i + 7][m][k] = original_xor[InvBij[currentinvbij + 7][m]][InvBij[currentinvbij + 15][k]];
                }
            }
            currentinvbij += 16;
        }
    }

    for(r = 0; r < 9; r++) {
        for(i = 0; i < 64; i++) {
            for(m = 0; m < 16; m++) {
                for(k = 0; k < 16; k++) {
                    XOR[r][96 + i][m][k] = Bij[currentbij][XOR[r][96 + i][m][k]]; // output encoding for those tables
                }
            }
            currentbij++;
        }
    }

    for(r = 0; r < 9; r++) {
        for(i = 0; i < 4; i++) {
            for(m = 0; m < 16; m++) {
                for(k = 0; k < 16; k++) {
                    XOR[r][160 + 8 * i][m][k] = original_xor[InvBij[currentinvbij][m]][InvBij[currentinvbij + 8][k]]; // input encoding for the last step XOR tables (32 per round)
                    XOR[r][160 + 8 * i + 1][m][k] = original_xor[InvBij[currentinvbij + 1][m]][InvBij[currentinvbij + 9][k]];
                    XOR[r][160 + 8 * i + 2][m][k] = original_xor[InvBij[currentinvbij + 2][m]][InvBij[currentinvbij + 10][k]];
                    XOR[r][160 + 8 * i + 3][m][k] = original_xor[InvBij[currentinvbij + 3][m]][InvBij[currentinvbij + 11][k]];
                    XOR[r][160 + 8 * i + 4][m][k] = original_xor[InvBij[currentinvbij + 4][m]][InvBij[currentinvbij + 12][k]];
                    XOR[r][160 + 8 * i + 5][m][k] = original_xor[InvBij[currentinvbij + 5][m]][InvBij[currentinvbij + 13][k]];
                    XOR[r][160 + 8 * i + 6][m][k] = original_xor[InvBij[currentinvbij + 6][m]][InvBij[currentinvbij + 14][k]];
                    XOR[r][160 + 8 * i + 7][m][k] = original_xor[InvBij[currentinvbij + 7][m]][InvBij[currentinvbij + 15][k]];
                }
            }
            currentinvbij += 16;
        }
    }

    int temparg[32] = {0,1,26,27,20,21,14,15,8,9,2,3,28,29,22,23,16,17,10,11,4,5,30,31,24,25,18,19,12,13,6,7}; // Used to tell which bijection must go where 
    for(r = 0; r < 9; r++) {
        for(i = 0; i < 32; i++) {
            for(m = 0; m < 16; m++) {
                for(k = 0; k < 16; k++) {
                    XOR[r][160 + i][m][k] = Bij[currentbij + temparg[i]][XOR[r][160 + i][m][k]]; // output encoding for those tables. More difficult because we must apply bjection while thinking of the next shiftrows
                }
            }
        }
        currentbij += 32;
    }

    unsigned int tempmixin1[9][16][256];
    unsigned int tempmixin2[9][16][256];
    for(r = 0; r < 9; r++) {
        for(i = 0; i < 16; i++) {
            concatenate_bij(InvBij[currentinvbij], InvBij[currentinvbij + 1], tempbij);
            for(j = 0; j < 256; j++) {
                tempmixin1[r][i][j] = tempbij[j];
                tempmixin2[r][i][j] = mixin[r + 1][i][j];
            }
            currentinvbij += 2;
        }
    }

    for(r = 0; r < 9; r++) {
        for(i = 0; i < 16; i++) {
            for(j = 0; j < 256; j++) {
                mixin[r + 1][i][j] = tempmixin2[r][i][tempmixin1[r][i][j]]; // Input encoding fo the mixout table (except for round 1)
            }
        }
    }
}


