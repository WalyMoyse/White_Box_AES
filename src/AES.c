#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <err.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "../include/AES.h"
#include "../include/cipher.h"
#include "../include/keyexp.h"
#include "../include/tools.h"
#include "../include/Bijection.h"
#include <limits.h>


int main() { 


    //Test without mixing bijection
    unsigned char Key[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F};
    unsigned char Text[16] = {0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff};
    unsigned char XOR[9][192][16][16];
    unsigned char tbox[16][256];
    unsigned int TBoxesTyiTables[9][16][256];
    construc_table(TBoxesTyiTables,tbox,XOR,Key);
    WBC_AES_No_MixBij(Text,tbox,TBoxesTyiTables,XOR);
    printf("Expected result according to NIST standard :\n");
    printf("69 C4 E0 D8 6A 7B 04 30 D8 CD B7 80 70 B4 C5 5A\n");
    printf("Result without mixing bijections : \n");
    for(int i = 0; i <16; i++) {
        printf("%02X ", Text[i]);
    }
    printf("\n");

    // same test but we add mixing bijections
    unsigned char Text2[16] = {0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff};
    unsigned int mixin[10][16][256];
    unsigned int mixout[9][16][256];
    construc_mixing_bijection(mixin,mixout,TBoxesTyiTables,tbox);
    WBC_AES_MixBij(Text2, XOR, mixin,mixout);
    printf("Result with mixing bijections :\n");
    for(int i = 0; i < 16; i++) {
        printf("%02X ", Text2[i]);
    }
    printf("\n");

    // Now we add internal encoding
    unsigned char Text3[16] = {0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff};
    encode_tables(mixin,mixout,XOR);
    WBC_AES_MixBij(Text3, XOR, mixin,mixout);
    printf("Result with internal encoding :\n");
    for(int i = 0; i < 16; i++) {
        printf("%02X ", Text3[i]);
    }
    printf("\n");
}
