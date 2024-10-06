#ifndef KEYEXP_H
#define KEYEXP_H
#include <stdio.h>

void RotWord(unsigned char* word);
void SubWord(unsigned char* word);
unsigned char* ExpandKey(unsigned char* key, int key_size);

#endif /* KEYEXP_H */