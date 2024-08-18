#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef uint16_t u16;
typedef uint8_t u8;

// #define PRINT

// #define NBROUND 32

#define KSIZE 80

static const u8 S0[16] = {14, 9, 15, 0, 13, 4, 10, 11, 1, 2, 8, 3, 7, 6, 12, 5};
static const u8 S1[16] = {4, 11, 14, 9, 15, 13, 0, 10, 7, 12, 5, 6, 2, 8, 1, 3};
static const u8 S2[16] = {1, 14, 7, 12, 15, 13, 0, 6, 11, 5, 9, 3, 2, 4, 8, 10};
static const u8 S3[16] = {7, 6, 8, 11, 0, 15, 3, 14, 9, 10, 12, 13, 5, 2, 4, 1};
static const u8 S4[16] = {14, 5, 15, 0, 7, 2, 12, 13, 1, 8, 4, 9, 11, 10, 6, 3};
static const u8 S5[16] = {2, 13, 11, 12, 15, 14, 0, 9, 7, 10, 6, 3, 1, 8, 4, 5};
static const u8 S6[16] = {11, 9, 4, 14, 0, 15, 10, 13, 6, 12, 5, 7, 3, 8, 1, 2};
static const u8 S7[16] = {13, 10, 15, 0, 14, 4, 9, 11, 2, 1, 8, 3, 7, 5, 12, 6};
static const u8 S8[16] = {8, 7, 14, 5, 15, 13, 0, 6, 11, 12, 9, 10, 2, 4, 1, 3};
static const u8 S9[16] = {11, 5, 15, 0, 7, 2, 9, 13, 4, 8, 1, 12, 14, 10, 3, 6};

void EncryptKeySchedule(int nrounds, u8 key[10], u8 output[][4]);
void Swap(u8 block[8]);
void OneRound(u8 x[8], u8 k[4]);
void Encrypt(int nrounds, u8 x[8], u8 subkey[][4]);
void OneRound_Inv(u8 y[8], u8 k[4]);
void Decrypt(int nrounds, u8 x[8], u8 subkey[][4]);

