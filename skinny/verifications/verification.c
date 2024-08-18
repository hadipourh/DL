/* 
MIT License

Copyright (c) 2024 Hosein Hadipour 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

 * Simmulation of dldistinguisher analysis for Skinny
 * Author: Hosein Hadipour
 * Contact: hsn.hadipour@gmail.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <sys/random.h>

// using namespace std;
typedef unsigned long long int UINT64;

// #define DEBUG 1
#define STEP ((1 << 10) - 1)
#define Nthreads 1         // Number of parallel threads utilized in this program
#define NumOfExperiments 5 // Number of independent experiments

// Table that encodes the parameters of the various Skinny versions:
// (block size, key size, number of rounds)
// Skinny-64-64: 32 rounds
// Skinny-64-128: 36 rounds
// Skinny-64-192: 40 rounds
// Skinny-128-128: 40 rounds
// Skinny-128-256: 48 rounds
// Skinny-128-384: 56 rounds
int versions[6][3] = {{64, 64, 32}, {64, 128, 36}, {64, 192, 40}, {128, 128, 40}, {128, 256, 48}, {128, 384, 56}};

// Packing of data is done as follows (state[i][j] stands for row i and column j):
// 0  1  2  3
// 4  5  6  7
// 8  9 10 11
// 12 13 14 15

// 4-bit Sbox
const unsigned char sbox_4[16] = {12, 6, 9, 0, 1, 10, 2, 11, 3, 8, 5, 13, 4, 14, 7, 15};
const unsigned char sbox_4_inv[16] = {3, 4, 6, 8, 12, 10, 1, 14, 9, 2, 5, 7, 0, 11, 13, 15};

// 8-bit Sbox
const unsigned char sbox_8[256] = {0x65, 0x4c, 0x6a, 0x42, 0x4b, 0x63, 0x43, 0x6b, 0x55, 0x75, 0x5a, 0x7a, 0x53, 0x73, 0x5b, 0x7b, 0x35, 0x8c, 0x3a, 0x81, 0x89, 0x33, 0x80, 0x3b, 0x95, 0x25, 0x98, 0x2a, 0x90, 0x23, 0x99, 0x2b, 0xe5, 0xcc, 0xe8, 0xc1, 0xc9, 0xe0, 0xc0, 0xe9, 0xd5, 0xf5, 0xd8, 0xf8, 0xd0, 0xf0, 0xd9, 0xf9, 0xa5, 0x1c, 0xa8, 0x12, 0x1b, 0xa0, 0x13, 0xa9, 0x05, 0xb5, 0x0a, 0xb8, 0x03, 0xb0, 0x0b, 0xb9, 0x32, 0x88, 0x3c, 0x85, 0x8d, 0x34, 0x84, 0x3d, 0x91, 0x22, 0x9c, 0x2c, 0x94, 0x24, 0x9d, 0x2d, 0x62, 0x4a, 0x6c, 0x45, 0x4d, 0x64, 0x44, 0x6d, 0x52, 0x72, 0x5c, 0x7c, 0x54, 0x74, 0x5d, 0x7d, 0xa1, 0x1a, 0xac, 0x15, 0x1d, 0xa4, 0x14, 0xad, 0x02, 0xb1, 0x0c, 0xbc, 0x04, 0xb4, 0x0d, 0xbd, 0xe1, 0xc8, 0xec, 0xc5, 0xcd, 0xe4, 0xc4, 0xed, 0xd1, 0xf1, 0xdc, 0xfc, 0xd4, 0xf4, 0xdd, 0xfd, 0x36, 0x8e, 0x38, 0x82, 0x8b, 0x30, 0x83, 0x39, 0x96, 0x26, 0x9a, 0x28, 0x93, 0x20, 0x9b, 0x29, 0x66, 0x4e, 0x68, 0x41, 0x49, 0x60, 0x40, 0x69, 0x56, 0x76, 0x58, 0x78, 0x50, 0x70, 0x59, 0x79, 0xa6, 0x1e, 0xaa, 0x11, 0x19, 0xa3, 0x10, 0xab, 0x06, 0xb6, 0x08, 0xba, 0x00, 0xb3, 0x09, 0xbb, 0xe6, 0xce, 0xea, 0xc2, 0xcb, 0xe3, 0xc3, 0xeb, 0xd6, 0xf6, 0xda, 0xfa, 0xd3, 0xf3, 0xdb, 0xfb, 0x31, 0x8a, 0x3e, 0x86, 0x8f, 0x37, 0x87, 0x3f, 0x92, 0x21, 0x9e, 0x2e, 0x97, 0x27, 0x9f, 0x2f, 0x61, 0x48, 0x6e, 0x46, 0x4f, 0x67, 0x47, 0x6f, 0x51, 0x71, 0x5e, 0x7e, 0x57, 0x77, 0x5f, 0x7f, 0xa2, 0x18, 0xae, 0x16, 0x1f, 0xa7, 0x17, 0xaf, 0x01, 0xb2, 0x0e, 0xbe, 0x07, 0xb7, 0x0f, 0xbf, 0xe2, 0xca, 0xee, 0xc6, 0xcf, 0xe7, 0xc7, 0xef, 0xd2, 0xf2, 0xde, 0xfe, 0xd7, 0xf7, 0xdf, 0xff};
const unsigned char sbox_8_inv[256] = {0xac, 0xe8, 0x68, 0x3c, 0x6c, 0x38, 0xa8, 0xec, 0xaa, 0xae, 0x3a, 0x3e, 0x6a, 0x6e, 0xea, 0xee, 0xa6, 0xa3, 0x33, 0x36, 0x66, 0x63, 0xe3, 0xe6, 0xe1, 0xa4, 0x61, 0x34, 0x31, 0x64, 0xa1, 0xe4, 0x8d, 0xc9, 0x49, 0x1d, 0x4d, 0x19, 0x89, 0xcd, 0x8b, 0x8f, 0x1b, 0x1f, 0x4b, 0x4f, 0xcb, 0xcf, 0x85, 0xc0, 0x40, 0x15, 0x45, 0x10, 0x80, 0xc5, 0x82, 0x87, 0x12, 0x17, 0x42, 0x47, 0xc2, 0xc7, 0x96, 0x93, 0x03, 0x06, 0x56, 0x53, 0xd3, 0xd6, 0xd1, 0x94, 0x51, 0x04, 0x01, 0x54, 0x91, 0xd4, 0x9c, 0xd8, 0x58, 0x0c, 0x5c, 0x08, 0x98, 0xdc, 0x9a, 0x9e, 0x0a, 0x0e, 0x5a, 0x5e, 0xda, 0xde, 0x95, 0xd0, 0x50, 0x05, 0x55, 0x00, 0x90, 0xd5, 0x92, 0x97, 0x02, 0x07, 0x52, 0x57, 0xd2, 0xd7, 0x9d, 0xd9, 0x59, 0x0d, 0x5d, 0x09, 0x99, 0xdd, 0x9b, 0x9f, 0x0b, 0x0f, 0x5b, 0x5f, 0xdb, 0xdf, 0x16, 0x13, 0x83, 0x86, 0x46, 0x43, 0xc3, 0xc6, 0x41, 0x14, 0xc1, 0x84, 0x11, 0x44, 0x81, 0xc4, 0x1c, 0x48, 0xc8, 0x8c, 0x4c, 0x18, 0x88, 0xcc, 0x1a, 0x1e, 0x8a, 0x8e, 0x4a, 0x4e, 0xca, 0xce, 0x35, 0x60, 0xe0, 0xa5, 0x65, 0x30, 0xa0, 0xe5, 0x32, 0x37, 0xa2, 0xa7, 0x62, 0x67, 0xe2, 0xe7, 0x3d, 0x69, 0xe9, 0xad, 0x6d, 0x39, 0xa9, 0xed, 0x3b, 0x3f, 0xab, 0xaf, 0x6b, 0x6f, 0xeb, 0xef, 0x26, 0x23, 0xb3, 0xb6, 0x76, 0x73, 0xf3, 0xf6, 0x71, 0x24, 0xf1, 0xb4, 0x21, 0x74, 0xb1, 0xf4, 0x2c, 0x78, 0xf8, 0xbc, 0x7c, 0x28, 0xb8, 0xfc, 0x2a, 0x2e, 0xba, 0xbe, 0x7a, 0x7e, 0xfa, 0xfe, 0x25, 0x70, 0xf0, 0xb5, 0x75, 0x20, 0xb0, 0xf5, 0x22, 0x27, 0xb2, 0xb7, 0x72, 0x77, 0xf2, 0xf7, 0x2d, 0x79, 0xf9, 0xbd, 0x7d, 0x29, 0xb9, 0xfd, 0x2b, 0x2f, 0xbb, 0xbf, 0x7b, 0x7f, 0xfb, 0xff};

// ShiftAndSwitchRows permutation
const unsigned char P[16] = {0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12};
const unsigned char P_inv[16] = {0, 1, 2, 3, 5, 6, 7, 4, 10, 11, 8, 9, 15, 12, 13, 14};

// Tweakey permutation
const unsigned char TWEAKEY_P[16] = {9, 15, 8, 13, 10, 14, 12, 11, 0, 1, 2, 3, 4, 5, 6, 7};
const unsigned char TWEAKEY_P_inv[16] = {8, 9, 10, 11, 12, 13, 14, 15, 2, 0, 4, 7, 6, 3, 5, 1};

// round constants
const unsigned char RC[62] = {
    0x01, 0x03, 0x07, 0x0F, 0x1F, 0x3E, 0x3D, 0x3B, 0x37, 0x2F,
    0x1E, 0x3C, 0x39, 0x33, 0x27, 0x0E, 0x1D, 0x3A, 0x35, 0x2B,
    0x16, 0x2C, 0x18, 0x30, 0x21, 0x02, 0x05, 0x0B, 0x17, 0x2E,
    0x1C, 0x38, 0x31, 0x23, 0x06, 0x0D, 0x1B, 0x36, 0x2D, 0x1A,
    0x34, 0x29, 0x12, 0x24, 0x08, 0x11, 0x22, 0x04, 0x09, 0x13,
    0x26, 0x0c, 0x19, 0x32, 0x25, 0x0a, 0x15, 0x2a, 0x14, 0x28,
    0x10, 0x20};

FILE *fic;

unsigned int init_prng(unsigned int offset)
{
    // int initial_seed = 0x30051991; My birthday!
    unsigned int initial_seed = 0;
    ssize_t temp;
    temp = getrandom(&initial_seed, sizeof(initial_seed), 0);
    if (temp == -1)
        perror("error!");
    initial_seed += offset;
    srand(initial_seed);
    printf("[+] PRNG initialized to 0x%08X\n", initial_seed);
    return initial_seed;
}

void string_state(unsigned char state[16], int ver)
{
    for (int i = 0; i < (versions[ver][0] >> 3); i++)
    {
        printf("%02x", state[i]);
    }
}

void string_tweak(unsigned char state[16], int ver)
{
    for (int i = 0; i < (versions[ver][1] >> 3); i++)
    {
        printf("%02x", state[i]);
    }
}

void display_matrix(unsigned char state[4][4], int ver)
{
    int i;
    unsigned char input[16];

    if (versions[ver][0] == 64)
    {
        for (i = 0; i < 8; i++)
            input[i] = ((state[(2 * i) >> 2][(2 * i) & 0x3] & 0xF) << 4) | (state[(2 * i + 1) >> 2][(2 * i + 1) & 0x3] & 0xF);
        for (i = 0; i < 8; i++)
            fprintf(fic, "%02x", input[i]);
    }
    else if (versions[ver][0] == 128)
    {
        for (i = 0; i < 16; i++)
            input[i] = state[i >> 2][i & 0x3] & 0xFF;
        for (i = 0; i < 16; i++)
            fprintf(fic, "%02x", input[i]);
    }
}

void display_cipher_state(unsigned char state[4][4], unsigned char keyCells[3][4][4], int ver)
{
    int k;

    fprintf(fic, "S = ");
    display_matrix(state, ver);
    for (k = 0; k < (int)(versions[ver][1] / versions[ver][0]); k++)
    {
        fprintf(fic, " - TK%i = ", k + 1);
        display_matrix(keyCells[k], ver);
    }
}

// Extract and apply the subtweakey to the internal state (must be the two top rows XORed together), then update the tweakey state
void AddKey(unsigned char state[4][4], unsigned char keyCells[3][4][4], int ver)
{
    int i, j, k;
    unsigned char pos;
    unsigned char keyCells_tmp[3][4][4];

    // apply the subtweakey to the internal state
    for (i = 0; i <= 1; i++)
    {
        for (j = 0; j < 4; j++)
        {
            state[i][j] ^= keyCells[0][i][j];
            if (2 * versions[ver][0] == versions[ver][1])
                state[i][j] ^= keyCells[1][i][j];
            else if (3 * versions[ver][0] == versions[ver][1])
                state[i][j] ^= keyCells[1][i][j] ^ keyCells[2][i][j];
        }
    }

    // update the subtweakey states with the permutation
    for (k = 0; k < (int)(versions[ver][1] / versions[ver][0]); k++)
    {
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < 4; j++)
            {
                // application of the TWEAKEY permutation
                pos = TWEAKEY_P[j + 4 * i];
                keyCells_tmp[k][i][j] = keyCells[k][pos >> 2][pos & 0x3];
            }
        }
    }

    // update the subtweakey states with the LFSRs
    for (k = 0; k < (int)(versions[ver][1] / versions[ver][0]); k++)
    {
        for (i = 0; i <= 1; i++)
        {
            for (j = 0; j < 4; j++)
            {
                // application of LFSRs for TK updates
                if (k == 1)
                {
                    if (versions[ver][0] == 64)
                        keyCells_tmp[k][i][j] = ((keyCells_tmp[k][i][j] << 1) & 0xE) ^ ((keyCells_tmp[k][i][j] >> 3) & 0x1) ^ ((keyCells_tmp[k][i][j] >> 2) & 0x1);
                    else
                        keyCells_tmp[k][i][j] = ((keyCells_tmp[k][i][j] << 1) & 0xFE) ^ ((keyCells_tmp[k][i][j] >> 7) & 0x01) ^ ((keyCells_tmp[k][i][j] >> 5) & 0x01);
                }
                else if (k == 2)
                {
                    if (versions[ver][0] == 64)
                        keyCells_tmp[k][i][j] = ((keyCells_tmp[k][i][j] >> 1) & 0x7) ^ ((keyCells_tmp[k][i][j]) & 0x8) ^ ((keyCells_tmp[k][i][j] << 3) & 0x8);
                    else
                        keyCells_tmp[k][i][j] = ((keyCells_tmp[k][i][j] >> 1) & 0x7F) ^ ((keyCells_tmp[k][i][j] << 7) & 0x80) ^ ((keyCells_tmp[k][i][j] << 1) & 0x80);
                }
            }
        }
    }

    for (k = 0; k < (int)(versions[ver][1] / versions[ver][0]); k++)
    {
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < 4; j++)
            {
                keyCells[k][i][j] = keyCells_tmp[k][i][j];
            }
        }
    }
}

// Extract and apply the subtweakey to the internal state (must be the two top rows XORed together), then update the tweakey state (inverse function}
void AddKey_inv(unsigned char state[4][4], unsigned char keyCells[3][4][4], int ver)
{
    int i, j, k;
    unsigned char pos;
    unsigned char keyCells_tmp[3][4][4];

    // update the subtweakey states with the permutation
    for (k = 0; k < (int)(versions[ver][1] / versions[ver][0]); k++)
    {
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < 4; j++)
            {
                // application of the inverse TWEAKEY permutation
                pos = TWEAKEY_P_inv[j + 4 * i];
                keyCells_tmp[k][i][j] = keyCells[k][pos >> 2][pos & 0x3];
            }
        }
    }

    // update the subtweakey states with the LFSRs
    for (k = 0; k < (int)(versions[ver][1] / versions[ver][0]); k++)
    {
        for (i = 2; i <= 3; i++)
        {
            for (j = 0; j < 4; j++)
            {
                // application of inverse LFSRs for TK updates
                if (k == 1)
                {
                    if (versions[ver][0] == 64)
                        keyCells_tmp[k][i][j] = ((keyCells_tmp[k][i][j] >> 1) & 0x7) ^ ((keyCells_tmp[k][i][j] << 3) & 0x8) ^ ((keyCells_tmp[k][i][j]) & 0x8);
                    else
                        keyCells_tmp[k][i][j] = ((keyCells_tmp[k][i][j] >> 1) & 0x7F) ^ ((keyCells_tmp[k][i][j] << 7) & 0x80) ^ ((keyCells_tmp[k][i][j] << 1) & 0x80);
                }
                else if (k == 2)
                {
                    if (versions[ver][0] == 64)
                        keyCells_tmp[k][i][j] = ((keyCells_tmp[k][i][j] << 1) & 0xE) ^ ((keyCells_tmp[k][i][j] >> 3) & 0x1) ^ ((keyCells_tmp[k][i][j] >> 2) & 0x1);
                    else
                        keyCells_tmp[k][i][j] = ((keyCells_tmp[k][i][j] << 1) & 0xFE) ^ ((keyCells_tmp[k][i][j] >> 7) & 0x01) ^ ((keyCells_tmp[k][i][j] >> 5) & 0x01);
                }
            }
        }
    }

    for (k = 0; k < (int)(versions[ver][1] / versions[ver][0]); k++)
    {
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < 4; j++)
            {
                keyCells[k][i][j] = keyCells_tmp[k][i][j];
            }
        }
    }

    // apply the subtweakey to the internal state
    for (i = 0; i <= 1; i++)
    {
        for (j = 0; j < 4; j++)
        {
            state[i][j] ^= keyCells[0][i][j];
            if (2 * versions[ver][0] == versions[ver][1])
                state[i][j] ^= keyCells[1][i][j];
            else if (3 * versions[ver][0] == versions[ver][1])
                state[i][j] ^= keyCells[1][i][j] ^ keyCells[2][i][j];
        }
    }
}

// Apply the constants: using a LFSR counter on 6 bits, we XOR the 6 bits to the first 6 bits of the internal state
void AddConstants(unsigned char state[4][4], int r)
{
    state[0][0] ^= (RC[r] & 0xf);
    state[1][0] ^= ((RC[r] >> 4) & 0x3);
    state[2][0] ^= 0x2;
}

// apply the 4-bit Sbox
void SubCell4(unsigned char state[4][4])
{
    int i, j;
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            state[i][j] = sbox_4[state[i][j]];
}

// apply the 4-bit inverse Sbox
void SubCell4_inv(unsigned char state[4][4])
{
    int i, j;
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            state[i][j] = sbox_4_inv[state[i][j]];
}

// apply the 8-bit Sbox
void SubCell8(unsigned char state[4][4])
{
    int i, j;
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            state[i][j] = sbox_8[state[i][j]];
}

// apply the 8-bit inverse Sbox
void SubCell8_inv(unsigned char state[4][4])
{
    int i, j;
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            state[i][j] = sbox_8_inv[state[i][j]];
}

// Apply the ShiftRows function
void ShiftRows(unsigned char state[4][4])
{
    int i, j, pos;

    unsigned char state_tmp[4][4];
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            // application of the ShiftRows permutation
            pos = P[j + 4 * i];
            state_tmp[i][j] = state[pos >> 2][pos & 0x3];
        }
    }

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            state[i][j] = state_tmp[i][j];
        }
    }
}

// Apply the inverse ShiftRows function
void ShiftRows_inv(unsigned char state[4][4])
{
    int i, j, pos;

    unsigned char state_tmp[4][4];
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            // application of the inverse ShiftRows permutation
            pos = P_inv[j + 4 * i];
            state_tmp[i][j] = state[pos >> 2][pos & 0x3];
        }
    }

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            state[i][j] = state_tmp[i][j];
        }
    }
}

// Apply the linear diffusion matrix
// M =
// 1 0 1 1
// 1 0 0 0
// 0 1 1 0
// 1 0 1 0
void MixColumn(unsigned char state[4][4])
{
    int j;
    unsigned char temp;

    for (j = 0; j < 4; j++)
    {
        state[1][j] ^= state[2][j];
        state[2][j] ^= state[0][j];
        state[3][j] ^= state[2][j];

        temp = state[3][j];
        state[3][j] = state[2][j];
        state[2][j] = state[1][j];
        state[1][j] = state[0][j];
        state[0][j] = temp;
    }
}

// Apply the inverse linear diffusion matrix
void MixColumn_inv(unsigned char state[4][4])
{
    int j;
    unsigned char temp;

    for (j = 0; j < 4; j++)
    {
        temp = state[3][j];
        state[3][j] = state[0][j];
        state[0][j] = state[1][j];
        state[1][j] = state[2][j];
        state[2][j] = temp;

        state[3][j] ^= state[2][j];
        state[2][j] ^= state[0][j];
        state[1][j] ^= state[2][j];
    }
}

// decryption function of Skinny
void dec(unsigned char *input, const unsigned char *userkey, int ver, int r)
{
    unsigned char state[4][4];
    unsigned char dummy[4][4] = {{0}};
    unsigned char keyCells[3][4][4];
    int i;

    memset(keyCells, 0, 48);
    for (i = 0; i < 16; i++)
    {
        if (versions[ver][0] == 64)
        {
            if (i & 1)
            {
                state[i >> 2][i & 0x3] = input[i >> 1] & 0xF;
                keyCells[0][i >> 2][i & 0x3] = userkey[i >> 1] & 0xF;
                if (versions[ver][1] >= 128)
                    keyCells[1][i >> 2][i & 0x3] = userkey[(i + 16) >> 1] & 0xF;
                if (versions[ver][1] >= 192)
                    keyCells[2][i >> 2][i & 0x3] = userkey[(i + 32) >> 1] & 0xF;
            }
            else
            {
                state[i >> 2][i & 0x3] = (input[i >> 1] >> 4) & 0xF;
                keyCells[0][i >> 2][i & 0x3] = (userkey[i >> 1] >> 4) & 0xF;
                if (versions[ver][1] >= 128)
                    keyCells[1][i >> 2][i & 0x3] = (userkey[(i + 16) >> 1] >> 4) & 0xF;
                if (versions[ver][1] >= 192)
                    keyCells[2][i >> 2][i & 0x3] = (userkey[(i + 32) >> 1] >> 4) & 0xF;
            }
        }
        else if (versions[ver][0] == 128)
        {
            state[i >> 2][i & 0x3] = input[i] & 0xFF;

            keyCells[0][i >> 2][i & 0x3] = userkey[i] & 0xFF;
            if (versions[ver][1] >= 256)
                keyCells[1][i >> 2][i & 0x3] = userkey[i + 16] & 0xFF;
            if (versions[ver][1] >= 384)
                keyCells[2][i >> 2][i & 0x3] = userkey[i + 32] & 0xFF;
        }
    }

    for (i = r - 1; i >= 0; i--)
    {
        AddKey(dummy, keyCells, ver);
    }

#ifdef DEBUG
    fprintf(fic, "DEC - initial state:                     ");
    display_cipher_state(state, keyCells, ver);
    fprintf(fic, "\n");
#endif

    for (i = r - 1; i >= 0; i--)
    {
        MixColumn_inv(state);
#ifdef DEBUG
        fprintf(fic, "DEC - round %.2i - after MixColumn_inv:    ", i);
        display_cipher_state(state, keyCells, ver);
        fprintf(fic, "\n");
#endif
        ShiftRows_inv(state);
#ifdef DEBUG
        fprintf(fic, "DEC - round %.2i - after ShiftRows_inv:    ", i);
        display_cipher_state(state, keyCells, ver);
        fprintf(fic, "\n");
#endif
        AddKey_inv(state, keyCells, ver);
#ifdef DEBUG
        fprintf(fic, "DEC - round %.2i - after AddKey_inv:       ", i);
        display_cipher_state(state, keyCells, ver);
        fprintf(fic, "\n");
#endif
        AddConstants(state, i);
#ifdef DEBUG
        fprintf(fic, "DEC - round %.2i - after AddConstants_inv: ", i);
        display_cipher_state(state, keyCells, ver);
        fprintf(fic, "\n");
#endif
        if (versions[ver][0] == 64)
            SubCell4_inv(state);
        else
            SubCell8_inv(state);
#ifdef DEBUG
        fprintf(fic, "DEC - round %.2i - after SubCell_inv:      ", i);
        display_cipher_state(state, keyCells, ver);
        fprintf(fic, "\n");
#endif
    }

#ifdef DEBUG
    fprintf(fic, "DEC - final state:                       ");
    display_cipher_state(state, keyCells, ver);
    fprintf(fic, "\n");
#endif

    if (versions[ver][0] == 64)
    {
        for (i = 0; i < 8; i++)
            input[i] = ((state[(2 * i) >> 2][(2 * i) & 0x3] & 0xF) << 4) | (state[(2 * i + 1) >> 2][(2 * i + 1) & 0x3] & 0xF);
    }
    else if (versions[ver][0] == 128)
    {
        for (i = 0; i < 16; i++)
            input[i] = state[i >> 2][i & 0x3] & 0xFF;
    }
}

// encryption function of Skinny
void enc(unsigned char *input, const unsigned char *userkey, int ver, int r)
{
    unsigned char state[4][4];
    unsigned char keyCells[3][4][4];
    int i;

    memset(keyCells, 0, 48);
    for (i = 0; i < 16; i++)
    {
        if (versions[ver][0] == 64)
        {
            if (i & 1)
            {
                state[i >> 2][i & 0x3] = input[i >> 1] & 0xF;
                keyCells[0][i >> 2][i & 0x3] = userkey[i >> 1] & 0xF;
                if (versions[ver][1] >= 128)
                    keyCells[1][i >> 2][i & 0x3] = userkey[(i + 16) >> 1] & 0xF;
                if (versions[ver][1] >= 192)
                    keyCells[2][i >> 2][i & 0x3] = userkey[(i + 32) >> 1] & 0xF;
            }
            else
            {
                state[i >> 2][i & 0x3] = (input[i >> 1] >> 4) & 0xF;
                keyCells[0][i >> 2][i & 0x3] = (userkey[i >> 1] >> 4) & 0xF;
                if (versions[ver][1] >= 128)
                    keyCells[1][i >> 2][i & 0x3] = (userkey[(i + 16) >> 1] >> 4) & 0xF;
                if (versions[ver][1] >= 192)
                    keyCells[2][i >> 2][i & 0x3] = (userkey[(i + 32) >> 1] >> 4) & 0xF;
            }
        }
        else if (versions[ver][0] == 128)
        {
            state[i >> 2][i & 0x3] = input[i] & 0xFF;
            keyCells[0][i >> 2][i & 0x3] = userkey[i] & 0xFF;
            if (versions[ver][1] >= 256)
                keyCells[1][i >> 2][i & 0x3] = userkey[i + 16] & 0xFF;
            if (versions[ver][1] >= 384)
                keyCells[2][i >> 2][i & 0x3] = userkey[i + 32] & 0xFF;
        }
    }

#ifdef DEBUG
    fprintf(fic, "ENC - initial state:                 ");
    display_cipher_state(state, keyCells, ver);
    fprintf(fic, "\n");
#endif

    for (i = 0; i < r; i++)
    {
        if (versions[ver][0] == 64)
            SubCell4(state);
        else
            SubCell8(state);
#ifdef DEBUG
        fprintf(fic, "ENC - round %.2i - after SubCell:      ", i);
        display_cipher_state(state, keyCells, ver);
        fprintf(fic, "\n");
#endif
        AddConstants(state, i);
#ifdef DEBUG
        fprintf(fic, "ENC - round %.2i - after AddConstants: ", i);
        display_cipher_state(state, keyCells, ver);
        fprintf(fic, "\n");
#endif
        AddKey(state, keyCells, ver);
#ifdef DEBUG
        fprintf(fic, "ENC - round %.2i - after AddKey:       ", i);
        display_cipher_state(state, keyCells, ver);
        fprintf(fic, "\n");
#endif
        ShiftRows(state);
#ifdef DEBUG
        fprintf(fic, "ENC - round %.2i - after ShiftRows:    ", i);
        display_cipher_state(state, keyCells, ver);
        fprintf(fic, "\n");
#endif
        MixColumn(state);
#ifdef DEBUG
        fprintf(fic, "ENC - round %.2i - after MixColumn:    ", i);
        display_cipher_state(state, keyCells, ver);
        fprintf(fic, "\n");
#endif
    } // The last subtweakey should not be added

#ifdef DEBUG
    fprintf(fic, "ENC - final state:                   ");
    display_cipher_state(state, keyCells, ver);
    fprintf(fic, "\n");
#endif

    if (versions[ver][0] == 64)
    {
        for (i = 0; i < 8; i++)
            input[i] = ((state[(2 * i) >> 2][(2 * i) & 0x3] & 0xF) << 4) | (state[(2 * i + 1) >> 2][(2 * i + 1) & 0x3] & 0xF);
    }
    else if (versions[ver][0] == 128)
    {
        for (i = 0; i < 16; i++)
            input[i] = state[i >> 2][i & 0x3] & 0xFF;
    }
}

// generate test vectors for all the versions of Skinny
void TestVectors(int ver)
{
    unsigned char p[16];
    unsigned char c[16];
    unsigned char k[48];
    int n;

    for (n = 1; n < 10; n++)
    {
        int i;
        for (i = 0; i < (versions[ver][0] >> 3); i++)
            c[i] = p[i] = rand() & 0xff;
        for (i = 0; i < (versions[ver][0] >> 3); i++)
            printf("%02x", p[i]);
        printf("\n");
        for (i = 0; i < (versions[ver][1] >> 3); i++)
            k[i] = rand() & 0xff;
        fprintf(fic, "TK = ");
        for (i = 0; i < (versions[ver][1] >> 3); i++)
            fprintf(fic, "%02x", k[i]);
        fprintf(fic, "\n");
        fprintf(fic, "P =  ");
        for (i = 0; i < (versions[ver][0] >> 3); i++)
            fprintf(fic, "%02x", p[i]);
        fprintf(fic, "\n");
        enc(c, k, ver, 10);
        fprintf(fic, "C =  ");
        for (i = 0; i < (versions[ver][0] >> 3); i++)
            fprintf(fic, "%02x", c[i]);
        fprintf(fic, "\n");
        dec(c, k, ver, 10);
        fprintf(fic, "P' = ");
        for (i = 0; i < (versions[ver][0] >> 3); i++)
            fprintf(fic, "%02x", c[i]);
        fprintf(fic, "\n\n");
    }
}

int dot_product(unsigned char a[16], unsigned char b[16])
{
    int result = 0;
    for (int i = 0; i < 16; i++) {
        result += __builtin_popcount(a[i] & b[i]);
    }
    return result & 0x01; // Check LSB for parity (even or odd)
}

UINT64 dldistinguisher(int r, int ver, unsigned long long N3, unsigned char *DP, unsigned char *LC, unsigned char *DK)
{
    unsigned char p1[16], p2[16];
    unsigned char k1[48], k2[48];
    unsigned char DC[16];

    // randomly choose k1
    for (int i = 0; i < (versions[ver][1] >> 3); i++)
        k1[i] = rand() & 0xff;
    // derive k2
    for (int i = 0; i < (versions[ver][1] >> 3); i++)
        k2[i] = k1[i] ^ DK[i];
    UINT64 counter_0 = 0;
    UINT64 counter_1 = 0;
    UINT64 absolute_correlation = 0;
    for (UINT64 t = 0; t < N3; t++)
    {
        // randomly choose p1
        for (int i = 0; i < (versions[ver][0] >> 3); i++)
        {
            p1[i] = rand() & 0xff;
        }
        // derive p2
        for (int i = 0; i < (versions[ver][0] >> 3); i++)
        {
            p2[i] = p1[i] ^ DP[i];
        }
        enc(p1, k1, ver, r);
        enc(p2, k2, ver, r);
        for (int i = 0; i < (versions[ver][0] >> 3); i++)
        {
            DC[i] = (p1[i] ^ p2[i]) & 0xff;
        }
        if (dot_product(DC, LC) == 0)
        {
            counter_0++;
        }
        else
        {
            counter_1++;
        }
    };
    if (counter_0 > counter_1)
    {
        absolute_correlation = counter_0 - counter_1;
    }
    else
    {
        absolute_correlation = counter_1 - counter_0;
    }
    return absolute_correlation;
}

double test_dldistinguishers(int R, int ver, int N1, UINT64 N2, UINT64 N3, unsigned char *DP, unsigned char *LC, unsigned char *DK)
{
    // Parallel execution
    UINT64 NUM[N1];
    printf("#Rounds: %d rounds\n", R);
    printf("#Total Queries = (#Parallel threads) * (#Bunches per thread) * (#Queries per bunch) = %d * %llu * %llu = 2^(%f)\n", N1, N2, N3, log(N1 * N2 * N3) / log(2));
    printf("#Queries per thread = (#Bunches per thread) * (#Queries per bunch) = %llu * %llu = 2^(%f)\n", N2, N3, log(N2 * N3) / log(2));
    clock_t clock_timer;
    double wall_timer;
    clock_timer = clock();
    wall_timer = omp_get_wtime();
    omp_set_num_threads(N1);
    #pragma omp parallel for
    for (int counter = 0; counter < N1; counter++)
    {
        UINT64 num = 0;
        int ID = omp_get_thread_num();
        for (UINT64 j = 0; j < N2; j++)
        {
            num += dldistinguisher(R, ver, N3, DP, LC, DK);
            if ((j & STEP) == 0)
            {
                printf("PID: %d  \t Bunch Number: %llu/%llu\n", ID, j, N2);
            }
        }
        NUM[ID] = num;
    }
    printf("%s: %0.4f\n", "time on clock", (double)(clock() - clock_timer) / CLOCKS_PER_SEC);
    printf("%s: %0.4f\n", "time on wall", omp_get_wtime() - wall_timer);
    double sum = 0;
    double sum_temp = 1;
    for (int i = 0; i < N1; i++)
        sum += NUM[i];
    printf("sum = %f\n", sum);
    sum_temp = (double)(N1 * N2 * N3) / sum;

    printf("2^(-%f)\n\n", log(sum_temp) / log(2));
    printf("####################################\n");
    return sum;
}

void convert_hexstr_to_statearray(int ver, char hex_str[], unsigned char dx[16])
{
    for (int i = 0; i < (versions[ver][0] >> 3); i++)
    {
        char hex[2];
        hex[0] = hex_str[2 * i];
        hex[1] = hex_str[2 * i + 1];
        dx[i] = (unsigned char)(strtol(hex, NULL, 16) & 0xff);
    }
}

void convert_hexstr_to_tweakarray(int ver, char hex_str[], unsigned char dt[48])
{
    for (int i = 0; i < (versions[ver][1] >> 3); i++)
    {
        char hex[2];
        hex[0] = hex_str[2 * i];
        hex[1] = hex_str[2 * i + 1];
        dt[i] = (unsigned char)(strtol(hex, NULL, 16) & 0xff);
    }
}

int main(int argc, char *argv[])
{
    // srand((unsigned)time(NULL)); // Initialization, should only be called once. int r = rand();
    init_prng(atoi(argv[1]));
    // //test all versions of Skinny
    // for (i = 0; i < (sizeof(versions) / sizeof(*versions)); i++)
    // {
    //     sprintf(name, "test_vectors_%i_%i.txt", versions[i][0], versions[i][1]);
    //     fic = fopen(name, "w");
    //     fprintf(fic, "\n\nSkinny-%i/%i: \n", versions[i][0], versions[i][1]);
    //     TestVectors(i);
    //     fclose(fic);
    //     printf("Generating test vectors for Skinny-%i/%i  -  saved in file test_vectors_%i_%i.txt \n", versions[i][0], versions[i][1], versions[i][0], versions[i][1]);
    // }
    unsigned char DP[16];
    unsigned char LC[16];
    unsigned char DK[48];
    // #######################################################################################################
    // #######################################################################################################
    // ############################## User must change only the following lines ##############################
    int R = 8;   // Number of rounds
    int ver = 2; // Determine the version:
                 // [0 = Skinny-64-64]
                 // [1 = Skinny-64-128]
                 // [2 = Skinny-64-192]
                 // [3 = Skinny-128-128]
                 // [4 = Skinny-128-256]
                 // [5 = Skinny-128-384]
    char dk_str[] = "000000000000000000000000000000000000000000000000";
    char dp_str[]  = "000000000000000a";
    char dc_str[]  = "0000300030003000";
    // ########################## Number of queries #########################
    int N1 = Nthreads; // Number of paralle threads : N1
    int deg1 = 0;
    int deg2 = 20;
    UINT64 N2 = 1 << deg1; // Number of bunches per threads: N2 = 2^(deg1)
    UINT64 N3 = 1 << deg2; // Number of queries per bunches: N3 = 2^(deg2)
    // #######################################################################################################
    // #######################################################################################################

    convert_hexstr_to_statearray(ver, dp_str, DP);
    convert_hexstr_to_statearray(ver, dc_str, LC);
    convert_hexstr_to_tweakarray(ver, dk_str, DK);
    // ################### Number of total queries : N1*N2*N3 ###############
    char all_results[NumOfExperiments][20];
    double sum = 0;
    double sum_temp = 0;
    for (int i = 0; i < NumOfExperiments; i++)
    {
        init_prng(init_prng(i));
        printf("Experiment Number %d:\n", i);
        sum_temp = test_dldistinguishers(R, ver, N1, N2, N3, DP, LC, DK);
        sum += sum_temp;
        sum_temp = (double)(N1 * N2 * N3) / sum_temp;
        sprintf(all_results[i], "2^(-%0.2f), ", log(sum_temp) / log(2));
    }
    printf("A summary of all results:\n");
    for (int i = 0; i < NumOfExperiments; i++)
    {
        printf("%s", all_results[i]);
    }
    printf("\n##########################\nAverage = 2^(-%0.4f)\n",
           (log(NumOfExperiments) + log(N1) + log(N2) + log(N3) - log(sum)) / log(2));
    return 0;
}
