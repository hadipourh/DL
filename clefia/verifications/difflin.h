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

Disclaimer: We acknowledge that the CLEFIA block cipher doesn't adhere to statistical assumptions 
in differential analysis, such as the random sub-key assumption
or Markov cipher assumption. The tool's primary function is to find some bounds
for the security of CLEFIA against differential and differential-linear cryptanalysis.
*/

#include <stdio.h>
#include <stdlib.h>  // rand(), srand()
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
#include <time.h>    // time()
#include <stdbool.h>
#include <sys/random.h>
#include "clefia.c"

#define Nthreads 1
#define STEP ((1 << 10) - 1)

typedef unsigned long long int UINT64;
void print_state(unsigned char *data, int bytelen);
unsigned char dot_product(unsigned char a[16], unsigned char b[16]);
UINT64 boomerang(int R, int N3, unsigned char *dp, unsigned char *dc);
double send_boomerangs(int R, int N1, UINT64 N2, UINT64 N3, unsigned char *dp, unsigned char *dc);
void convert_hexstr_to_statearray(char hex_str[], unsigned char dx[8]);
unsigned int init_prng(unsigned int offset);

// #######################################################################################################
// #######################################################################################################
// ############################## User must change only the following lines ##############################
const int DEG1 = 0;
const int DEG2 = 27;
int NUMBER_OF_EXPERIMENTS = 5;   // Number of independent experiments
int NUMBER_OF_ROUNDS = 7;   // Number of rounds
char DP_STR[] = "000000000000000000080000d77e2bfc";
char LC_STR[] = "381c8e920000000000000000f5000000";

// 4 rounds
// char DP_STR[] = "00000000000000100000000000000000";
// char LC_STR[] = "00000000000010000000000000000000";

// 4 rounds
// char DP_STR[] = "00000000000000000000000000000001";
// char LC_STR[] = "00000000000000000000000000000010";
// #######################################################################################################
// #######################################################################################################
