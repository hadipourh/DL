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

Disclaimer: We acknowledge that the TWINE block cipher doesn't adhere to statistical assumptions 
in differential analysis, such as the random sub-key assumption
or Markov cipher assumption. The tool's primary function is to find some bounds
for the security of TWINE against differential and differential-linear cryptanalysis.
*/


#include <stdio.h>
#include <stdlib.h>  // rand(), srand()
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
#include <inttypes.h>
#include <time.h>    // time()
#include <stdbool.h>
#include <sys/random.h>
#include "twine.h"

#define Nthreads 1
#define STEP ((1 << 10) - 1)

typedef unsigned long long int UINT64;
void print_state(u8 *m);
int dot_product(u8 mask[], u8 data[]);
UINT64 dldistinguisher(int R, UINT64 N3, u8 *dp, u8 *dc);
double run_bunch_of_dldistinguishers(int R, int N1, UINT64 N2, UINT64 N3, u8 *dp, u8 *lc);
void convert_hexstr_to_statearray(char hex_str[], u8 dx[16]);
unsigned int init_prng(unsigned int offset);

// #######################################################################################################
// #######################################################################################################
// ############################## User must change only the following lines ##############################
const int DEG1 = 0;              // Number of bunches per thread: N2 = 2^(DEG1)
const int DEG2 = 23;             // Number of queries per bunch:  N3 = 2^(DEG2)
int NUMBER_OF_EXPERIMENTS = 10;   // Number of independent experiments
int NUMBER_OF_ROUNDS = 10;        // Number of rounds

char DP_STR[] = "0300000000000000";
char LC_STR[] = "0000000c00000000";
// #######################################################################################################
// #######################################################################################################
