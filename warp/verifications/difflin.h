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

#define Nthreads 1
#define STEP ((1 << 9) - 1)

typedef unsigned long long int UINT64;

unsigned int init_prng(unsigned int offset);
int dot_product(const int mask[], const int data[]);
void print_state(int *m);
bool test();
UINT64 bunch_of_diff_lin_tests(int R, UINT64 N3, int* dp, int* lc);
UINT64 parallel_diff_lin_tests(int R, int N1, UINT64 N2, UINT64 N3, int *dp, int *lc);
void convert_hexstr_to_statearray(char hex_str[], int dx[32]);

// #######################################################################################################
// #######################################################################################################
// ############################## User must change only the following lines ##############################
const int DEG1 = 0;
const int DEG2 = 20;
const int NUMBER_OF_EXPERIMENTS = 10;   // Number of independent experiments
const int NUMBER_OF_ROUNDS = 11;       // Number of rounds

char DP_STR[] = "00000000000000a00000000000000000";
char DC_STR[] = "00000000000000020000000000000000";
// #######################################################################################################
// #######################################################################################################




// // #######################################################################################################
// // #######################################################################################################
// // ############################## User must change only the following lines ##############################
// const int DEG1 = 0;
// const int DEG2 = 20;
// const int NUMBER_OF_EXPERIMENTS = 3;   // Number of independent experiments
// const int NUMBER_OF_ROUNDS = 11;       // Number of rounds

// char DP_STR[] = "00000000000000020000000000000000";
// char DC_STR[] = "00000000000000020000000000000000";
// // Expected correlation: 1
// // ##################################################
