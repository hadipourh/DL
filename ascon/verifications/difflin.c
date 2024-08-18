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

#include "permutations.h"
#include <sys/random.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>

//----------------------------------
// Initialize random generator (PRNG)
//----------------------------------

unsigned int init_prng(unsigned int offset) {
    unsigned int initial_seed = 0;
    ssize_t temp;
    temp = getrandom(&initial_seed, sizeof(initial_seed), 0);
    if (temp == -1) perror("error!");
    initial_seed += offset;
    srand(initial_seed);
	printf("[+] PRNG initialized to 0x%08X\n", initial_seed);
    return initial_seed;
}

//----------------------------------
// Fill Ascon state with random bytes
//----------------------------------

forceinline void fill_ascon_state_with_random(ascon_state_t* s) {
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 8; j++) {
      s->b[i][j] = (uint8_t)rand();
    }
  }
}

//----------------------------------
// Dot product ot two Ascon states
//----------------------------------

uint64_t dot_product(const ascon_state_t* state1, const ascon_state_t* state2) {
  uint64_t result = 0;

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 8; j++) {
      result += __builtin_popcount(state1->b[i][j] & state2->b[i][j]);
    }
  }

  return result & 0x01;
}

//----------------------------------
// Print Ascon state in 5 rows
//----------------------------------

void print_ascon_state(const ascon_state_t* s) {
  for (int i = 0; i < 5; i++) {
    printf("Row %d: ", i);
    for (int j = sizeof(uint64_t) - 1; j >= 0; j--) {
      printf("%02" PRIx8, s->b[i][j]);
    }
    printf("\n");
  }
}

int main(int argc, char** argv)
{
    init_prng(time(NULL));    
    ascon_state_t input_diff;
    ascon_state_t output_mask;
    ascon_state_t state_1;
    ascon_state_t state_2;
    uint64_t output_parity_1;
    uint64_t output_parity_2;
    //#######################################################################
    //#######################################################################
    //#######################################################################
    int nrounds = 5;
    input_diff.x[0] = 0x0000000000000080;
    input_diff.x[1] = 0x0000000000000000;
    input_diff.x[2] = 0x0000000000000000;
    input_diff.x[3] = 0x0000000000000080;
    input_diff.x[4] = 0x0000000000000080;

    output_mask.x[0] = 0x6da496ddb4932449;
    output_mask.x[1] = 0x7110f752d23e65d3;
    output_mask.x[2] = 0x0000000000000000;
    output_mask.x[3] = 0x0000000000000000;
    output_mask.x[4] = 0xe631e6e25c7f614b;

    int deg = 22; // num_of_experiments = 2^deg
    //#######################################################################
    //#######################################################################
    //#######################################################################
    uint64_t num_of_experiments = 1ULL << deg;
    uint64_t counter_0 = 0;
    uint64_t counter_1 = 0;
    uint64_t absolute_correlation = 0;
    double log_correlation = 0;
    clock_t start_time, end_time;
    double execution_time;
    start_time = clock();
    for(uint64_t experiment_number = 0; experiment_number < num_of_experiments; experiment_number++)
    {
        fill_ascon_state_with_random(&state_1);        
        for(int i = 0; i < 5; i++)
        {
            state_2.x[i] = state_1.x[i] ^ input_diff.x[i];
        } 
        ascon_permutation(&state_1, nrounds);
        ascon_permutation(&state_2, nrounds);
        output_parity_1 = dot_product(&output_mask, &state_1);
        output_parity_2 = dot_product(&output_mask, &state_2);
        if (output_parity_1 == output_parity_2)
        {
            counter_0++;
        }
        else
        {
            counter_1++;
        }
    }
    // print counter_0 - counter_1
    printf("\n");
    printf("counter_0 - counter_1 = %ld\n", counter_0 - counter_1);
    if (counter_0 > counter_1)
			absolute_correlation = counter_0 - counter_1;
		else
			absolute_correlation = counter_1 - counter_0;
    end_time = clock();
    execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %.2f seconds\n", execution_time);    
    printf("\nInput diff:\n");
    print_ascon_state(&input_diff);
    printf("\nOutput mask:\n");
    print_ascon_state(&output_mask);    
    printf("\nNumber of experiments = %lu\n", num_of_experiments);
    printf("\nAbsolute correlation = %lu\n", absolute_correlation);
    log_correlation = (log(absolute_correlation) / log(2)) - deg;
		printf("\nCorrelation = 2^(%0.2f)\n", log_correlation);
    return 0;
}
