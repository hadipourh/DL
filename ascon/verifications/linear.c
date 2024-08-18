/*
Implementor: 
Date: Jul 6, 2023
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

int dot_product(const ascon_state_t* state1, const ascon_state_t* state2) {
  int result = 0;

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 8; j++) {
      result ^= __builtin_popcount(state1->b[i][j] & state2->b[i][j]);
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
    ascon_state_t input_mask;
    ascon_state_t output_mask;
    ascon_state_t state;
    int input_parity;
    int output_parity;
    //#######################################################################
    //#######################################################################
    //#######################################################################
    int nrounds = 2;
    input_mask.x[0] = 0x0000000000000000;
    input_mask.x[1] = 0x4000000000008100;
    input_mask.x[2] = 0x4000000000008100;
    input_mask.x[3] = 0x0000000000000000;
    input_mask.x[4] = 0x0000000000000000;

    output_mask.x[0] = 0x0000000000000000;
    output_mask.x[1] = 0x0000000000000000;
    output_mask.x[2] = 0x7f04314f4725bb35;
    output_mask.x[3] = 0xa908e54eef7984b5;
    output_mask.x[4] = 0x0000000000000000;
    int deg = 20; // num_of_experiments = 2^deg
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
        fill_ascon_state_with_random(&state);
        input_parity = dot_product(&input_mask, &state);
        ascon_permutation(&state, nrounds);
        output_parity = dot_product(&output_mask, &state);
        if (input_parity == output_parity)
        {
            counter_0++;
        }
        else
        {
            counter_1++;
        }
    }
    if (counter_0 > counter_1)
			absolute_correlation = counter_0 - counter_1;
		else
			absolute_correlation = counter_1 - counter_0;
    end_time = clock();
    execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %.2f seconds\n", execution_time);    
    printf("\nInput mask:\n");
    print_ascon_state(&input_mask);
    printf("\nOutput mask:\n");
    print_ascon_state(&output_mask);    
    printf("\nNumber of experiments = %lu\n", num_of_experiments);
    printf("\nAbsolute correlation = %lu\n", absolute_correlation);
    log_correlation = (log(absolute_correlation) / log(2)) - deg;
		printf("\nProbability = 2^(%0.2f)\n", log_correlation);
    return 0;
}


// //#######################################################################
// //#######################################################################
// //#######################################################################
// int nrounds = 1;
// input_mask.x[0] = 0x0000000100000000;
// input_mask.x[1] = 0x0000000000000000;
// input_mask.x[2] = 0x0000000000000000;
// input_mask.x[3] = 0x0000000000000000;
// input_mask.x[4] = 0x0000000100000000;

// output_mask.x[0] = 0x0000000000000000;
// output_mask.x[1] = 0x7ba9691f32e9b888;
// output_mask.x[2] = 0x0000000000000000;
// output_mask.x[3] = 0x0000000000000000;
// output_mask.x[4] = 0x0000000000000000;
// int deg = 20; // num_of_experiments = 2^deg
//// exected correlation is 2^-1


// //#######################################################################
// //#######################################################################
// //#######################################################################
// int nrounds = 2;
// input_mask.x[0] = 0x0000000000000000;
// input_mask.x[1] = 0x4000000000008100;
// input_mask.x[2] = 0x4000000000008100;
// input_mask.x[3] = 0x0000000000000000;
// input_mask.x[4] = 0x0000000000000000;

// output_mask.x[0] = 0x0000000000000000;
// output_mask.x[1] = 0x0000000000000000;
// output_mask.x[2] = 0x7f04314f4725bb35;
// output_mask.x[3] = 0xa908e54eef7984b5;
// output_mask.x[4] = 0x0000000000000000;
// int deg = 20; // num_of_experiments = 2^deg
// // the expected correlation is 2^-4