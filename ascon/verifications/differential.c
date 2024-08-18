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
    ascon_state_t input_difference;
    ascon_state_t output_difference;
    ascon_state_t output_mask; // Used for truncated trails
    ascon_state_t state_1;
    ascon_state_t state_2;

    //#######################################################################
    //#######################################################################
    //#######################################################################
    int nrounds = 1;
    input_difference.x[0] = 0x8000000000000000;
    input_difference.x[1] = 0x0000000000000000;
    input_difference.x[2] = 0x0000000000000000;
    input_difference.x[3] = 0x8000000000000000;
    input_difference.x[4] = 0x8000000000000000;

    output_difference.x[0] = 0x0000000000000000;
    output_difference.x[1] = 0x0000000000000000;
    output_difference.x[2] = 0xc200000000000000;
    output_difference.x[3] = 0x0000000000000000;
    output_difference.x[4] = 0x0000000000000000;

    // output_mask.x[0] = 0x0000000000010000;
    // output_mask.x[1] = 0x0000000000010000;
    // output_mask.x[2] = 0x0000000000000000;
    // output_mask.x[3] = 0x0000000000010000;
    // output_mask.x[4] = 0x0000000000010000;

    output_mask.x[0] = 0xFFFFFFFFFFFFFFFF;
    output_mask.x[1] = 0xFFFFFFFFFFFFFFFF;
    output_mask.x[2] = 0xFFFFFFFFFFFFFFFF;
    output_mask.x[3] = 0xFFFFFFFFFFFFFFFF;
    output_mask.x[4] = 0xFFFFFFFFFFFFFFFF;

    int deg = 21; // num_of_experiments = 2^deg
    //#######################################################################
    //#######################################################################
    //#######################################################################
    uint64_t num_of_experiments = 1ULL << deg;
    uint64_t sum = 0;
    double log_probability = 0;
    clock_t start_time, end_time;
    double execution_time;
    start_time = clock();
    for(uint64_t experiment_number = 0; experiment_number < num_of_experiments; experiment_number++)
    {
        fill_ascon_state_with_random(&state_1);        
        for(int i = 0; i < 5; i++)
        {
            state_2.x[i] = state_1.x[i] ^ input_difference.x[i];
        }
        ascon_permutation(&state_1, nrounds);
        ascon_permutation(&state_2, nrounds);        
        int flag = 1;
        for(int i = 0; i < 5; i++)
        {
            if(((state_1.x[i] ^ state_2.x[i]) & output_mask.x[i]) != (output_difference.x[i] & output_mask.x[i]))
            {
                flag = 0;
                break;
            }
        }
        if (flag == 1)
        {
            sum++;
        }
    }
    end_time = clock();
    execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;    
    printf("Execution time: %.2f seconds\n", execution_time);    
    printf("\nInput difference:\n");
    print_ascon_state(&input_difference);
    printf("\nOutput difference:\n");
    print_ascon_state(&output_difference);
    log_probability = (log(sum) / log(2)) - deg;
		printf("\nProbability = 2^(%0.2f)\n", log_probability);
    return 0;
}

// //#######################################################################
// //#######################################################################
// //#######################################################################
// // Some tests:
// int nrounds = 2;
// input_difference.x[0] = 0x0000004000000000;
// input_difference.x[1] = 0x0000000000000000;
// input_difference.x[2] = 0x0000000000000000;
// input_difference.x[3] = 0x0000004000000000;
// input_difference.x[4] = 0x0000004000000000;

// output_difference.x[0] = 0x0000000000000000;
// output_difference.x[1] = 0x0000000000000000;
// output_difference.x[2] = 0x0000005004000000;
// output_difference.x[3] = 0x0000006118708000;
// output_difference.x[4] = 0x0000000000000000;
// int deg = 20; // num_of_experiments = 2^deg

