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

#include "warp.h"
#include "difflin.h"

FILE *fic;

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


int dot_product(const int mask[], const int data[]) {
    int result = 0;
    for (int i = 0; i < 32; i++) {
        result += __builtin_popcount(mask[i] & data[i]);
    }
    return result & 0x01; // Check LSB for parity (even or odd)
}

void print_state(int *m)
{
    for (int i = 0; i < BR; i++)
    {
        printf("%x ", m[i]&0xf);
    }
    printf("\n");
};


bool test()
{
    int k[32], p[32], c[32];
    int temp[32];
    int R = 10;
    for(int i = 0; i < 32; i++) k[i] = rand() & 0xf;
    for(int i = 0; i < 32; i++) p[i] = rand() & 0xf;
    enc(p, c, k, R);
    for(int i = 0; i < 32; i++) temp[i] = p[i];
    dec(p, c, k, R);
    bool flag = true;
    for(int i = 0; i < 32; i++)
    {
        if (p[i] != temp[i])
        {
            flag = false;
            break;
        }
    }
    return flag;
}

UINT64 bunch_of_diff_lin_tests(int R, UINT64 N3, int* dp, int* lc)
{
    UINT64 counter_0 = 0;
    UINT64 counter_1 = 0;
    for(int i = 0; i < 32; i++) dp[i] = dp[i] & 0xf;
    for(int i = 0; i < 32; i++) lc[i] = lc[i] & 0xf;
    int k[32];
    int p1[32],p2[32];
    int c1[32],c2[32];
    // Randomly choose the master key
    for(int i = 0; i < 32; i++) k[i] = rand() & 0xf;

	for (UINT64 t = 0; t < N3; t++){
		// randomly choose p1
		for(int i = 0; i < 32; i++) p1[i] = rand() & 0xf;
        // compute p2
        for(int i = 0; i < 32; i++) p2[i] = p1[i] ^ dp[i];
        // compute c1 and c2 by encrypting p1 and p2, respectively
        enc(p1, c1, k, R);
        enc(p2, c2, k, R);
        int output_parity_1 = 0;
        int output_parity_2 = 0;
        output_parity_1 = dot_product(lc, c1);
        output_parity_2 = dot_product(lc, c2);
        if (output_parity_1 == output_parity_2) {
            counter_0++;
        }
        else {
            counter_1++;
        }
	}
    UINT64 absolute_correlation;
    if (counter_0 > counter_1) {
        absolute_correlation = counter_0 - counter_1;
    }
    else {
        absolute_correlation = counter_1 - counter_0;
    }
    return absolute_correlation;
}

UINT64 parallel_diff_lin_tests(int R, int N1, UINT64 N2, UINT64 N3, int *dp, int *lc)
{
    // Parallel execution
    UINT64 absolute_correlation[N1];
    printf("#Rounds: %d rounds\n", R);
    printf("#Total Queries = (#Threads)*(#Bunces)*(#Queries) = %d * %llu * %llu = 2^(%0.2f)\n", N1, N2, N3, log(N1 * N2 * N3) / log(2));
    printf("#Queries per thread = (#Bunches)*(#Queries) = %llu * %llu = 2^(%0.2f)\n", N2, N3, log(N2 * N3) / log(2));
    clock_t clock_timer;
    double wall_timer;
    clock_timer = clock();
    wall_timer = omp_get_wtime();
    omp_set_num_threads(N1);
    #pragma omp parallel for
    for (int counter = 0; counter < N1; counter++)
    {
        UINT64 sac = 0;
        int ID = omp_get_thread_num();
        // init_prng(ID);
        for (UINT64 j = 0; j < N2; j++)
        {
            sac += bunch_of_diff_lin_tests(R, N3, dp, lc);
            if ((j & STEP) == 0){
                printf("PID: %d  \t Bunch Number: %llu/%llu\n", ID, j, N2);
            }
        } 
        absolute_correlation[ID] = sac;
    }
    printf("%s: %0.4f\n", "time on clock", (double)(clock() - clock_timer) / CLOCKS_PER_SEC);
    printf("%s: %0.4f\n", "time on wall", omp_get_wtime() - wall_timer);
    UINT64 total_absolute_correlation = 0;
    double sum_temp = 1;
    for (int i = 0; i < N1; i++)
    {
        total_absolute_correlation += absolute_correlation[i];
    }
    printf("Absolute correlation: %lld\n", total_absolute_correlation);
    sum_temp = (double)(N1 * N2 * N3) / (double)(total_absolute_correlation);

    printf("Correlation         : 2^(-%0.2f)\n", log(sum_temp) / log(2));
    printf("#################################################################################\n");
    return total_absolute_correlation;
}

void convert_hexstr_to_statearray(char hex_str[], int dx[32])
{
    for (int i = 0; i < 32; i++)
    {
        char hex[2];
        hex[0] = hex_str[i];
        hex[1] = '\0';
        dx[i] = (int)(strtol(hex, NULL, 16) & 0xf);
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage: %s <task_id>\n", argv[0]);
        return 1;
    }

    unsigned int task_id = atoi(argv[1]);
    unsigned int initial_seed = init_prng(task_id);

    bool check = test();

    printf("%s: %s\n", "Check decryption", check ? "true" : "false");

    int dp[32]; 
    int lc[32];
    convert_hexstr_to_statearray(DP_STR, dp);
    convert_hexstr_to_statearray(DC_STR, lc);
    // Number of queries
    int N1 = Nthreads;                    // Number of parallel threads: N1
    UINT64 N2 = (UINT64)1 << DEG1;        // Number of bunches per thread: N2 = 2^(DEG1)
    UINT64 N3 = (UINT64)1 << DEG2;        // Number of queries per bunch: N3 = 2^(DEG2)
                                          // Number of total queries: N1 * N2 * N3
    UINT64 sum = 0;

    for (int i = 0; i < NUMBER_OF_EXPERIMENTS; i++) {
        sum += parallel_diff_lin_tests(NUMBER_OF_ROUNDS, N1, N2, N3, dp, lc);
    }

    double temp = log(NUMBER_OF_EXPERIMENTS) + log(N1) + log(N2) + log(N3);
    char name[30];
    sprintf(name, "result_%d_%d.txt", NUMBER_OF_ROUNDS, task_id);
    FILE* fic = fopen(name, "w");

    if (fic == NULL) {
        printf("Error opening file for writing.\n");
        return 1;
    }

    double avg_pr = (temp - log(sum)) / log(2);

    fprintf(fic, "Initial seed 0x%08X\n", initial_seed);
    fprintf(fic, "Diff-Lin distinguisher for %d rounds of WARP\n", NUMBER_OF_ROUNDS);
    fprintf(fic, "Input difference  : \t %s\n", DP_STR);
    fprintf(fic, "Output linear mask: \t %s\n", DC_STR);
    fprintf(fic, "Average correlation  = 2^(-%0.2f)\n", avg_pr);
    fprintf(fic, "Number of pairs      = 2^%d\n", (int)(temp / log(2)));
    fprintf(fic, "Number of satisfying = %llu\n", sum);

    fclose(fic);

    printf("\nAverage correlation = 2^(-%0.2f)\n", avg_pr);

    return 0;
}
