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

UINT64 difference(int R, int N3, int* dp, int* dc)
{
    UINT64 num = 0;
    for(int i = 0; i < 32; i++) dp[i] = dp[i] & 0xf;
    for(int i = 0; i < 32; i++) dc[i] = dc[i] & 0xf;
    int k[32];
    int p1[32],p2[32];
    int c1[32],c2[32];
    // Randomly choose the master key
    for(int i = 0; i < 32; i++) k[i] = rand() & 0xf;

	for (int t = 0; t < N3; t++){
        bool flag;
		// randomly choose p1
		for(int i = 0; i < 32; i++) p1[i] = rand() & 0xf;
        // compute p2
        for(int i = 0; i < 32; i++) p2[i] = p1[i] ^ dp[i];
        // compute c1 and c2 by encrypting p1 and p2, respectively
        enc(p1, c1, k, R);
        enc(p2, c2, k, R);
		flag = true;
        for (int i = 0; i < 32; i++)
        {
            if ((c1[i] ^ c2[i]) != dc[i])
            {
                flag = 0;
            }
        }
		if (flag) {num ++;}
	}
	return num;
}

double send_differences(int R, int N1, UINT64 N2, UINT64 N3, int *dp, int *dc)
{
    // Parallel execution
    int NUM[N1];
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
        int num = 0;
        int ID = omp_get_thread_num();
        //init_prng(ID);
        for (UINT64 j = 0; j < N2; j++)
        {
            num += difference(R, N3, dp, dc);
            if ((j & STEP) == 0){
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

    printf("2^(-%f)\n", log(sum_temp) / log(2));
    printf("####################################\n");
    return sum;
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
    unsigned int task_id;
    task_id = atoi(argv[1]);
    unsigned int initial_seed;
    initial_seed = init_prng(task_id);
    int dp[32];
    int dc[32];
    bool check;
    check = test();
    printf("%s: %s\n", "Check decryption", check ? "true" : "false");
    convert_hexstr_to_statearray(DP_STR, dp);
    convert_hexstr_to_statearray(DC_STR, dc);
    //########################## Number of queries #########################
    int N1 = Nthreads;             // Number of parallel threads :  N1
    UINT64 N2 = (UINT64)1 << DEG1; // Number of bunches per thread: N2 = 2^(DEG1)
    UINT64 N3 = (UINT64)1 << DEG2; // Number of queries per bunch:  N3 = 2^(DEG2)
    //################### Number of total queries : N1*N2*N3 ###############
    double sum = 0;
    for (int i = 0; i < NUMBER_OF_EXPERIMENTS; i++)
    {
        sum += send_differences(NUMBER_OF_ROUNDS, N1, N2, N3, dp, dc);
    }
    double temp = log(NUMBER_OF_EXPERIMENTS) + log(N1) + log(N2) + log(N3);
    char name[30];
    sprintf(name, "result_%d_%d.txt", NUMBER_OF_ROUNDS, task_id);
    fic = fopen(name, "w");
    fprintf(fic, "Initial seed 0x%08X\n", initial_seed);
    fprintf(fic, "Boomerang distinguisher for %d rounds of WARP\n", NUMBER_OF_ROUNDS);
    fprintf(fic, "Input difference: \t %s\n", DP_STR);
    fprintf(fic, "Output difference: \t %s\n", DC_STR);
    double avg_pr = (temp - log(sum))/log(2);
    fprintf(fic, "Average probability = 2^(-%0.4f)\n", avg_pr);
    fprintf(fic, "Number of boomerangs thrown = 2^%d\n", (int)(temp/log(2)));
    fprintf(fic, "Number of boomerangs returned = %d\n", (int)sum);
    fclose(fic);
    printf("\nAverage probability = 2^(-%0.4f)\n", avg_pr);
    return 0;
}
