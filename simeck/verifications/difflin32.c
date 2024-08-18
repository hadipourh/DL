/*
MIT License

Copyright (c) 2024 

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


//#################################################################################
#include "difflin32.h"
//#################################################################################


FILE *fic;

bool dot_product(const uint16_t *x, const uint16_t *y) {
    bool result = false;
    for (size_t i = 0; i < 2; ++i) {
        result ^= __builtin_popcount(x[i] & y[i]) % 2;
    }
    return result;
}

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

// Function to generate a random 16-bit integer
uint16_t generate_random_16bit() {
    uint16_t result = 0;
    int remainingBits = 16;

    while (remainingBits > 0) {
        // Generate a random byte
        uint16_t randomByte = rand() & 0xFF;

        // Determine how many bits to shift
        int shift = remainingBits >= 8 ? 8 : remainingBits;

        // Shift the random byte to the correct position and OR it with the result
        result <<= shift;
        result |= randomByte;

        // Update the remaining bits
        remainingBits -= shift;
    }

    return result;
}

// Function to convert binary string to hex string
void binaryToHex(const char *binary, char *hex, int hex_size) {
    uint16_t decimal = strtol(binary, NULL, 2);
    snprintf(hex, hex_size, "%0*X", hex_size / 2, decimal);
}

// Function to split binary string into left and right parts and convert to hex
void splitAndConvert(const char *binaryString, uint16_t output[]) {
    size_t binaryStringLength = strlen(binaryString);
    int leftSize, rightSize;
    
    // Determine left and right sizes
    if (binaryStringLength == 32) {
        leftSize = 16;
        rightSize = 16;
    } else if (binaryStringLength == 48) {
        leftSize = 24;
        rightSize = 24;
    } else if (binaryStringLength == 64) {
        leftSize = 32;
        rightSize = 32;
    } else {
        printf("Invalid binary string size.\n");
        return;
    }
    
    // Allocate memory for left and right strings
    char *leftBinary = (char *)malloc((leftSize + 1) * sizeof(char));
    char *rightBinary = (char *)malloc((rightSize + 1) * sizeof(char));
    
    // Split binary string into left and right parts
    strncpy(leftBinary, binaryString, leftSize);
    leftBinary[leftSize] = '\0';
    strncpy(rightBinary, binaryString + leftSize, rightSize);
    rightBinary[rightSize] = '\0';
    
    // Convert left and right parts to hex strings
    char leftHex[leftSize + 1];
    char rightHex[rightSize + 1];
    binaryToHex(leftBinary, leftHex, leftSize + 1);
    binaryToHex(rightBinary, rightHex, rightSize + 1);
    
    // Convert hex strings to uint16_t and store in output array
    output[1] = strtol(leftHex, NULL, 16);
    output[0] = strtol(rightHex, NULL, 16);
    
    // Free allocated memory
    free(leftBinary);
    free(rightBinary);
}

uint64_t dldistinguisher(int R, uint64_t N3, uint16_t *dp, uint16_t *lc)
{ 
    // Randomly choose the master key
    const uint16_t key[] = {generate_random_16bit(), generate_random_16bit(), generate_random_16bit(), generate_random_16bit()};
    uint16_t p1[] = {0x0, 0x0}, p2[] = {0x0, 0x0};
    uint64_t counter_0 = 0;
    uint64_t counter_1 = 0;
	for (uint64_t t = 0; t < N3; t++){
        uint16_t c1[] = {0x0, 0x0}, c2[] = {0x0, 0x0};
        uint16_t dc[] = {0x0, 0x0};
		// randomly choose p1
        p1[0] = generate_random_16bit();
        p1[1] = generate_random_16bit();        
        // compute p2
        for (int i = 0; i < 2; ++i) {
            p2[i] = (p1[i] ^ dp[i]) & 0xffff;
        }
        // compute c1 and c2 by encrypting p1 and p2, respectively
        simeck_32_64(R, key, p1, c1);
        simeck_32_64(R, key, p2, c2);
        // compute parity of c1 and c2
        for (int i = 0; i < 2; ++i) {
            dc[i] = (c1[i] ^ c2[i]) & 0xffff;
        }
        if (dot_product(lc, dc))
        {
            counter_1++;
        } else {
            counter_0++;
        }
	}
    uint64_t absolute_correlation;
    if (counter_0 > counter_1)
    {
        absolute_correlation = counter_0 - counter_1;
    } else {
        absolute_correlation = counter_1 - counter_0;
    }
	return absolute_correlation;
}

double run_bunch_of_dldistinguishers(int R, int N1, uint64_t N2, uint64_t N3, uint16_t *dp, uint16_t *lc)
{
    // Parallel execution
    int NUM[N1];
    printf("#Rounds: %d rounds\n", R);
    printf("#Total Queries = (#Parallel threads) * (#Bunches per thread) * (#Queries per bunch) = %d * %lu * %lu = 2^(%f)\n", N1, N2, N3, log(N1 * N2 * N3) / log(2));
    printf("#Queries per thread = (#Bunches per thread) * (#Queries per bunch) = %lu * %lu = 2^(%f)\n", N2, N3, log(N2 * N3) / log(2));
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
        for (uint64_t j = 0; j < N2; j++)
        {
            num += dldistinguisher(R, N3, dp, lc);
            if ((j & STEP) == 0){
                printf("PID: %d  \t Bunch Number: %lu/%lu\n", ID, j, N2);
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

int main(int argc, char **argv)
{
    unsigned int task_id;
    task_id = atoi(argv[1]);
    unsigned int initial_seed;
    initial_seed = init_prng(task_id);
    uint16_t dp[2];
    uint16_t lc[2];
    //########################## Number of queries #########################
    int N1 = NUMBER_OF_THREADS;                 // Number of parallel threads :  N1
    uint64_t N2 = (uint64_t)1 << DEG1; // Number of bunches per thread: N2 = 2^(deg1)
    uint64_t N3 = (uint64_t)1 << DEG2; // Number of queries per bunch:  N3 = 2^(deg2)    
    splitAndConvert(DP_STR, dp);
    splitAndConvert(LC_STR, lc);
    printf("DP: %x, %x\n", dp[1], dp[0]);
    printf("LC: %x, %x\n", lc[1], lc[0]);
    //################### Number of total queries : N1*N2*N3 ###############
    double sum = 0;
    for (int i = 0; i < NUMBER_OF_EXPERIMENTS; i++)
    {
        sum += run_bunch_of_dldistinguishers(NUMBER_OF_ROUNDS, N1, N2, N3, dp, lc);
    }
    double temp = log(NUMBER_OF_EXPERIMENTS) + log(N1) + log(N2) + log(N3);
    char name[30];
    sprintf(name, "result_%d_%d.txt", NUMBER_OF_ROUNDS, task_id);
    fic = fopen(name, "w");
    fprintf(fic, "Initial seed 0x%08X\n", initial_seed);
    fprintf(fic, "DL distinguisher for %d rounds of TWINE\n", NUMBER_OF_ROUNDS);
    // fprintf(fic, "Input difference: \t %s\n", DP_STR);
    // fprintf(fic, "Output difference: \t %s\n", LC_STR);
    double avg_pr = (temp - log(sum))/log(2);
    fprintf(fic, "Average probability = 2^(-%0.4f)\n", avg_pr);
    fprintf(fic, "Number of experiments thrown = 2^%d\n", (int)(temp/log(2)));
    fprintf(fic, "Number of successes returned = %d\n", (int)sum);
    fclose(fic);
    printf("\nAverage probability = 2^(-%0.4f)\n", avg_pr);
    return 0;
}
