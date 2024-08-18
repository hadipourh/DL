#include "aesni.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/random.h>
#include <stdint.h>
#include <inttypes.h>

#define NUM_OF_ENCRYPTIONS_IN_TIMING (1ULL << 22)

void init_prng(unsigned int offset);
void print_state(uint8_t *state);
double speed();
unsigned char dot_product(unsigned char a[16], unsigned char b[16]);
void convert_hexstr_to_statearray(char hex_str[], unsigned char dx[16]);
uint64_t dldistinguisher(uint8_t* master_key, uint8_t* input_difference, uint8_t* output_mask, int round_count, uint64_t N2);
int main();

void init_prng(unsigned int offset) {
    unsigned int initial_seed = 0;
    ssize_t temp;
    temp = getrandom(&initial_seed, sizeof(initial_seed), 0);
    if (temp == -1) perror("error!");
    initial_seed += offset;
    srand(initial_seed);
    printf("[+] PRNG initialized to 0x%08X\n", initial_seed);    
}

void print_state(uint8_t *state) {
    int i;
    for (i = 0; i < 16; i++) {
        printf("%02X", state[i]);
    }
    printf("\n");
}

double speed(){
    // generate a random master key and plaintext
    uint8_t master_key[16];
    uint8_t plaintext[16];
    uint8_t ciphertext[16];
    for (int i = 0; i < 16; i++){
        master_key[i] = rand() & 0xff;
        plaintext[i] = rand() & 0xff;       
    }
    // derive the subkeys
    __m128i key_schedule[20];
    aes128_load_key(master_key, key_schedule);
    clock_t start, end;
    double cpu_time_used;
    double rate;
    unsigned long long int i = 0;
    start = clock();
    for (i = 0; i < NUM_OF_ENCRYPTIONS_IN_TIMING; i++){        
        aes128_enc(key_schedule, plaintext, ciphertext);        
    }    
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    rate = (NUM_OF_ENCRYPTIONS_IN_TIMING * 16 ) / (cpu_time_used*1000000000);
    return rate;
}

unsigned char dot_product(unsigned char a[16], unsigned char b[16])
{
    int result = 0;
    for (int i = 0; i < 16; i++) {
        result += __builtin_popcount(a[i] & b[i]);
    }
    return result & 0x01; // Check LSB for parity (even or odd)
}

void convert_hexstr_to_statearray(char hex_str[], unsigned char dx[16])
{
    for (int i = 0; i < 16; i++)
    {
        char hex[3];
        hex[0] = hex_str[2*i];
        hex[1] = hex_str[2*i + 1];
        hex[2] = '\0';
        dx[i] = (unsigned char)(strtol(hex, NULL, 16) & 0xff);
    }
}

uint64_t dldistinguisher(uint8_t* master_key, uint8_t* input_difference, uint8_t* output_mask, int round_count, uint64_t N2)
{
    __m128i* key_schedule = malloc(20*sizeof(__m128));
    uint8_t* p1 = malloc(16*sizeof(uint8_t));
    uint8_t* p2 = malloc(16*sizeof(uint8_t));
    uint8_t* c1 = malloc(16*sizeof(uint8_t));
    uint8_t* c2 = malloc(16*sizeof(uint8_t));
    aes128_load_key_enc_only(master_key, key_schedule);
    uint64_t counter_0 = 0;
    uint64_t counter_1 = 0;    
    for (uint64_t i = 0; i < N2; i++){
        for (int j = 0; j < 16; j++){
            p1[j] = rand() & 0xff;
            p2[j] = p1[j] ^ input_difference[j];
        }     
        aes_encrypt_block(key_schedule, p1, c1, round_count);
        aes_encrypt_block(key_schedule, p2, c2, round_count);
        if (dot_product(c1, output_mask) == dot_product(c2, output_mask)){
            counter_0++;
        } else {
            counter_1++;
        }
    }
    free(key_schedule);
    free(p1);
    free(p2);
    free(c1);
    free(c2);
    uint64_t absolute_correlation;
    absolute_correlation = counter_0 - counter_1;
    if (counter_0 > counter_1){
        absolute_correlation = counter_0 - counter_1;
    } else {
        absolute_correlation = counter_1 - counter_0;
    }
    int64_t difference = (int64_t)(counter_0) - (int64_t)(counter_1);        
    printf("Difference = %ld\n", difference);
    return absolute_correlation;
}

int main(int argc, char *argv[])
{
    // Initialize the PRNG
    unsigned int task_id;
    task_id = atoi(argv[1]);    
    init_prng(task_id);
    // Test the AES implementation
    int status;
    status = aes128_self_test();
    if (status != 0){
        printf("AES does not work correctly!\n");
        return 0;
    } else {
        printf("AES works correctly!\n");
    }
    double cpu_time;
    cpu_time = speed();
    printf("average speed over %llu times of encryption\t: %0.02f (Gigabytes/Second)\n", NUM_OF_ENCRYPTIONS_IN_TIMING, cpu_time);
    
    // Check the distinguisher
    //##########################################################################################################################
    //##########################################################################################################################
    //##########################################################################################################################
    int DEG1 = 0;
    int DEG2 = 25;
    uint64_t N1 = 1ULL << DEG1;
    uint64_t N2 = 1ULL << DEG2;
    int NUMBER_OF_EXPERIMENTS = 10;   // Number of independent experiments
    int NUMBER_OF_ROUNDS = 3;   // Number of rounds    
    char DP_STR[] = "0000000000000000000000b400000000";
    char LC_STR[] = "0000000032ab66980000000000000000";
    //##########################################################################################################################
    
    uint8_t* input_difference = malloc(16*sizeof(uint8_t));
    uint8_t* output_mask = malloc(16*sizeof(uint8_t));
    convert_hexstr_to_statearray(DP_STR, input_difference);
    convert_hexstr_to_statearray(LC_STR, output_mask);

   
    double sum = 0;
    for(int n = 0; n < NUMBER_OF_EXPERIMENTS; n++)
    {
        double num = 0;
        clock_t clock_timer;        
        clock_timer = clock();
        for(uint64_t i = 0; i < N1; i++){
            uint8_t* master_key = malloc(16*sizeof(uint8_t));
            for (int j = 0; j < 16; j++){
                master_key[j] = rand() & 0xff;
            }
            num += dldistinguisher(master_key, input_difference, output_mask, NUMBER_OF_ROUNDS, N2);
            free(master_key);
        }
        double elapsed_time = (double)(clock() - clock_timer) / CLOCKS_PER_SEC;
        printf("Execution time: %0.2f\n", elapsed_time);
        sum += num;
        double temp = log(N1) + log(N2);
        double avg_cr = (temp - log(num))/log(2);
        printf("\nCorrelation = 2^(-%0.4f)\n", avg_cr);
        printf("####################################\n");
    }
    double temp = log(NUMBER_OF_EXPERIMENTS) + log(N1) + log(N2);
    double avg_cr = (temp - log(sum))/log(2);    
    printf("\nAverage correlation = 2^(-%0.4f)\n", avg_cr);
    printf("####################################\n");
    return 0;
}
