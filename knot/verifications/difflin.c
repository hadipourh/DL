/*
Optimized implementation of KNOT-AEAD-128-256
State size = 256 bits
Key size = 128 bits
The original code is provided by the designers of KNOT:
http://www.sklois.ac.cn/kycg1/dbcg/202006/t20200617_565008.html

This implementation was slightly modified by  to apply the differential-linear analysis.
*/
#include "api.h"
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/random.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <inttypes.h>

typedef unsigned char u8;
typedef unsigned long long u64;
typedef long long i64;
#define RATE (64 / 8)

#define PR0_ROUNDS 52
#define PR_ROUNDS 28
#define PRF_ROUNDS 32

#define ROTR(x,n) (((x)>>(n))|((x)<<(64-(n))))
#define LOTR64(x,n) (((x)<<(n))|((x)>>(64-(n))))

#define EXT_BYTE(x,n) ((u8)((u64)(x)>>(8*(n))))
#define INS_BYTE(x,n) ((u64)(x)<<(8*(n)))
#define U64BIG(x) (x)
static const u8 constant6[63] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x21, 0x03,
		0x06, 0x0c, 0x18, 0x31, 0x22, 0x05, 0x0a, 0x14, 0x29, 0x13, 0x27, 0x0f,
		0x1e, 0x3d, 0x3a, 0x34, 0x28, 0x11, 0x23, 0x07, 0x0e, 0x1c, 0x39, 0x32,
		0x24, 0x09, 0x12, 0x25, 0x0b, 0x16, 0x2d, 0x1b, 0x37, 0x2e, 0x1d, 0x3b,
		0x36, 0x2c, 0x19, 0x33, 0x26, 0x0d, 0x1a, 0x35, 0x2a, 0x15, 0x2b, 0x17,
		0x2f, 0x1f, 0x3f, 0x3e, 0x3c, 0x38, 0x30, 0x20 };

#define ARR_SIZE(a) (sizeof((a))/sizeof((a[0])))
#define sbox(a, b, c, d, e, f, g, h)                                                                            \
{                                                                                                                             \
	t1 = ~a; t2 = b & t1;t3 = c ^ t2; h = d ^ t3; t5 = b | c; t6 = d ^ t1; g = t5 ^ t6; t8 = b ^ d; t9 = t3 & t6; e = t8 ^ t9; t11 = g & t8; f = t3 ^ t11; \
}

#define ROUND256(i) ({\
	x0^=constant6[i];\
	sbox(x0, x1, x2, x3,x4, x5, x6, x7);\
	x0=x4;\
	x1=LOTR64(x5,1);\
	x2=LOTR64(x6,8);\
	x3=LOTR64(x7,25);\
})

int crypto_aead_encrypt(unsigned char *c, unsigned long long *clen,
	const unsigned char *m, unsigned long long mlen,
	const unsigned char *ad, unsigned long long adlen,
	const unsigned char *nsec, const unsigned char *npub,
	const unsigned char *k) {

	u64 K0 = U64BIG(((u64*)k)[0]);
	u64 K1 = U64BIG(((u64*)k)[1]);
	u64 N0 = U64BIG(((u64*)npub)[0]);
	u64 N1 = U64BIG(((u64*)npub)[1]);

	u64 t1, t2, t3, t5, t6, t8, t9, t11;
	u64 x3, x2, x1, x0, x7, x6, x5, x4;
	u64 rlen, i;

	// initialization
	x0 = N0;
	x1 = N1;
	x2 = K0;
	x3 = K1;

	for (i = 0; i < PR0_ROUNDS; i++) {
		ROUND256(i);
	}
	// process associated data
	if (adlen) {
		rlen = adlen;
		while (rlen >= RATE) {
			x0 ^= U64BIG(*(u64*)ad);
			for (i = 0; i < PR_ROUNDS; i++) {
				ROUND256(i);
			}
			rlen -= RATE;
			ad += RATE;
		}
		for (i = 0; i < rlen; ++i, ++ad)
			x0 ^= INS_BYTE(*ad, i);
		x0 ^= INS_BYTE(0x01, rlen);

		for (i = 0; i < PR_ROUNDS; i++) {
			ROUND256(i);
		}
	}
	x3 ^= 0x8000000000000000;
	// process plaintext

	if (mlen) {
		rlen = mlen;
		while (rlen >= RATE) {
			x0 ^= U64BIG(*(u64*)m);
			*(u64*)c = U64BIG(x0);

			for (i = 0; i < PR_ROUNDS; i++) {
				ROUND256(i);
			}
			rlen -= RATE;
			m += RATE;
			c += RATE;
		}
		for (i = 0; i < rlen; ++i, ++m, ++c) {
			x0 ^= INS_BYTE(*m, i);
			*c = EXT_BYTE(x0, i);
		}
		x0 ^= INS_BYTE(0x01, rlen);
	}
	// finalization
	for (i = 0; i < PRF_ROUNDS; i++) {
		ROUND256(i);
	}
	// return tag

	for (i = 0; i < 8; ++i) {
		*c = EXT_BYTE(U64BIG(x0), i);
		c++;
	}
	for (i = 0; i < 8; ++i) {
		*c = EXT_BYTE(U64BIG(x1), i);
		c++;
	}
	*clen = mlen + CRYPTO_KEYBYTES;
	return 0;
}

int crypto_aead_decrypt(unsigned char *m, unsigned long long *mlen,
	unsigned char *nsec, const unsigned char *c, unsigned long long clen,
	const unsigned char *ad, unsigned long long adlen,
	const unsigned char *npub, const unsigned char *k) {

	*mlen = 0;
	if (clen < CRYPTO_KEYBYTES)
		return -1;
	u64 K0 = U64BIG(((u64*)k)[0]);
	u64 K1 = U64BIG(((u64*)k)[1]);
	u64 N0 = U64BIG(((u64*)npub)[0]);
	u64 N1 = U64BIG(((u64*)npub)[1]);

	u64 t1, t2, t3, t5, t6, t8, t9, t11;
	u64 x3, x2, x1, x0, x7, x6, x5, x4;
	u64 rlen, i;

	// initialization
	x0 = N0;
	x1 = N1;
	x2 = K0;
	x3 = K1;

	for (i = 0; i < PR0_ROUNDS; i++) {
		ROUND256(i);
	}
	// process associated data
	if (adlen) {
		rlen = adlen;
		while (rlen >= RATE) {
			x0 ^= U64BIG(*(u64*)ad);
			for (i = 0; i < PR_ROUNDS; i++) {
				ROUND256(i);
			}
			rlen -= RATE;
			ad += RATE;
		}
		for (i = 0; i < rlen; ++i, ++ad)
			x0 ^= INS_BYTE(*ad, i);
		x0 ^= INS_BYTE(0x01, rlen);

		for (i = 0; i < PR_ROUNDS; i++) {
			ROUND256(i);
		}
	}
	x3 ^= 0x8000000000000000;
	// process plaintext

	rlen = clen - CRYPTO_KEYBYTES;

	if (rlen) {
		while (rlen >= RATE) {
			*(u64*)m = U64BIG(x0) ^ *(u64*)c;
			x0 = U64BIG(*((u64*)c));

			for (i = 0; i < PR_ROUNDS; i++) {
				ROUND256(i);
			}
			rlen -= RATE;
			m += RATE;
			c += RATE;
		}
		for (i = 0; i < rlen; ++i, ++m, ++c) {
			*m = EXT_BYTE(x0, i) ^ *c;
			x0 &= ~INS_BYTE(0xff, i);
			x0 |= INS_BYTE(*c, i);
		}
		x0 ^= INS_BYTE(0x01, rlen);
	}
	// finalization
	for (i = 0; i < PRF_ROUNDS; i++) {
		ROUND256(i);
	}
	// return -1 if verification fails
	t1 = *((u64*)c);
	c += RATE;
	t2 = *((u64*)c);

	if (t1 != U64BIG(x0) || t2 != U64BIG(x1)) {
		return -1;
	}
	*mlen = clen - CRYPTO_KEYBYTES;
	return 0;
}


int check_implementation() {
	printf("Checking the implementation...\n");
    // Define the key, nonce, plaintext, and associated data as uint8_t arrays
    uint8_t key[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                       0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F};
    uint8_t nonce[16] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                         0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F};
    uint8_t plaintext[1] = {0x00}; // PT = 00

    // Data set 1
    uint8_t associated_data_1[4] = {0x00, 0x01, 0x02, 0x03}; // AD = 00010203
    uint8_t ciphertext_1[32]; // Increased buffer size to handle the extra byte in the ciphertext
    unsigned long long ciphertext_len_1;

    // Provided ciphertext for data set 1
    uint8_t provided_ciphertext_1[17] = {
        0x97, 0xEB, 0xC5, 0x22, 0x73, 0xB3, 0xD4, 0x0E,
        0xAF, 0x77, 0x4C, 0xE8, 0x8C, 0xCE, 0x94, 0x69, 0xFD
    };

    // Perform encryption for data set 1
    int result_encryption_1 = crypto_aead_encrypt(ciphertext_1, &ciphertext_len_1,
                                                 plaintext, sizeof(plaintext),
                                                 associated_data_1, sizeof(associated_data_1),
                                                 NULL, nonce, key);

    if (result_encryption_1 == 0) {
        printf("Data Set 1:\nEncryption successful.\nCiphertext: ");
        for (unsigned long long i = 0; i < ciphertext_len_1; i++) {
            printf("%02X", ciphertext_1[i]);
        }
        printf("\nCiphertext length: %llu\n", ciphertext_len_1);

        // Verify the generated ciphertext with the provided ciphertext for data set 1
        if (memcmp(ciphertext_1, provided_ciphertext_1, ciphertext_len_1) == 0) {
            printf("\nCiphertext matches the provided ciphertext for data set 1.\n");
        } else {
            printf("\nCiphertext does not match the provided ciphertext for data set 1.\n");
        }
    } else {
        printf("Data Set 1:\nEncryption failed.\n");
        return -1;
    }

    // Data set 1 decryption
    uint8_t decrypted_text_1[1]; // Adjust the size as needed
    unsigned long long decrypted_text_len_1;
    int result_decryption_1 = crypto_aead_decrypt(decrypted_text_1, &decrypted_text_len_1,
                                                 NULL, ciphertext_1, ciphertext_len_1,
                                                 associated_data_1, sizeof(associated_data_1),
                                                 nonce, key);

    if (result_decryption_1 == 0) {
        printf("\nDecryption successful for data set 1.\nDecrypted Text: ");
        for (unsigned long long i = 0; i < decrypted_text_len_1; i++) {
            printf("%02X", decrypted_text_1[i]);
        }
        printf("\nDecrypted Text length: %llu\n", decrypted_text_len_1);

        // Verify the decrypted text with the original plaintext
        if (memcmp(decrypted_text_1, plaintext, decrypted_text_len_1) == 0) {
            printf("\nDecrypted text matches the original plaintext for data set 1.\n");
        } else {
            printf("\nDecrypted text does not match the original plaintext for data set 1.\n");
        }
    } else {
        printf("\nDecryption failed for data set 1.\n");
        return -1;
    }
    return 0;
}


// Function to apply one round of operations on a 256-bit state
static void MYROUND256(u64* state, int i) {
    u64 t1, t2, t3, t5, t6, t8, t9, t11;
    u64 x0 = state[0];
    u64 x1 = state[1];
    u64 x2 = state[2];
    u64 x3 = state[3];
    u64 x4, x5, x6, x7;

    x0 ^= constant6[i];
    sbox(x0, x1, x2, x3, x4, x5, x6, x7);
    x0 = x4;
    x1 = LOTR64(x5, 1);
    x2 = LOTR64(x6, 8);
    x3 = LOTR64(x7, 25);

    state[0] = x0;
    state[1] = x1;
    state[2] = x2;
    state[3] = x3;
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

// Function to fill the knot256_state with random values
void generate_random_state(u64 state[4]) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 8; ++j) {
			state[i] = (state[i] << 8) | (rand() & 0xff);
		}
	}
}

// Function to print the knot256_state
void print_knot256_state(const u64* state) {
    for (int i = 0; i < 4; i++) {
        printf("Row %d: %016llx\n", i, state[i]);
    }
}

uint64_t count_set_bits(uint64_t n) {
    uint64_t count = 0;
    while (n) {
        n &= (n - 1);
        count++;
    }
    return count;
}

// Function to compute the dot product of two knot256_states
uint64_t dot_product(const u64 state1[4], const u64 state2[4]) {
    uint64_t result = 0;
    for (int i = 0; i < 4; i++) {
        // result += __builtin_popcountll(state1[i] & state2[i]);
		result += count_set_bits(state1[i] & state2[i]);
    }
    return result & 0x01;
}

// Function to check if ROUND256 and MYROUND256 produce the same output
void check_my_round_function() {
	printf("Checking the MYROUND256 function...\n");
    u64 state[4];    
	u64 t1, t2, t3, t5, t6, t8, t9, t11;
	u64 x3, x2, x1, x0, x7, x6, x5, x4;
	generate_random_state(state);
	print_knot256_state(state);
	x0 = state[0];
    x1 = state[1];
    x2 = state[2];
    x3 = state[3];
    for (int i = 0; i < PR0_ROUNDS; i++) {
        ROUND256(i);
    }

    // Apply MYROUND256 on the same state
    for (int i = 0; i < PR0_ROUNDS; i++) {
        MYROUND256(state, i);
    }

    // Compare the outputs
    if (x0 == state[0] && x1 == state[1] && x2 == state[2] && x3 == state[3]) {
        printf("ROUND256 and MYROUND256 produce the same output.\n");
    } else {
        printf("ROUND256 and MYROUND256 do not produce the same output.\n");
    }
}

int main() {
	check_implementation();
	check_my_round_function();
	init_prng(time(NULL));
    u64 input_diff[4];
    u64 output_mask[4];
    u64 state_1[4] = {0x0, 0x0, 0x0, 0x0};
    u64 state_2[4];
    uint64_t output_parity_1;
    uint64_t output_parity_2;
    //#######################################################################
    //#######################################################################
    //#######################################################################
    int nrounds = 9;
	input_diff[0] = 0x4000000000000000;
	input_diff[1] = 0x8000000000000000;
	input_diff[2] = 0x0000000040000040;
	input_diff[3] = 0x0000800000000000;



	output_mask[0] = 0x0000010000000000;
	output_mask[1] = 0x0000000000000000;
	output_mask[2] = 0x0000000000000000;
	output_mask[3] = 0x0000010000000000;

    int deg = 25; // num_of_experiments = 2^deg
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
        generate_random_state(state_1);        
        for(int i = 0; i < 4; i++)
        {
            state_2[i] = (state_1[i] ^ input_diff[i]) & 0xffffffffffffffff;
        } 	
        for (int r = 0; r < nrounds; r++)
		{
			MYROUND256(state_1, r);
			MYROUND256(state_2, r);
		}
        output_parity_1 = dot_product(output_mask, state_1);
        output_parity_2 = dot_product(output_mask, state_2);
        if (output_parity_1 == output_parity_2)
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
    printf("\nInput diff:\n");
    print_knot256_state(input_diff);
    printf("\nOutput mask:\n");
    print_knot256_state(output_mask);
    printf("\nNumber of experiments = %lu = 2^(%02d)\n", num_of_experiments, deg);
    printf("\nAbsolute correlation = %lu\n", absolute_correlation);
    log_correlation = (log((double)absolute_correlation) / log(2)) - deg;
	printf("\nCorrelation = 2^(%0.2f)\n", log_correlation);
	return 0;
}
