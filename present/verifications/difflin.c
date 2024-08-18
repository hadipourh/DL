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

#include <stdio.h>    //Standard C headers...
#include <stdint.h> 
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <getopt.h>
#include <sys/random.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#include "boxes.inc"  //S-Boxes and P-Boxes
#include "verbose.inc"		

//----------------------------------
// Macros for bit manipulation
//----------------------------------                       returns...
#define high45_64(h45in) 		( (uint64_t)h45in >> 9 ) //45 msb as lsb
#define high61_64(h4in) 		( (uint64_t)h4in >> 3 )	//61 msb as lsb
#define high4_64(h4in) 			( (uint64_t)h4in >> 60 ) //4 msb as lsb
#define high8to4_64(h8in) 	( ((uint64_t)h8in >> 56)&0x0F ) //4 msb as 2. lsb
#define high16_64(h16in) 		( (uint64_t)h16in >> 48 ) //16 msb as lsb
#define high1_64(h1in) 			( (uint64_t)h1in >> 63 ) //msb as lsb
#define low4_64(l4in) 			( (uint64_t)l4in << 60 ) //4 lsb as msb
#define low8to4_64(l4in) 		( (uint64_t)l4in << 56 ) //4 lsb as 2. msb
#define low16_64(l4in) 			( (uint64_t)l4in << 48 ) //4 lsb as msb
#define rotate1l_64(r1lin)	( high1_64(r1lin) | ( r1lin << 1 ) ) //input rotated left (1x)
#define rotate1r_64(r1rin)	( high1_64(r1rin) | ( r1rin >> 1 ) ) //input rotated right (1x)
#define rotate4l_64(r4lin)	( high4_64(r4lin) | ( r4lin << 4 ) ) //input rotated left (4x)
#define rotate4r_64(r4rin)	( high4_64(r4rin) | ( r4rin >> 4 ) )	

//----------------------------------
// Function prototypes
//----------------------------------
uint64_t encrypt( uint64_t, uint64_t*, uint16_t, _Bool );
uint64_t decrypt( uint64_t, uint64_t*, uint16_t, _Bool );
uint64_t* key_schedule( uint64_t, uint64_t, uint16_t, _Bool, _Bool );
unsigned int init_prng(unsigned int offset);
int dot_prod(uint64_t A, uint64_t B);
uint64_t insertsl(uint64_t mask, uint64_t sm, uint64_t lg);
uint64_t decrypt_one_round( uint64_t in, uint64_t subkey, uint16_t Rounds, _Bool Roundwise );
int hamming_weight(uint64_t x);



//----------------------------------
// Start of code
//----------------------------------
int main()
{
    // Define the input/output masks and other parameters for linear analysis
	//#######################################################################
	int num_of_rounds = 8;
	uint64_t inputdiff = 0x0000000009000900;
	uint64_t outputmask = 0x0001000000010001;
	int DEG = 27;
	int N = 4;
    //#######################################################################
    uint64_t N1 = 1ULL << DEG; // Number of queries:  N1 = 2^(DEG)
	unsigned int initial_seed;
	uint64_t CORR;
	uint64_t sum = 0;

    clock_t clock_timer;
    clock_timer = clock();
	for (int num_of_experiments = 0; num_of_experiments < N; num_of_experiments++)
	{
		initial_seed = init_prng(143);
		double corr_log = 0;
		uint64_t counter0 = 0ULL;
		uint64_t counter1 = 0ULL;
		uint64_t key_high = 0ULL;
		uint64_t key_low = 0ULL;
		uint64_t plaintext1 = 0ULL;
		uint64_t plaintext2 = 0ULL;
		uint64_t ciphertext1 = 0ULL;
		uint64_t ciphertext2 = 0ULL;
		uint64_t *subkey;
		key_high = 0ULL;
		key_low = 0ULL;
		for(int i = 0; i < 16; i++) 
		{
			key_high |= (uint64_t)(rand() & 0xf) << i*4;
			key_low |= (uint64_t)(rand() & 0xf) << i*4;
		}			
		// Perform the key-schedule
		subkey = key_schedule(key_high, key_low, num_of_rounds + 1, 0, 0);
		for (uint64_t loopcnt = 0; loopcnt < N1; loopcnt++)
		{
			plaintext1 = 0ULL;
			for(int i = 0; i < 16; i++) plaintext1 |= (uint64_t)(rand() & 0xf) << (i*4);
			plaintext2 = plaintext1 ^ inputdiff;
			ciphertext1 = encrypt(plaintext1, subkey, num_of_rounds, 0);
			ciphertext2 = encrypt(plaintext2, subkey, num_of_rounds, 0);
			if (dot_prod(ciphertext1, outputmask) == dot_prod(ciphertext2, outputmask))
				counter0 += 1;
			else
				counter1 += 1;
		}
		if (counter0 > counter1)
			CORR = counter0 - counter1;
		else
			CORR = counter1 - counter0;
		printf("Exp No. %d \t Initial seed: 0x%X\n", num_of_experiments, initial_seed);
		printf("%s: %0.4f\n", "time on clock", (double)(clock() - clock_timer) / CLOCKS_PER_SEC);
		corr_log = (log(CORR) / log(2)) - DEG;
		printf("Correlation = 2^(%0.2f)\n", corr_log);
		printf("#############################################################\n");
		sum += CORR;
	}
	double corr = 0;
	corr = ((log(sum) - log(N)) / log(2)) - DEG;
	printf("Average correlation: 2^(%0.2f)\n", corr);
	printf("#############################################################\n");	
}  

//----------------------------------
// Init PRNG
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
// Dot product
//----------------------------------

int dot_prod(uint64_t A, uint64_t B)
{
	int output = 0;
	uint64_t C = (A & B);    
    while(C > 0){
        output ^= (C & 0x1);
        C >>= 1;
    }
    return output;
}

//----------------------------------
// Compute Hamming weight
//----------------------------------

int hamming_weight(uint64_t x)
{
	int output = 0;
	while(x > 0)
	{
		output += x & 0x1;
		x >>= 1;
	}
	return output;
}

//----------------------------------
// Insert a small binary vector within a larger one
//----------------------------------

uint64_t insertsl(uint64_t mask, uint64_t sm, uint64_t lg)
{
	for (int i = 0; i < 64; i++)
	{
		if (((mask >> (64-i)) & 0x1) == 1)
		{
			lg &= 0xFFFFFFFFFFFFFFFE; // kill the lsb
			lg |= (sm & 0x1); // insert the lsb of sm in lsb of lg
			sm >>= 1;
		}
		lg = rotate1l_64(lg); // rotate lg 
	}
	return lg;
}

//----------------------------------
// Key Scheduling
//----------------------------------
uint64_t* key_schedule( uint64_t key_high, uint64_t key_low, uint16_t Rounds, _Bool KeySize80, _Bool Output )
{

	uint64_t temp64;
	uint64_t i;

    uint64_t *subkey = (uint64_t *)malloc(Rounds*sizeof(uint64_t));

	if(subkey != NULL)
	{
		if(Output) v_key_start();

		if (KeySize80)
		{
			key_low &= 0xFFFF;

			if(Output) v_k80_init(key_high, key_low);

			for ( i=0; i<Rounds; i++)
			{
				subkey[i] = key_high;

				//----------------------------------
				// Shift
				//----------------------------------
				temp64 = key_high;
				key_high <<= 61;
				key_high |= (key_low<<45);
				key_high |= (temp64>>19);
				key_low = (temp64>>3)&0xFFFF;

				if(Output && (i+2 <= Rounds) ) v_k80_shift(key_high, key_low);

				//----------------------------------
				// S-Box
				//----------------------------------
				temp64 = high4_64(key_high); //get highest nibble
				key_high &=	0x0FFFFFFFFFFFFFFF;	//kill highest nibble
				temp64 = Sbox[temp64];
				key_high |=	low4_64(temp64); //put new value to highest nibble (sbox)

				if(Output && (i+2 <= Rounds) ) v_k80_sbox(key_high, key_low);

				//----------------------------------
				// Round Salt
				//----------------------------------
				key_low ^= ( ( (i+1) & 0x01 ) << 15  );
				key_high ^= ( (i+1) >> 1 );

				if(Output && (i+2 <= Rounds) ) v_k80_round(key_high, key_low, i);

			}
		}
		else //128 Bit
		{
			if(Output) v_k128_init(key_high, key_low);

			for ( i=0; i<Rounds; i++)
				{
				subkey[i] = key_high;

				//----------------------------------
				// Shift
				//----------------------------------
				temp64 = high61_64(key_high);
				key_high <<= 61;
				key_high |= high61_64(key_low);
				key_low <<= 61;
				key_low |= temp64;

				if(Output && (i+2 <= Rounds) ) v_k128_shift(key_high, key_low);

				//----------------------------------
				// S-Box
				//----------------------------------
				temp64 = high4_64(key_high); //get highest nibble
				key_high &=	0x0FFFFFFFFFFFFFFF;	//kill highest nibble
				temp64 = Sbox[temp64];
				key_high |=	low4_64(temp64); //put new value to highest nibble (sbox)
				temp64 = high8to4_64(key_high);	//get 2. highest nibble
				key_high &=	0xF0FFFFFFFFFFFFFF;	//kill 2. highest nibble
				temp64 = Sbox[temp64];
				key_high |=low8to4_64(temp64);	//put new value to 2. highest nibble (sbox)

				if(Output && (i+2 <= Rounds) ) v_k128_sbox(key_high, key_low);

				//----------------------------------
				// Round Salt
				//----------------------------------
				key_low ^= ( ( (i+1) & 0x03 ) << 62 );	//add counter to lower key part
				key_high ^= (  (i+1)  >> 2 );	//add counter to higher key part

				if(Output && (i+2 <= Rounds) ) v_k128_round(key_high, key_low, i);
			}
		}
		if(Output) v_final();
	}
	else
	{
      printf("RAM problem!\n");
      exit(0);
	}
	return subkey;
}




//----------------------------------
// Encryption
//----------------------------------
uint64_t encrypt( uint64_t in, uint64_t *subkey, uint16_t Rounds, _Bool Roundwise )
{	//Start encryption

	#define out in
	uint16_t RoundNr;
	uint64_t text;

	if (Roundwise) v_enc_start(in);

	for (RoundNr=1; RoundNr <= Rounds; RoundNr++)
	{	//Start "for"
		uint16_t temp;
		#define SboxNr temp
		#define PBit temp

		if (Roundwise) v_roundstart(RoundNr, subkey[RoundNr-1]);


		//----------------------------------
		// Xor with roundkey
		//----------------------------------
		text = in ^ subkey[RoundNr-1];

		if (Roundwise) v_after_xor(text);


		//----------------------------------
		// S-Boxes
		//----------------------------------
		for ( SboxNr=0; SboxNr<16; SboxNr++ )
		{
			uint16_t SboxVal;

			SboxVal	=	text & 0x0F; //get lowest nibble
			text &=	0xFFFFFFFFFFFFFFF0;	//kill lowest nibble
			text |=	Sbox[SboxVal];	//put new value to lowest nibble (sbox)
			text = rotate4l_64(text); //next(rotate by one nibble)
		}

		if (Roundwise) v_after_s(text);


		//----------------------------------
		// P-Box
		//----------------------------------
		for ( PBit = 0, out = 0; PBit<64; PBit++ )
		{
			out = rotate1l_64(out);	//next(rotate by one bit)
			out |= ( ( text >> (63-Pbox[PBit]) ) & 1 );  //put new value to lowest bit (pbox)
		}

		if (Roundwise) v_after_p(in);

	} //End "for"

	text = in ^ subkey[RoundNr - 1];

	if (Roundwise) v_enc_final(text, subkey[RoundNr - 1]);

	return text;

} //End encryption


//----------------------------------
// Decryption
//----------------------------------
uint64_t decrypt( uint64_t in, uint64_t *subkey, uint16_t Rounds, _Bool Roundwise )
{	//Start decryption
	#define out in
	uint16_t RoundNr;
	uint64_t text = 0;


	if (Roundwise) v_dec_start(in);

	for ( RoundNr=1; RoundNr<=Rounds; RoundNr++)
	{	//Start "for"
		uint16_t temp;
		#define SboxNr temp
		#define PBit temp

		if (Roundwise) v_roundstart(RoundNr, subkey[Rounds-RoundNr]);


		//----------------------------------
		// Xor with roundkey
		//----------------------------------
		text = in ^ subkey[Rounds-RoundNr];

		if (Roundwise) v_after_xor(text);


		//----------------------------------
		// P-Box
		//----------------------------------
		for ( PBit = 0, out = 0; PBit<64; PBit++ )
		{
			out = rotate1l_64(out);	//next(rotate by one bit)
			out |= ( ( text >> (63-PboxInv[PBit]) ) & 1 ); //put new value to lowest bit (pbox)
		}

		if (Roundwise) v_after_p(out);


		//----------------------------------
		// S-Boxes
		//----------------------------------
		for ( SboxNr=0; SboxNr<16; SboxNr++ )
		{
			uint16_t SboxVal;

			SboxVal	=	out & 0x0F;
			out &=	0xFFFFFFFFFFFFFFF0;
			out |=	SboxInv[SboxVal];
			out = rotate4l_64(out);
		}

		if (Roundwise) v_after_s(out);

	} //End "for"

	if (Roundwise) v_final();

	return text;

} //End decryption


uint64_t decrypt_one_round( uint64_t in, uint64_t subkey, uint16_t Rounds, _Bool Roundwise )
{	//Start decryption
	#define out in
	uint16_t RoundNr;
	uint64_t text = 0;

	for ( RoundNr=1; RoundNr<=Rounds; RoundNr++)
	{	//Start "for"
		uint16_t temp;
		#define SboxNr temp
		#define PBit temp

		//----------------------------------
		// Xor with roundkey
		//----------------------------------
		text = in ^ subkey;

		//----------------------------------
		// P-Box
		//----------------------------------
		for ( PBit = 0, out = 0; PBit<64; PBit++ )
		{
			out = rotate1l_64(out);	//next(rotate by one bit)
			out |= ( ( text >> (63-PboxInv[PBit]) ) & 1 ); //put new value to lowest bit (pbox)
		}

		//----------------------------------
		// S-Boxes
		//----------------------------------
		for ( SboxNr=0; SboxNr<16; SboxNr++ )
		{
			uint16_t SboxVal;

			SboxVal	=	out & 0x0F;
			out &=	0xFFFFFFFFFFFFFFF0;
			out |=	SboxInv[SboxVal];
			out = rotate4l_64(out);
		}

	} //End "for"

	return out;

} //End decryption
