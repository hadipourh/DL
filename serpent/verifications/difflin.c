/*
#################################################################################
#################################################################################
#################################################################################
Modified on Sep 19, 2023, to compute the correlation
of differential-linear distinguishers of SERPENT block cipher.
Author: Hosein Hadipour

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

#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/random.h>
#include <time.h>
#include <math.h>

typedef uint32_t u32;
typedef unsigned char byte;

/* Number of rounds per Serpent encrypt/decrypt operation.  */
#define ROUNDS 32

/* Magic number, used during generating of the subkeys.  */
#define PHI 0x9E3779B9

/* Serpent works on 128 bit blocks.  */
typedef u32 serpent_block_t[4];

/* Serpent key, provided by the user.  If the original key is shorter
   than 256 bits, it is padded.  */
typedef u32 serpent_key_t[8];

/* The key schedule consists of 33 128 bit subkeys.  */
typedef u32 serpent_subkeys_t[ROUNDS + 1][4];

/* A Serpent context.  */
typedef struct serpent_context
{
  serpent_subkeys_t keys;	/* Generated subkeys.  */
} serpent_context_t;

/* Functions for loading and storing unaligned u32 values of different
   endianness.  */
static inline u32 buf_get_be32(const void *_buf)
{
  const byte *in = _buf;
  return ((u32)in[0] << 24) | ((u32)in[1] << 16) | \
		 ((u32)in[2] << 8) | (u32)in[3];
}

static inline u32 buf_get_le32(const void *_buf)
{
  const byte *in = _buf;
  return ((u32)in[3] << 24) | ((u32)in[2] << 16) | \
		 ((u32)in[1] << 8) | (u32)in[0];
}

static inline void buf_put_be32(void *_buf, u32 val)
{
  byte *out = _buf;
  out[0] = val >> 24;
  out[1] = val >> 16;
  out[2] = val >> 8;
  out[3] = val;
}

static inline void buf_put_le32(void *_buf, u32 val)
{
  byte *out = _buf;
  out[3] = val >> 24;
  out[2] = val >> 16;
  out[1] = val >> 8;
  out[0] = val;
}

static inline u32 rol(u32 x, int n)
{
	return ( (x << (n&(32-1))) | (x >> ((32-n)&(32-1))) );
}

static inline u32 ror(u32 x, int n)
{
	return ( (x >> (n&(32-1))) | (x << ((32-n)&(32-1))) );
}

// Function to securely clear memory
void secure_memzero(void *ptr, size_t size) {
	volatile unsigned char *p = ptr;
	while (size--) {
		*p++ = 0;
	}
}

/*
 * These are the S-Boxes of Serpent from following research paper.
 *
 *  D. A. Osvik, “Speeding up Serpent,” in Third AES Candidate Conference,
 *   (New York, New York, USA), p. 317–329, National Institute of Standards and
 *   Technology, 2000.
 *
 * Paper is also available at: http://www.ii.uib.no/~osvik/pub/aes3.pdf
 *
 */

#define SBOX0(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r3 ^= r0; r4 =  r1; \
	r1 &= r3; r4 ^= r2; \
	r1 ^= r0; r0 |= r3; \
	r0 ^= r4; r4 ^= r3; \
	r3 ^= r2; r2 |= r1; \
	r2 ^= r4; r4 = ~r4; \
	r4 |= r1; r1 ^= r3; \
	r1 ^= r4; r3 |= r0; \
	r1 ^= r3; r4 ^= r3; \
	\
	w = r1; x = r4; y = r2; z = r0; \
  }

#define SBOX0_INVERSE(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r2 = ~r2; r4 =  r1; \
	r1 |= r0; r4 = ~r4; \
	r1 ^= r2; r2 |= r4; \
	r1 ^= r3; r0 ^= r4; \
	r2 ^= r0; r0 &= r3; \
	r4 ^= r0; r0 |= r1; \
	r0 ^= r2; r3 ^= r4; \
	r2 ^= r1; r3 ^= r0; \
	r3 ^= r1; \
	r2 &= r3; \
	r4 ^= r2; \
	\
	w = r0; x = r4; y = r1; z = r3; \
  }

#define SBOX1(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r0 = ~r0; r2 = ~r2; \
	r4 =  r0; r0 &= r1; \
	r2 ^= r0; r0 |= r3; \
	r3 ^= r2; r1 ^= r0; \
	r0 ^= r4; r4 |= r1; \
	r1 ^= r3; r2 |= r0; \
	r2 &= r4; r0 ^= r1; \
	r1 &= r2; \
	r1 ^= r0; r0 &= r2; \
	r0 ^= r4; \
	\
	w = r2; x = r0; y = r3; z = r1; \
  }

#define SBOX1_INVERSE(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r4 =  r1; r1 ^= r3; \
	r3 &= r1; r4 ^= r2; \
	r3 ^= r0; r0 |= r1; \
	r2 ^= r3; r0 ^= r4; \
	r0 |= r2; r1 ^= r3; \
	r0 ^= r1; r1 |= r3; \
	r1 ^= r0; r4 = ~r4; \
	r4 ^= r1; r1 |= r0; \
	r1 ^= r0; \
	r1 |= r4; \
	r3 ^= r1; \
	\
	w = r4; x = r0; y = r3; z = r2; \
  }

#define SBOX2(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r4 =  r0; r0 &= r2; \
	r0 ^= r3; r2 ^= r1; \
	r2 ^= r0; r3 |= r4; \
	r3 ^= r1; r4 ^= r2; \
	r1 =  r3; r3 |= r4; \
	r3 ^= r0; r0 &= r1; \
	r4 ^= r0; r1 ^= r3; \
	r1 ^= r4; r4 = ~r4; \
	\
	w = r2; x = r3; y = r1; z = r4; \
  }

#define SBOX2_INVERSE(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r2 ^= r3; r3 ^= r0; \
	r4 =  r3; r3 &= r2; \
	r3 ^= r1; r1 |= r2; \
	r1 ^= r4; r4 &= r3; \
	r2 ^= r3; r4 &= r0; \
	r4 ^= r2; r2 &= r1; \
	r2 |= r0; r3 = ~r3; \
	r2 ^= r3; r0 ^= r3; \
	r0 &= r1; r3 ^= r4; \
	r3 ^= r0; \
	\
	w = r1; x = r4; y = r2; z = r3; \
  }

#define SBOX3(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r4 =  r0; r0 |= r3; \
	r3 ^= r1; r1 &= r4; \
	r4 ^= r2; r2 ^= r3; \
	r3 &= r0; r4 |= r1; \
	r3 ^= r4; r0 ^= r1; \
	r4 &= r0; r1 ^= r3; \
	r4 ^= r2; r1 |= r0; \
	r1 ^= r2; r0 ^= r3; \
	r2 =  r1; r1 |= r3; \
	r1 ^= r0; \
	\
	w = r1; x = r2; y = r3; z = r4; \
  }

#define SBOX3_INVERSE(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r4 =  r2; r2 ^= r1; \
	r0 ^= r2; r4 &= r2; \
	r4 ^= r0; r0 &= r1; \
	r1 ^= r3; r3 |= r4; \
	r2 ^= r3; r0 ^= r3; \
	r1 ^= r4; r3 &= r2; \
	r3 ^= r1; r1 ^= r0; \
	r1 |= r2; r0 ^= r3; \
	r1 ^= r4; \
	r0 ^= r1; \
	\
	w = r2; x = r1; y = r3; z = r0; \
  }

#define SBOX4(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r1 ^= r3; r3 = ~r3; \
	r2 ^= r3; r3 ^= r0; \
	r4 =  r1; r1 &= r3; \
	r1 ^= r2; r4 ^= r3; \
	r0 ^= r4; r2 &= r4; \
	r2 ^= r0; r0 &= r1; \
	r3 ^= r0; r4 |= r1; \
	r4 ^= r0; r0 |= r3; \
	r0 ^= r2; r2 &= r3; \
	r0 = ~r0; r4 ^= r2; \
	\
	w = r1; x = r4; y = r0; z = r3; \
  }

#define SBOX4_INVERSE(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r4 =  r2; r2 &= r3; \
	r2 ^= r1; r1 |= r3; \
	r1 &= r0; r4 ^= r2; \
	r4 ^= r1; r1 &= r2; \
	r0 = ~r0; r3 ^= r4; \
	r1 ^= r3; r3 &= r0; \
	r3 ^= r2; r0 ^= r1; \
	r2 &= r0; r3 ^= r0; \
	r2 ^= r4; \
	r2 |= r3; r3 ^= r0; \
	r2 ^= r1; \
	\
	w = r0; x = r3; y = r2; z = r4; \
  }

#define SBOX5(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r0 ^= r1; r1 ^= r3; \
	r3 = ~r3; r4 =  r1; \
	r1 &= r0; r2 ^= r3; \
	r1 ^= r2; r2 |= r4; \
	r4 ^= r3; r3 &= r1; \
	r3 ^= r0; r4 ^= r1; \
	r4 ^= r2; r2 ^= r0; \
	r0 &= r3; r2 = ~r2; \
	r0 ^= r4; r4 |= r3; \
	r2 ^= r4; \
	\
	w = r1; x = r3; y = r0; z = r2; \
  }

#define SBOX5_INVERSE(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r1 = ~r1; r4 =  r3; \
	r2 ^= r1; r3 |= r0; \
	r3 ^= r2; r2 |= r1; \
	r2 &= r0; r4 ^= r3; \
	r2 ^= r4; r4 |= r0; \
	r4 ^= r1; r1 &= r2; \
	r1 ^= r3; r4 ^= r2; \
	r3 &= r4; r4 ^= r1; \
	r3 ^= r4; r4 = ~r4; \
	r3 ^= r0; \
	\
	w = r1; x = r4; y = r3; z = r2; \
  }

#define SBOX6(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r2 = ~r2; r4 =  r3; \
	r3 &= r0; r0 ^= r4; \
	r3 ^= r2; r2 |= r4; \
	r1 ^= r3; r2 ^= r0; \
	r0 |= r1; r2 ^= r1; \
	r4 ^= r0; r0 |= r3; \
	r0 ^= r2; r4 ^= r3; \
	r4 ^= r0; r3 = ~r3; \
	r2 &= r4; \
	r2 ^= r3; \
	\
	w = r0; x = r1; y = r4; z = r2; \
  }

#define SBOX6_INVERSE(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r0 ^= r2; r4 =  r2; \
	r2 &= r0; r4 ^= r3; \
	r2 = ~r2; r3 ^= r1; \
	r2 ^= r3; r4 |= r0; \
	r0 ^= r2; r3 ^= r4; \
	r4 ^= r1; r1 &= r3; \
	r1 ^= r0; r0 ^= r3; \
	r0 |= r2; r3 ^= r1; \
	r4 ^= r0; \
	\
	w = r1; x = r2; y = r4; z = r3; \
  }

#define SBOX7(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r4 =  r1; r1 |= r2; \
	r1 ^= r3; r4 ^= r2; \
	r2 ^= r1; r3 |= r4; \
	r3 &= r0; r4 ^= r2; \
	r3 ^= r1; r1 |= r4; \
	r1 ^= r0; r0 |= r4; \
	r0 ^= r2; r1 ^= r4; \
	r2 ^= r1; r1 &= r0; \
	r1 ^= r4; r2 = ~r2; \
	r2 |= r0; \
	r4 ^= r2; \
	\
	w = r4; x = r3; y = r1; z = r0; \
  }

#define SBOX7_INVERSE(r0, r1, r2, r3, w, x, y, z) \
  { \
	u32 r4; \
	\
	r4 =  r2; r2 ^= r0; \
	r0 &= r3; r4 |= r3; \
	r2 = ~r2; r3 ^= r1; \
	r1 |= r0; r0 ^= r2; \
	r2 &= r4; r3 &= r4; \
	r1 ^= r2; r2 ^= r0; \
	r0 |= r2; r4 ^= r1; \
	r0 ^= r3; r3 ^= r4; \
	r4 |= r0; r3 ^= r2; \
	r4 ^= r2; \
	\
	w = r3; x = r0; y = r1; z = r4; \
  }

/* XOR BLOCK1 into BLOCK0.  */
#define BLOCK_XOR(block0, block1) \
  {                               \
	block0[0] ^= block1[0];       \
	block0[1] ^= block1[1];       \
	block0[2] ^= block1[2];       \
	block0[3] ^= block1[3];       \
  }

/* Copy BLOCK_SRC to BLOCK_DST.  */
#define BLOCK_COPY(block_dst, block_src) \
  {                                      \
	block_dst[0] = block_src[0];         \
	block_dst[1] = block_src[1];         \
	block_dst[2] = block_src[2];         \
	block_dst[3] = block_src[3];         \
  }

/* Apply SBOX number WHICH to to the block found in ARRAY0, writing
   the output to the block found in ARRAY1.  */
#define SBOX(which, array0, array1)                         \
  SBOX##which (array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]);

#define SBOX_CASE(which, array0, array1)                         \
  switch(which) {                                                \
	case 0: SBOX0(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 1: SBOX1(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 2: SBOX2(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 3: SBOX3(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 4: SBOX4(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 5: SBOX5(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 6: SBOX6(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 7: SBOX7(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
  }

/* Apply inverse SBOX number WHICH to to the block found in ARRAY0, writing
   the output to the block found in ARRAY1.  */
#define SBOX_INVERSE(which, array0, array1)                           \
  SBOX##which##_INVERSE (array0[0], array0[1], array0[2], array0[3],  \
						 array1[0], array1[1], array1[2], array1[3]);
						
#define SBOX_INVERSE_CASE(which, array0, array1)                           \
  switch(which) {                                           \
	case 0: SBOX0_INVERSE(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 1: SBOX1_INVERSE(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 2: SBOX2_INVERSE(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 3: SBOX3_INVERSE(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 4: SBOX4_INVERSE(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 5: SBOX5_INVERSE(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 6: SBOX6_INVERSE(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
	case 7: SBOX7_INVERSE(array0[0], array0[1], array0[2], array0[3],  \
			   array1[0], array1[1], array1[2], array1[3]); break; \
  }

/* Apply the linear transformation to BLOCK.  */
#define LINEAR_TRANSFORMATION(block)                  \
  {                                                   \
	block[0] = rol (block[0], 13);                    \
	block[2] = rol (block[2], 3);                     \
	block[1] = block[1] ^ block[0] ^ block[2];        \
	block[3] = block[3] ^ block[2] ^ (block[0] << 3); \
	block[1] = rol (block[1], 1);                     \
	block[3] = rol (block[3], 7);                     \
	block[0] = block[0] ^ block[1] ^ block[3];        \
	block[2] = block[2] ^ block[3] ^ (block[1] << 7); \
	block[0] = rol (block[0], 5);                     \
	block[2] = rol (block[2], 22);                    \
  }

/* Apply the inverse linear transformation to BLOCK.  */
#define LINEAR_TRANSFORMATION_INVERSE(block)          \
  {                                                   \
	block[2] = ror (block[2], 22);                    \
	block[0] = ror (block[0] , 5);                    \
	block[2] = block[2] ^ block[3] ^ (block[1] << 7); \
	block[0] = block[0] ^ block[1] ^ block[3];        \
	block[3] = ror (block[3], 7);                     \
	block[1] = ror (block[1], 1);                     \
	block[3] = block[3] ^ block[2] ^ (block[0] << 3); \
	block[1] = block[1] ^ block[0] ^ block[2];        \
	block[2] = ror (block[2], 3);                     \
	block[0] = ror (block[0], 13);                    \
  }

/* Apply a Serpent round to BLOCK, using the SBOX number WHICH and the
   subkeys contained in SUBKEYS.  Use BLOCK_TMP as temporary storage.
   This macro increments `round'.  */
#define ROUND(which, subkeys, block, block_tmp) \
  {                                             \
	BLOCK_XOR (block, subkeys[round]);          \
	round++;                                    \
	SBOX (which, block, block_tmp);             \
	LINEAR_TRANSFORMATION (block_tmp);          \
	BLOCK_COPY (block, block_tmp);              \
  }

#define ROUND_CASE(which, subkeys, block, block_tmp, rno) \
  {                                             \
	BLOCK_XOR (block, subkeys[rno]);            \
	SBOX_CASE(which, block, block_tmp);         \
	LINEAR_TRANSFORMATION (block_tmp);          \
	BLOCK_COPY (block, block_tmp);              \
  }

/* Apply the last Serpent round to BLOCK, using the SBOX number WHICH
   and the subkeys contained in SUBKEYS.  Use BLOCK_TMP as temporary
   storage.  The result will be stored in BLOCK_TMP.  This macro
   increments `round'.  */
#define ROUND_LAST(which, subkeys, block, block_tmp) \
  {                                                  \
	BLOCK_XOR (block, subkeys[round]);               \
	round++;                                         \
	SBOX (which, block, block_tmp);                  \
	BLOCK_XOR (block_tmp, subkeys[round]);           \
	round++;                                         \
  }

/* Apply an inverse Serpent round to BLOCK, using the SBOX number
   WHICH and the subkeys contained in SUBKEYS.  Use BLOCK_TMP as
   temporary storage.  This macro increments `round'.  */
#define ROUND_INVERSE(which, subkey, block, block_tmp) \
  {                                                    \
	LINEAR_TRANSFORMATION_INVERSE (block);             \
	SBOX_INVERSE (which, block, block_tmp);            \
	BLOCK_XOR (block_tmp, subkey[round]);              \
	round--;                                           \
	BLOCK_COPY (block, block_tmp);                     \
  }

#define ROUND_INVERSE_CASE(which, subkey, block, block_tmp, rno) \
  {                                                    \
	LINEAR_TRANSFORMATION_INVERSE (block);             \
	SBOX_INVERSE_CASE (which, block, block_tmp);       \
	BLOCK_XOR (block_tmp, subkey[rno]);                \
	BLOCK_COPY (block, block_tmp);                     \
  }

/* Apply the first Serpent round to BLOCK, using the SBOX number WHICH
   and the subkeys contained in SUBKEYS.  Use BLOCK_TMP as temporary
   storage.  The result will be stored in BLOCK_TMP.  This macro
   increments `round'.  */
#define ROUND_FIRST_INVERSE(which, subkeys, block, block_tmp) \
  {                                                           \
	BLOCK_XOR (block, subkeys[round]);                        \
	round--;                                                  \
	SBOX_INVERSE (which, block, block_tmp);                   \
	BLOCK_XOR (block_tmp, subkeys[round]);                    \
	round--;                                                  \
  }

/* Convert the user provided key KEY of KEY_LENGTH bytes into the
   internally used format.  */
static void
serpent_key_prepare (const byte *key, unsigned int key_length,
			 serpent_key_t key_prepared)
{
  int i;

  /* Copy key.  */
  key_length /= 4;
  for (i = 0; i < key_length; i++)
	key_prepared[i] = buf_get_le32 (key + i * 4);

  if (i < 8)
	{
	  /* Key must be padded according to the Serpent
	 specification.  */
	  key_prepared[i] = 0x00000001;

	  for (i++; i < 8; i++)
	key_prepared[i] = 0;
	}
}

/* Derive the 33 subkeys from KEY and store them in SUBKEYS.  */
static void
serpent_subkeys_generate (serpent_key_t key, serpent_subkeys_t subkeys)
{
  u32 w[8];		/* The `prekey'.  */
  u32 ws[4];
  u32 wt[4];

  /* Initialize with key values.  */
  w[0] = key[0];
  w[1] = key[1];
  w[2] = key[2];
  w[3] = key[3];
  w[4] = key[4];
  w[5] = key[5];
  w[6] = key[6];
  w[7] = key[7];

  /* Expand to intermediate key using the affine recurrence.  */
#define EXPAND_KEY4(wo, r)                                                     \
  wo[0] = w[(r+0)%8] =                                                         \
	rol (w[(r+0)%8] ^ w[(r+3)%8] ^ w[(r+5)%8] ^ w[(r+7)%8] ^ PHI ^ (r+0), 11); \
  wo[1] = w[(r+1)%8] =                                                         \
	rol (w[(r+1)%8] ^ w[(r+4)%8] ^ w[(r+6)%8] ^ w[(r+0)%8] ^ PHI ^ (r+1), 11); \
  wo[2] = w[(r+2)%8] =                                                         \
	rol (w[(r+2)%8] ^ w[(r+5)%8] ^ w[(r+7)%8] ^ w[(r+1)%8] ^ PHI ^ (r+2), 11); \
  wo[3] = w[(r+3)%8] =                                                         \
	rol (w[(r+3)%8] ^ w[(r+6)%8] ^ w[(r+0)%8] ^ w[(r+2)%8] ^ PHI ^ (r+3), 11);

#define EXPAND_KEY(r)       \
	EXPAND_KEY4(ws, (r));     \
	EXPAND_KEY4(wt, (r + 4));

	/* Calculate subkeys via S-Boxes, in bitslice mode.  */
	EXPAND_KEY (0); SBOX (3, ws, subkeys[0]); SBOX (2, wt, subkeys[1]);
	EXPAND_KEY (8); SBOX (1, ws, subkeys[2]); SBOX (0, wt, subkeys[3]);
	EXPAND_KEY (16); SBOX (7, ws, subkeys[4]); SBOX (6, wt, subkeys[5]);
	EXPAND_KEY (24); SBOX (5, ws, subkeys[6]); SBOX (4, wt, subkeys[7]);
	EXPAND_KEY (32); SBOX (3, ws, subkeys[8]); SBOX (2, wt, subkeys[9]);
	EXPAND_KEY (40); SBOX (1, ws, subkeys[10]); SBOX (0, wt, subkeys[11]);
	EXPAND_KEY (48); SBOX (7, ws, subkeys[12]); SBOX (6, wt, subkeys[13]);
	EXPAND_KEY (56); SBOX (5, ws, subkeys[14]); SBOX (4, wt, subkeys[15]);
	EXPAND_KEY (64); SBOX (3, ws, subkeys[16]); SBOX (2, wt, subkeys[17]);
	EXPAND_KEY (72); SBOX (1, ws, subkeys[18]); SBOX (0, wt, subkeys[19]);
	EXPAND_KEY (80); SBOX (7, ws, subkeys[20]); SBOX (6, wt, subkeys[21]);
	EXPAND_KEY (88); SBOX (5, ws, subkeys[22]); SBOX (4, wt, subkeys[23]);
	EXPAND_KEY (96); SBOX (3, ws, subkeys[24]); SBOX (2, wt, subkeys[25]);
	EXPAND_KEY (104); SBOX (1, ws, subkeys[26]); SBOX (0, wt, subkeys[27]);
	EXPAND_KEY (112); SBOX (7, ws, subkeys[28]); SBOX (6, wt, subkeys[29]);
	EXPAND_KEY (120); SBOX (5, ws, subkeys[30]); SBOX (4, wt, subkeys[31]);
	EXPAND_KEY4 (ws, 128); SBOX (3, ws, subkeys[32]);

	secure_memzero (ws, sizeof (ws));
	secure_memzero (wt, sizeof (wt));
	secure_memzero (w, sizeof (w));
}

/* Initialize CONTEXT with the key KEY of KEY_LENGTH bits.  */
static int
serpent_setkey_internal (serpent_context_t *context,
			 const byte *key, unsigned int key_length)
{
  serpent_key_t key_prepared;

  if (key_length > 32)
	return 0;

  serpent_key_prepare (key, key_length, key_prepared);
  serpent_subkeys_generate (key_prepared, context->keys);
  secure_memzero (key_prepared, sizeof(key_prepared));
  return 0;
}

static void
serpent_encrypt_internal (serpent_context_t *context,
			  const byte *input, byte *output)
{
	serpent_block_t b, b_next;
	int round = 0;

	b[0] = buf_get_le32 (input + 0);
	b[1] = buf_get_le32 (input + 4);
	b[2] = buf_get_le32 (input + 8);
	b[3] = buf_get_le32 (input + 12);

	ROUND (0, context->keys, b, b_next);
	ROUND (1, context->keys, b, b_next);
	ROUND (2, context->keys, b, b_next);
	ROUND (3, context->keys, b, b_next);
	ROUND (4, context->keys, b, b_next);
	ROUND (5, context->keys, b, b_next);
	ROUND (6, context->keys, b, b_next);
	ROUND (7, context->keys, b, b_next);
	ROUND (0, context->keys, b, b_next);
	ROUND (1, context->keys, b, b_next);
	ROUND (2, context->keys, b, b_next);
	ROUND (3, context->keys, b, b_next);
	ROUND (4, context->keys, b, b_next);
	ROUND (5, context->keys, b, b_next);
	ROUND (6, context->keys, b, b_next);
	ROUND (7, context->keys, b, b_next);
	ROUND (0, context->keys, b, b_next);
	ROUND (1, context->keys, b, b_next);
	ROUND (2, context->keys, b, b_next);
	ROUND (3, context->keys, b, b_next);
	ROUND (4, context->keys, b, b_next);
	ROUND (5, context->keys, b, b_next);
	ROUND (6, context->keys, b, b_next);
	ROUND (7, context->keys, b, b_next);
	ROUND (0, context->keys, b, b_next);
	ROUND (1, context->keys, b, b_next);
	ROUND (2, context->keys, b, b_next);
	ROUND (3, context->keys, b, b_next);
	ROUND (4, context->keys, b, b_next);
	ROUND (5, context->keys, b, b_next);
	ROUND (6, context->keys, b, b_next);

	ROUND_LAST (7, context->keys, b, b_next);

	buf_put_le32 (output + 0, b_next[0]);
	buf_put_le32 (output + 4, b_next[1]);
	buf_put_le32 (output + 8, b_next[2]);
	buf_put_le32 (output + 12, b_next[3]);
}

static void
serpent_decrypt_internal (serpent_context_t *context,
			  const byte *input, byte *output)
{
	serpent_block_t b, b_next;
	int round = ROUNDS;

	b_next[0] = buf_get_le32 (input + 0);
	b_next[1] = buf_get_le32 (input + 4);
	b_next[2] = buf_get_le32 (input + 8);
	b_next[3] = buf_get_le32 (input + 12);

	ROUND_FIRST_INVERSE (7, context->keys, b_next, b);

	ROUND_INVERSE (6, context->keys, b, b_next);
	ROUND_INVERSE (5, context->keys, b, b_next);
	ROUND_INVERSE (4, context->keys, b, b_next);
	ROUND_INVERSE (3, context->keys, b, b_next);
	ROUND_INVERSE (2, context->keys, b, b_next);
	ROUND_INVERSE (1, context->keys, b, b_next);
	ROUND_INVERSE (0, context->keys, b, b_next);
	ROUND_INVERSE (7, context->keys, b, b_next);
	ROUND_INVERSE (6, context->keys, b, b_next);
	ROUND_INVERSE (5, context->keys, b, b_next);
	ROUND_INVERSE (4, context->keys, b, b_next);
	ROUND_INVERSE (3, context->keys, b, b_next);
	ROUND_INVERSE (2, context->keys, b, b_next);
	ROUND_INVERSE (1, context->keys, b, b_next);
	ROUND_INVERSE (0, context->keys, b, b_next);
	ROUND_INVERSE (7, context->keys, b, b_next);
	ROUND_INVERSE (6, context->keys, b, b_next);
	ROUND_INVERSE (5, context->keys, b, b_next);
	ROUND_INVERSE (4, context->keys, b, b_next);
	ROUND_INVERSE (3, context->keys, b, b_next);
	ROUND_INVERSE (2, context->keys, b, b_next);
	ROUND_INVERSE (1, context->keys, b, b_next);
	ROUND_INVERSE (0, context->keys, b, b_next);
	ROUND_INVERSE (7, context->keys, b, b_next);
	ROUND_INVERSE (6, context->keys, b, b_next);
	ROUND_INVERSE (5, context->keys, b, b_next);
	ROUND_INVERSE (4, context->keys, b, b_next);
	ROUND_INVERSE (3, context->keys, b, b_next);
	ROUND_INVERSE (2, context->keys, b, b_next);
	ROUND_INVERSE (1, context->keys, b, b_next);
	ROUND_INVERSE (0, context->keys, b, b_next);

	buf_put_le32 (output + 0, b_next[0]);
	buf_put_le32 (output + 4, b_next[1]);
	buf_put_le32 (output + 8, b_next[2]);
	buf_put_le32 (output + 12, b_next[3]);
}

// static unsigned int
// serpent_encrypt (void *ctx, byte *buffer_out, const byte *buffer_in)
// {
// 	serpent_context_t *context = ctx;

// 	serpent_encrypt_internal (context, buffer_in, buffer_out);
// 	return /*burn_stack*/ (2 * sizeof (serpent_block_t));
// }

// static unsigned int
// serpent_decrypt (void *ctx, byte *buffer_out, const byte *buffer_in)
// {
// 	serpent_context_t *context = ctx;

// 	serpent_decrypt_internal (context, buffer_in, buffer_out);
// 	return /*burn_stack*/ (2 * sizeof (serpent_block_t));
// }

void init_prng(unsigned int offset) {
    unsigned int initial_seed = 0;
    ssize_t temp;
    temp = getrandom(&initial_seed, sizeof(initial_seed), 0);
    if (temp == -1) perror("error!");
    initial_seed += offset;
    srand(initial_seed);
    printf("[+] PRNG initialized to 0x%08X\n", initial_seed);    
}

static void
encrypt (serpent_context_t *context,
			  const byte *input, byte *output, int offset, int nr)
{
	serpent_block_t b, b_next;
	b_next[0] = 0;
	b_next[1] = 0;
	b_next[2] = 0;
	b_next[3] = 0;

	b[0] = buf_get_le32 (input + 0);
	b[1] = buf_get_le32 (input + 4);
	b[2] = buf_get_le32 (input + 8);
	b[3] = buf_get_le32 (input + 12);

	for(int r = 0; r < nr; r++)
	{
		int RN = (offset + r)%8;
		ROUND_CASE(RN, context->keys, b, b_next, r);
	}

	buf_put_le32 (output + 0, b_next[0]);
	buf_put_le32 (output + 4, b_next[1]);
	buf_put_le32 (output + 8, b_next[2]);
	buf_put_le32 (output + 12, b_next[3]);
}

static void
decrypt (serpent_context_t *context,
			  const byte *input, byte *output, int offset, int nr)
{
	serpent_block_t b, b_next;
	b_next[0] = 0;
	b_next[1] = 0;
	b_next[2] = 0;
	b_next[3] = 0;

	b[0] = buf_get_le32 (input + 0);
	b[1] = buf_get_le32 (input + 4);
	b[2] = buf_get_le32 (input + 8);
	b[3] = buf_get_le32 (input + 12);
	
	for(int r = nr - 1; r > -1; r--)
	{
		int RN = (offset + r)%8;
		ROUND_INVERSE_CASE(RN, context->keys, b, b_next, r);
	}

	buf_put_le32 (output + 0, b_next[0]);
	buf_put_le32 (output + 4, b_next[1]);
	buf_put_le32 (output + 8, b_next[2]);
	buf_put_le32 (output + 12, b_next[3]);
}

void print_state(byte state[], int length)
{
	for (int i = length - 1; i > -1; i--)
		printf("%02x", state[i]);
	printf("\n");
}

byte dot_product(byte a[16], byte b[16])
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
        dx[15 - i] = (unsigned char)(strtol(hex, NULL, 16) & 0xff);
    }
}

uint64_t dldistinguisher(byte* master_key, byte* input_difference, byte* output_mask, int offset, int num_of_rounds, uint64_t num_of_tries)
{
	serpent_context_t ctx;
	serpent_setkey_internal (&ctx, master_key, 32);
	byte p1[16];
	byte p2[16];
	byte c1[16];
	byte c2[16];
	uint64_t counter_0 = 0;
	uint64_t counter_1 = 0;
	for(uint64_t i = 0; i < num_of_tries; i++){
		for(int j = 0; j < 16; j++){
			p1[j] = rand() & 0xff;
			p2[j] = p1[j] ^ input_difference[j];
		}
		encrypt(&ctx, p1, c1, offset, num_of_rounds);
		encrypt(&ctx, p2, c2, offset, num_of_rounds);
		if(dot_product(c1, output_mask) == dot_product(c2, output_mask))
			counter_0++;
		else
			counter_1++;
	}
	uint64_t absolute_correlation;
	if(counter_0 > counter_1)
		absolute_correlation = counter_0 - counter_1;
	else
		absolute_correlation = counter_1 - counter_0;
	return absolute_correlation;
}

int main(int argc, char *argv[])
{
	serpent_context_t ctx;
	byte key[32];
	byte plaintext[16];
	byte ciphertext[16];
	byte temp[16];
	for (int i = 0; i < 32; i++)
		key[i] = 0;
	for (int i = 0; i < 16; i++)
		plaintext[i] = 0;
	plaintext[15] = 0x80;
	serpent_setkey_internal (&ctx, key, 32);
	printf("key        : ");
	print_state(key, 32);
	printf("plaintext  : ");
	print_state(plaintext, 16);
	serpent_encrypt_internal (&ctx, plaintext, ciphertext);
	printf("ciphertext : ");
	print_state(ciphertext, 16);
	serpent_decrypt_internal (&ctx, ciphertext, temp);
	printf("decrypted  : ");
	print_state(temp, 16);
	if (memcmp(plaintext, temp, 16) == 0)
	{
		printf("decryption successful\n");
	} else
	{
		printf("decryption failed\n");
		return 0;
	}
	/*
	//#########################################################################################################
	Check the encrypt and decrypt functions
	*/
    printf("//#################################################################################################\n");
	printf("Check the encrypt and decrypt functions\n");
	printf("key        : ");
	print_state(key, 32);
	printf("plaintext  : ");
	print_state(plaintext, 16);
	encrypt(&ctx, plaintext, ciphertext, 1, 10);
	printf("ciphertext : ");
	print_state(ciphertext, 16);
	decrypt(&ctx, ciphertext, temp, 1, 10);
	printf("decrypted  : ");
	print_state(temp, 16);
	if (memcmp(plaintext, temp, 16) == 0)
		printf("decryption successful\n");
	else
	{
		printf("decryption failed\n");
		return 0;
	}
	printf("//#################################################################################################\n");
    // Check the distinguisher
    //#########################################################################################################
    //#########################################################################################################
    //#########################################################################################################
    int DEG1 = 0;
    int DEG2 = 25;
    uint64_t N1 = 1ULL << DEG1;
    uint64_t N2 = 1ULL << DEG2;
    int NUMBER_OF_EXPERIMENTS = 7;   // Number of independent experiments

    int NUMBER_OF_ROUNDS = 3;   // Number of rounds
	int offset = 4;   
    char DP_STR[] = "00000010040000004000000000000208";
    char LC_STR[] = "00100000000000000010000002000000";
    //##########################################################################################################
	// unsigned int task_id;
    // task_id = atoi(argv[1]);   
    init_prng(0); // we can feed init_prng with task_id
	// printf("Task ID: %d\n", task_id);
	
	byte input_difference[16];
	byte output_mask[16];
	convert_hexstr_to_statearray(DP_STR, input_difference);
    convert_hexstr_to_statearray(LC_STR, output_mask);
	printf("Input difference: ");
	print_state(input_difference, 16);
	printf("Output mask     : ");
	print_state(output_mask, 16);

	double sum = 0;
	for(int n = 0; n < NUMBER_OF_EXPERIMENTS; n++)
	{
		double num = 0;
		clock_t clock_timer;        
        clock_timer = clock();
		for(uint64_t i = 0; i < N1; i++)
		{
			byte key[32];
			for(int j = 0; j < 32; j++)
				key[j] = rand() & 0xff;
			num += dldistinguisher(key, input_difference, output_mask, offset, NUMBER_OF_ROUNDS, N2);
		}
		double elapsed_time = ((double) (clock() - clock_timer)) / CLOCKS_PER_SEC;
		printf("Execution time: %0.2f\n", elapsed_time);
		sum += num;
		double temp = log((double)N1) + log((double)N2);
		double avg_cr = (temp - log((double)num))/log(2.0);
		printf("\nCorrelation = 2^(-%0.4f)\n", avg_cr);
        printf("####################################\n");
	}
	double result = log((double)NUMBER_OF_EXPERIMENTS) + log((double)N1) + log((double)N2);
    double avg_cr = (result - log((double)sum))/log(2.0);
    printf("\nAverage correlation = 2^(-%0.4f)\n", avg_cr);
    printf("####################################\n");
	return 0;  
}
