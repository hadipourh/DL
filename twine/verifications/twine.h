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

Disclaimer: We acknowledge that the TWINE block cipher doesn't adhere to statistical assumptions 
in differential analysis, such as the random sub-key assumption
or Markov cipher assumption. The tool's primary function is to find some bounds
for the security of TWINE against differential and differential-linear cryptanalysis.
*/

#include <stdint.h>
#include <stdio.h>

typedef uint16_t u16;
typedef uint8_t u8;

#define KSIZE 80
//#define KSIZE 128


static const u8 S[16] = {0x0c, 0x00, 0x0f, 0x0a,
                         0x02, 0x0b, 0x09, 0x05,
                         0x08, 0x03, 0x0d, 0x07,
                         0x01, 0x0e, 0x06, 0x04};

static const u8 Pi[16] = {0x05, 0x00, 0x01, 0x04,
                          0x07, 0x0c, 0x03, 0x08,
                          0x0d, 0x06, 0x09, 0x02,
                          0x0f, 0x0a, 0x0b, 0x0e};

static const u8 Pi_Inv[16] = {0x01, 0x02, 0x0b, 0x06,
                              0x03, 0x00, 0x09, 0x04,
                              0x07, 0x0a, 0x0d, 0x0e,
                              0x05, 0x08, 0x0f, 0x0c};

static const u8 CON[35] = {0x01, 0x02, 0x04, 0x08,
                          0x10, 0x20, 0x03, 0x06,
                          0x0c, 0x18, 0x30, 0x23,
                          0x05, 0x0a, 0x14, 0x28,
                          0x13, 0x26, 0x0f, 0x1e,
                          0x3c, 0x3b, 0x35, 0x29,
                          0x11, 0x22, 0x07, 0x0e,
                          0x1c, 0x38, 0x33, 0x25,
                          0x09, 0x12, 0x24
                          };

void OneRound(u8 x[16], u8 k[8]);
void OneRound_Inv(u8 x[16], u8 k[8]);
void KeySch(int nrounds, const u16 *key, u8 output[36][8]);
void OneRound(u8 x[16], u8 k[8]);
void Encrypt(int nrounds, u8 x[16], u8 Subkey[36][8]);
void Decrypt(int nrounds, u8 x[16], u8 Subkey[36][8]);