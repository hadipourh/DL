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

#define PRINT

#include "twine.h"

void KeySch(int nrounds, const u16 *key, u8 output[][8])
{
    int i;
    u8 KeyR[KSIZE/4],temp, temp1, temp2, temp3;
          
    for(i=0;i<(KSIZE/4);i++)
    {
      KeyR[i] = (key[(i/4)]>>(4*(i&0x03))) & 0x0F;
    }

    for(i=0;i < nrounds + 1;i++)
    {
#if KSIZE == 80                      
      output[i][0] = KeyR[1];
      output[i][1] = KeyR[3];
      output[i][2] = KeyR[4];
      output[i][3] = KeyR[6];
      output[i][4] = KeyR[13];
      output[i][5] = KeyR[14];
      output[i][6] = KeyR[15];
      output[i][7] = KeyR[16];        
      
      KeyR[1] = KeyR[1] ^ S[KeyR[0]];
      KeyR[4] = KeyR[4] ^ S[KeyR[16]];
      KeyR[7] = KeyR[7] ^ (CON[i]>>3);
      KeyR[19] = KeyR[19] ^ (CON[i]&0x07);
      
      temp = KeyR[0];
      KeyR[0] = KeyR[1];
      KeyR[1] = KeyR[2];      
      KeyR[2] = KeyR[3];      
      KeyR[3] = temp;            
      
      temp = KeyR[0];
      temp1 = KeyR[1];
      temp2 = KeyR[2];
      temp3 = KeyR[3];
      
      KeyR[0] = KeyR[4];
      KeyR[1] = KeyR[5];      
      KeyR[2] = KeyR[6];      
      KeyR[3] = KeyR[7];            
      
      KeyR[4] = KeyR[8];
      KeyR[5] = KeyR[9];      
      KeyR[6] = KeyR[10];      
      KeyR[7] = KeyR[11];           
      
      KeyR[8] = KeyR[12];
      KeyR[9] = KeyR[13];      
      KeyR[10] = KeyR[14];  
      KeyR[11] = KeyR[15];

      KeyR[12] = KeyR[16];
      KeyR[13] = KeyR[17];      
      KeyR[14] = KeyR[18];      
      KeyR[15] = KeyR[19];
      
      KeyR[16] = temp;
      KeyR[17] = temp1;      
      KeyR[18] = temp2;      
      KeyR[19] = temp3;     
      
#elif KSIZE == 128
output[i][0] = KeyR[2];
      output[i][1] = KeyR[3];
      output[i][2] = KeyR[12];
      output[i][3] = KeyR[15];
      output[i][4] = KeyR[17];
      output[i][5] = KeyR[18];
      output[i][6] = KeyR[28];
      output[i][7] = KeyR[31];        
      
      KeyR[1]=KeyR[1] ^S[KeyR[0]];
      KeyR[4]=KeyR[4] ^S[KeyR[16]];
      KeyR[23]=KeyR[23] ^S[KeyR[30]];
      
      KeyR[7]=KeyR[7] ^(CON[i]>>3);
      KeyR[19]=KeyR[19] ^(CON[i]&0x07);
      
      temp=KeyR[0];
      KeyR[0]=KeyR[1];
      KeyR[1]=KeyR[2];      
      KeyR[2]=KeyR[3];      
      KeyR[3]=temp;  
                
      temp=KeyR[0];
      temp1=KeyR[1];
      temp2=KeyR[2];
      temp3=KeyR[3];
      
      KeyR[0]=KeyR[4];
      KeyR[1]=KeyR[5];      
      KeyR[2]=KeyR[6];      
      KeyR[3]=KeyR[7];            
      
      KeyR[4]=KeyR[8];
      KeyR[5]=KeyR[9];      
      KeyR[6]=KeyR[10];      
      KeyR[7]=KeyR[11];           
      
      KeyR[8]=KeyR[12];
      KeyR[9]=KeyR[13];      
      KeyR[10]=KeyR[14];      
      KeyR[11]=KeyR[15];     

      KeyR[12]=KeyR[16];
      KeyR[13]=KeyR[17];      
      KeyR[14]=KeyR[18];      
      KeyR[15]=KeyR[19];

      KeyR[16]=KeyR[20];
      KeyR[17]=KeyR[21];      
      KeyR[18]=KeyR[22];      
      KeyR[19]=KeyR[23];
      
      KeyR[20]=KeyR[24];
      KeyR[21]=KeyR[25];      
      KeyR[22]=KeyR[26];      
      KeyR[23]=KeyR[27];

      KeyR[24]=KeyR[28];
      KeyR[25]=KeyR[29];      
      KeyR[26]=KeyR[30];      
      KeyR[27]=KeyR[31];

      KeyR[28]=temp;
      KeyR[29]=temp1;      
      KeyR[30]=temp2;      
      KeyR[31]=temp3;
#endif             
  }
}

void OneRound(u8 x[16], u8 k[8])
{
  u8 t[16];
  u8 i;

  for(i = 0; i < 8; i++)
  {
  x[2*i+1] = (S[x[2*i]^k[i]]^x[2*i+1]) & 0x0F;
  }

  for(i=0; i < 16; i++)
  {
    t[Pi[i]]=x[i];
  }

  for(i = 0; i < 16; i++)
  {
  x[i] = t[i];
  }
}

void Encrypt(int nrounds, u8 x[16], u8 Subkey[36][8])
{
  u8 i;
  // apply round function nrounds times
  for(i = 0; i < nrounds; i++)
  {
    OneRound(x,Subkey[i]);
  }
}

void OneRound_Inv(u8 x[16], u8 k[8])
{
  u8 t[16];
  u8 i;

  for(i=0;i<16;i++)
  {
    t[i] = x[Pi[i]];
  }

  for(i=0;i<8;i++)
  {    
    x[2*i+1] = (S[t[2*i]^k[i]]^t[2*i+1]) & 0x0F;
    x[2*i] = t[2*i];
  }
}


void Decrypt(int nrounds, u8 x[16], u8 Subkey[36][8])
{  
  // apply the inverse of round function nrounds times
  for(int i = (nrounds - 1); i >= 0; i--)
  {
    OneRound_Inv(x, Subkey[i]);
  }
}
