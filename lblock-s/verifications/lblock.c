#include "lblock.h"

void EncryptKeySchedule(int nrounds, u8 key[10], u8 output[][4])
{
    u8 i, KeyR[4];
     
    output[0][3] = key[9];
    output[0][2] = key[8];
    output[0][1] = key[7];
    output[0][0] = key[6];
     
    for(i = 1; i < nrounds; i++)
    {
    // K <<< 29                 
    KeyR[3]=key[9];
    KeyR[2]=key[8];
    KeyR[1]=key[7];     
    KeyR[0]=key[6];     
    
    key[9]=(((key[6] & 0x07)<<5)&0xE0) ^ (((key[5]& 0xF8)>>3) & 0x1F);
    key[8]=(((key[5] & 0x07)<<5)&0xE0) ^ (((key[4]& 0xF8)>>3) & 0x1F);
    key[7]=(((key[4] & 0x07)<<5)&0xE0) ^ (((key[3]& 0xF8)>>3) & 0x1F);
    key[6]=(((key[3] & 0x07)<<5)&0xE0) ^ (((key[2]& 0xF8)>>3) & 0x1F);
    key[5]=(((key[2] & 0x07)<<5)&0xE0) ^ (((key[1]& 0xF8)>>3) & 0x1F);
    key[4]=(((key[1] & 0x07)<<5)&0xE0) ^ (((key[0]& 0xF8)>>3) & 0x1F);
    key[3]=(((key[0] & 0x07)<<5)&0xE0) ^ (((KeyR[3]& 0xF8)>>3) & 0x1F);
    key[2]=(((KeyR[3] & 0x07)<<5)&0xE0) ^ (((KeyR[2]& 0xF8)>>3) & 0x1F);
    key[1]=(((KeyR[2] & 0x07)<<5)&0xE0) ^ (((KeyR[1]& 0xF8)>>3) & 0x1F);
    key[0]=(((KeyR[1] & 0x07)<<5)&0xE0) ^ (((KeyR[0]& 0xF8)>>3) & 0x1F);         
                    
    // reste du keyschedule                 
    key[9]=(S9[((key[9]>>4) & 0x0F)]<<4) ^ S8[(key[9]& 0x0F)];
     
    key[6]=key[6] ^ ((i>>2) & 0x07);
    key[5]=key[5] ^ ((i & 0x03)<<6);
        
    output[i][3] = key[9];
    output[i][2] = key[8];
    output[i][1] = key[7];
    output[i][0] = key[6];                      
    }                          
}

void Swap(u8 block[8])
{
    u8 tmp[4];

    tmp[0] = block[0];
    tmp[1] = block[1];
    tmp[2] = block[2];
    tmp[3] = block[3];

    block[0] = block[4];
    block[1] = block[5];
    block[2] = block[6];
    block[3] = block[7];

    block[4] = tmp[0];
    block[5] = tmp[1];
    block[6] = tmp[2];
    block[7] = tmp[3];
}

void OneRound(u8 x[8], u8 k[4])
{
	u8 t[4],tmp[4];

	// AJOUT CLE
    tmp[0]=x[4]^k[0];
    tmp[1]=x[5]^k[1];
    tmp[2]=x[6]^k[2];
    tmp[3]=x[7]^k[3];         

    // PASSAGE DANS LES BOITES S
    tmp[0] = ((S1[((tmp[0])>>4) & 0x0F])<<4)^S0[(tmp[0] & 0x0F)];
    tmp[1] = ((S3[((tmp[1])>>4) & 0x0F])<<4)^S2[(tmp[1] & 0x0F)];
    tmp[2] = ((S5[((tmp[2])>>4) & 0x0F])<<4)^S4[(tmp[2] & 0x0F)];
    tmp[3] = ((S7[((tmp[3])>>4) & 0x0F])<<4)^S6[(tmp[3] & 0x0F)];          
    
    // PASSAGE DE LA PERMUTATION P
	t[0] =((tmp[0]>>4) & 0x0F)^(tmp[1] & 0xF0);
	t[1] = (tmp[0] & 0x0F) ^ ((tmp[1]& 0x0F)<<4);
	t[2] = ((tmp[2]>>4) & 0x0F)^(tmp[3] & 0xF0);
	t[3] = (tmp[2] & 0x0F) ^ ((tmp[3]& 0x0F)<<4);
    // FIN DE LA FONCTION F

    // PARTIE GAUCHE AVEC DECALAGE DE 8 SUR LA GAUCHE  
    tmp[0]=x[3]^t[0]; 
    tmp[1]=x[0]^t[1]; 
    tmp[2]=x[1]^t[2]; 
    tmp[3]=x[2]^t[3]; 
    
	// PARTIE DROITE
    x[0]=tmp[0];
    x[1]=tmp[1];
    x[2]=tmp[2];
    x[3]=tmp[3];      
}

void Encrypt(int nrounds, u8 x[8], u8 subkey[][4])
{
    for(int i = 0; i < nrounds; i++)
    {
       OneRound(x, subkey[i]);        
       Swap(x);
    }    
}

void OneRound_Inv(u8 y[8], u8 k[4])
{
    u8 t[4],tmp[4];

	// FAIRE PASSER Y_0, Y_1, Y_2, Y_3 dans F
	// AJOUT CLE
    tmp[0]=y[4]^k[0]; 
    tmp[1]=y[5]^k[1]; 
    tmp[2]=y[6]^k[2]; 
    tmp[3]=y[7]^k[3];  
     
 
    // PASSAGE DANS LES BOITES S
    tmp[0] = ((S1[((tmp[0])>>4) & 0x0F])<<4)^S0[(tmp[0] & 0x0F)];
    tmp[1] = ((S3[((tmp[1])>>4) & 0x0F])<<4)^S2[(tmp[1] & 0x0F)];
    tmp[2] = ((S5[((tmp[2])>>4) & 0x0F])<<4)^S4[(tmp[2] & 0x0F)];
    tmp[3] = ((S7[((tmp[3])>>4) & 0x0F])<<4)^S6[(tmp[3] & 0x0F)];    
 
    // PASSAGE DE LA PERMUTATION P
	t[0] =((tmp[0]>>4) & 0x0F)^(tmp[1] & 0xF0);
	t[1] = (tmp[0] & 0x0F) ^ ((tmp[1]& 0x0F)<<4);
	t[2] = ((tmp[2]>>4) & 0x0F)^(tmp[3] & 0xF0);
	t[3] = (tmp[2] & 0x0F) ^ ((tmp[3]& 0x0F)<<4);
    // FIN DE LA FONCTION F
    
    // PARTIE DROITE AVEC DECALAGE DE 8 SUR LA DROITE
	tmp[0]= y[0]^t[0]; 
    tmp[1]= y[1]^t[1]; 
    tmp[2]= y[2]^t[2]; 
    tmp[3]= y[3]^t[3]; 
 
 	// PARTIE GAUCHE
    y[0]=tmp[1];
    y[1]=tmp[2];
    y[2]=tmp[3];
    y[3]=tmp[0];
 
}

void Decrypt(int nrounds, u8 x[8], u8 subkey[][4])
{
     for(int i = nrounds - 1; i >= 0; i--)
     {
        Swap(x);
        OneRound_Inv(x, subkey[i]);   
     }
}
