/*********************************************
 * Reference implementation by WARP Team     *
**********************************************/
/*
Modified version of the refernce implementation
to simulation boomerang attack on WARP
Modifier: 
 */

#include "warp.h"

void printState(int *state)
{
    printf("L: ");
    for (int x = 0; x < BR_HALF; x++)
    {
        printf("%x ", state[2 * x + 0]);
    }
    printf("R: ");
    for (int x = 0; x < BR_HALF; x++)
    {
        printf("%x ", state[2 * x + 1]);
    }
    printf("\n");
}

void sboxkey(int *state, int *k, int r)
{
    for (int i = 0; i < BR_HALF; i++)
    {
        state[i] = Sbox[state[i]] ^ k[(r % 2) * 16 + i];
    }
}

void permutation(int *state)
{
    int temp[BR];
    for (int j = 0; j < BR; j++)
    {
        temp[j] = state[j];
    }
    for (int j = 0; j < BR; j++)
    {
        state[perm[j]] = temp[j];
    }
}

void inv_permutation(int *state)
{
    int temp[BR];
    for (int j = 0; j < BR; j++)
    {
        temp[j] = state[j];
    }
    for (int j = 0; j < BR; j++)
    {
        state[j] = temp[perm[j]];
    }
}

void enc(int *m, int *c, int *k, int R)
{
    /*intermediate value*/
    int state[BR];

    /*left half intermediate value*/
    int temp[BR_HALF];

    for (int i = 0; i < BR; i++)
    {
        state[i] = m[i];
    }

    /*round function(1 to 40 round)*/
    for (int i = 0; i < R; i++)
    {
        #if PRINT_INTER
        printf("%d round\n", i + 1);
        printState(state);
        #endif

        for (int j = 0; j < BR_HALF; j++)
        {
            temp[j] = state[j * 2];
        }
        /*insert key and Sbox*/
        sboxkey(temp, k, i);
        /*XOR*/
        for (int j = 0; j < BR_HALF; j++)
        {
            state[2 * j + 1] = state[2 * j + 1] ^ temp[j];
        }
        /*add round constants*/
        state[1] = state[1] ^ RC0[i];
        state[3] = state[3] ^ RC1[i];

        /*permutation*/
        permutation(state);
    }

    /*last round function */
    #if PRINT_INTER
    printf("%d round\n", R);
    printState(state);
    #endif

    #if PRINT_INTER
    printState(state);
    #endif

    /*no permutation in the last round*/

    /*copy ciphertext*/
    for (int i = 0; i < BR; i++)
    {
        c[i] = state[i];
    }

}

void dec(int *m, int *c, int *k, int R)
{
    /*intermediate value*/
    int state[BR];

    /*left half intermediate value*/
    int temp[BR_HALF];

    for (int i = 0; i < BR; i++)
    {
        state[i] = c[i];
    }

    /*round function(1 to 40 round)*/
    for (int i = R - 1; i >= 0; i--)
    {
        #if PRINT_INTER
        printf("%d round\n", i + 1);
        printState(state);
        #endif

        /*inverse of permutation*/
        inv_permutation(state);
        /*add round constants*/
        state[1] = state[1] ^ RC0[i];
        state[3] = state[3] ^ RC1[i];

        for (int j = 0; j < BR_HALF; j++)
        {
            temp[j] = state[j * 2];
        }
        /*insert key and Sbox*/
        sboxkey(temp, k, i);
        /*XOR*/
        for (int j = 0; j < BR_HALF; j++)
        {
            state[2 * j + 1] = state[2 * j + 1] ^ temp[j];
        }
    }

    /*last round function */
    #if PRINT_INTER
    printf("%d round\n", R);
    printState(state);
    #endif

    #if PRINT_INTER
    printState(state);
    #endif

    /*no permutation in the last round*/

    /*copy plaintext*/
    for (int i = 0; i < BR; i++)
    {
        m[i] = state[i];
    }

}