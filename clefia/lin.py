#!/usr/env/bin python3
#-*- coding: UTF-8 -*-

"""
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

Disclaimer: We acknowledge that the CLEFIA block cipher doesn't adhere to statistical assumptions 
in linear analysis, such as the random sub-key assumption
or Markov cipher assumption. The tool's primary function is to find some bounds
for the security of CLEFIA against linear and differential-linear cryptanalysis.
"""

from argparse import ArgumentParser, RawTextHelpFormatter
import yaml
import time
from gurobipy import *
import math
import uuid
import os

class Lin:
    """
    This class is used to find linear trail as well as
    computing the linear effect of CLEFIA block cipher.

    x_roundNumber_branchNumber_byteNumber_bitNumber

    Variable mapping in differential cryptanalysis:

    ... x_r_0                                                    --- x_r_1 ...
    ... |                                                        |     |
    ... |--- x_r_0 ----> | S |--- y_r_0 --->|4*4 MDS|--- z_r_0---+---->+   ...
    ... |                                                              |   ...

    Variable mapping in linear cryptanalysis:

    ... x_r_0                                                     --- x_r_1 ...
    ... |                                                         |     |
    ... |--- y_r_0 ----> | S |--- z_r_0 --->|4*4 MDS|--- x_r_1 ---+---->+   ...
    ... |                                                               |   ...

    """

    diff_count = 0

    def __init__(self, params) -> None:
        Lin.diff_count += 1

        self.nrounds = params["nrounds"]
        self.time_limit = params["timelimit"]
        self.start_weight = params["startweight"]
        self.end_weight = params["endweight"]
        self.fixed_variables = params['fixedVariables']
        self.mode = params['mode']
        self.number_of_trails = params["numberoftrails"]
        self.eps = 1e-3
        self.big_m = 2*8        
        self.lp_file_name = f"clefia_nr_{self.nrounds}_{uuid.uuid4()}.lp"        
        self.result_file_name = f"clefia_nr_{self.nrounds}_{uuid.uuid4()}.result"
        self.binary_variables = []
        self.integer_variables = []

        """
        We used S-box Analyzer to encode the linear behavior of S-boxes:
        https://github.com//sboxanalyzer
        Encoding linear behavior of S-boxes
        a0, ..., a7 (a0: msb of input difference)
        b0, ..., b7 (b0: msb of output difference)
        p0, p1, p2, p3, p4 encode the squared correlation of each possible linear transition:

        For S0:
        p0 <--> -10
        p1 <--> -8
        p2 <--> -6.830075
        p3 <--> -6
        p4 <--> -5.356144
        p5 <--> -4.830075
        p6 <--> -4.385290

        For S1:
        p0 <--> -12
        p1 <--> -10
        p2 <--> -8.830075
        p3 <--> -8
        p4 <--> -7.356144
        p5 <--> -6.830075
        p6 <--> -6.385290
        p7 <--> -6
        """

        self.sbox_inequalities = [0 for _ in range(2)]
        self.sbox_inequalities[0] = dict()
        self.sbox_inequalities[1] = dict()
        self.s0_correlations = {"10": "10", "8": "8", "6_830075": "6.830075", "6": "6", "5_356144": "5.356144", "4_830075": "4.830075", "4_385290": "4.385290"}
        self.s1_correlations = {"12": "12", "10": "10", "8_830075": "8.830075", "8": "8", "7_356144": "7.356144", "6_830075": "6.830075", "6_385290": "6.385290", "6": "6"}
        for pr in self.s0_correlations.keys():
            with open(os.path.join("lat-encoding", f"s0_{pr}.txt")) as fileobj:
                self.sbox_inequalities[0][pr] = fileobj.readlines()
        for pr in self.s1_correlations.keys():
            with open(os.path.join("lat-encoding", f"s1_{pr}.txt")) as fileobj:
                self.sbox_inequalities[1][pr] = fileobj.readlines()

        self.mds0_constraints_template = [
            'a_00000 + a_01011 + a_01100 + a_01101 + a_01111 + a_10010 + a_10011 + a_10100 + a_10110 + a_11010 + a_11101 + a_11110 + a_11111 + b_00000 - 2 u_00000 = 0',
            'u_00000 <= 7',
            'u_00000 >= 0',
            'a_00001 + a_01000 + a_10011 + a_10100 + a_10101 + a_10111 + a_11000 + a_11011 + a_11100 + a_11101 + a_11111 + b_00001 - 2 u_00001 = 0',
            'u_00001 <= 6',
            'u_00001 >= 0',
            'a_00010 + a_01001 + a_10000 + a_11000 + a_11001 + b_00010 - 2 u_00010 = 0',
            'u_00010 <= 3',
            'u_00010 >= 0',
            'a_00011 + a_01010 + a_10001 + a_11001 + a_11010 + b_00011 - 2 u_00011 = 0',
            'u_00011 <= 3',
            'u_00011 >= 0',
            'a_00100 + a_01011 + a_10010 + a_11010 + a_11011 + b_00100 - 2 u_00100 = 0',
            'u_00100 <= 3',
            'u_00100 >= 0',
            'a_00101 + a_01100 + a_10011 + a_11011 + a_11100 + b_00101 - 2 u_00101 = 0',
            'u_00101 <= 3',
            'u_00101 >= 0',
            'a_00110 + a_01101 + a_10100 + a_11100 + a_11101 + b_00110 - 2 u_00110 = 0',
            'u_00110 <= 3',
            'u_00110 >= 0',
            'a_00111 + a_01110 + a_10101 + a_11101 + a_11110 + b_00111 - 2 u_00111 = 0',
            'u_00111 <= 3',
            'u_00111 >= 0',
            'a_00011 + a_00100 + a_00101 + a_00111 + a_01000 + a_10010 + a_10101 + a_10110 + a_10111 + a_11010 + a_11011 + a_11100 + a_11110 + b_01000 - 2 u_01000 = 0',
            'u_01000 <= 7',
            'u_01000 >= 0',
            'a_00000 + a_01001 + a_10000 + a_10011 + a_10100 + a_10101 + a_10111 + a_11011 + a_11100 + a_11101 + a_11111 + b_01001 - 2 u_01001 = 0',
            'u_01001 <= 6',
            'u_01001 >= 0',
            'a_00001 + a_01010 + a_10000 + a_10001 + a_11000 + b_01010 - 2 u_01010 = 0',
            'u_01010 <= 3',
            'u_01010 >= 0',
            'a_00010 + a_01011 + a_10001 + a_10010 + a_11001 + b_01011 - 2 u_01011 = 0',
            'u_01011 <= 3',
            'u_01011 >= 0',
            'a_00011 + a_01100 + a_10010 + a_10011 + a_11010 + b_01100 - 2 u_01100 = 0',
            'u_01100 <= 3',
            'u_01100 >= 0',
            'a_00100 + a_01101 + a_10011 + a_10100 + a_11011 + b_01101 - 2 u_01101 = 0',
            'u_01101 <= 3',
            'u_01101 >= 0',
            'a_00101 + a_01110 + a_10100 + a_10101 + a_11100 + b_01110 - 2 u_01110 = 0',
            'u_01110 <= 3',
            'u_01110 >= 0',
            'a_00110 + a_01111 + a_10101 + a_10110 + a_11101 + b_01111 - 2 u_01111 = 0',
            'u_01111 <= 3',
            'u_01111 >= 0',
            'a_00010 + a_00011 + a_00100 + a_00110 + a_01010 + a_01101 + a_01110 + a_01111 + a_10000 + a_11011 + a_11100 + a_11101 + a_11111 + b_10000 - 2 u_10000 = 0',
            'u_10000 <= 7',
            'u_10000 >= 0',
            'a_00011 + a_00100 + a_00101 + a_00111 + a_01000 + a_01011 + a_01100 + a_01101 + a_01111 + a_10001 + a_11000 + b_10001 - 2 u_10001 = 0',
            'u_10001 <= 6',
            'u_10001 >= 0',
            'a_00000 + a_01000 + a_01001 + a_10010 + a_11001 + b_10010 - 2 u_10010 = 0',
            'u_10010 <= 3',
            'u_10010 >= 0',
            'a_00001 + a_01001 + a_01010 + a_10011 + a_11010 + b_10011 - 2 u_10011 = 0',
            'u_10011 <= 3',
            'u_10011 >= 0',
            'a_00010 + a_01010 + a_01011 + a_10100 + a_11011 + b_10100 - 2 u_10100 = 0',
            'u_10100 <= 3',
            'u_10100 >= 0',
            'a_00011 + a_01011 + a_01100 + a_10101 + a_11100 + b_10101 - 2 u_10101 = 0',
            'u_10101 <= 3',
            'u_10101 >= 0',
            'a_00100 + a_01100 + a_01101 + a_10110 + a_11101 + b_10110 - 2 u_10110 = 0',
            'u_10110 <= 3',
            'u_10110 >= 0',
            'a_00101 + a_01101 + a_01110 + a_10111 + a_11110 + b_10111 - 2 u_10111 = 0',
            'u_10111 <= 3',
            'u_10111 >= 0',
            'a_00010 + a_00101 + a_00110 + a_00111 + a_01010 + a_01011 + a_01100 + a_01110 + a_10011 + a_10100 + a_10101 + a_10111 + a_11000 + b_11000 - 2 u_11000 = 0',
            'u_11000 <= 7',
            'u_11000 >= 0',
            'a_00000 + a_00011 + a_00100 + a_00101 + a_00111 + a_01011 + a_01100 + a_01101 + a_01111 + a_10000 + a_11001 + b_11001 - 2 u_11001 = 0',
            'u_11001 <= 6',
            'u_11001 >= 0',
            'a_00000 + a_00001 + a_01000 + a_10001 + a_11010 + b_11010 - 2 u_11010 = 0',
            'u_11010 <= 3',
            'u_11010 >= 0',
            'a_00001 + a_00010 + a_01001 + a_10010 + a_11011 + b_11011 - 2 u_11011 = 0',
            'u_11011 <= 3',
            'u_11011 >= 0',
            'a_00010 + a_00011 + a_01010 + a_10011 + a_11100 + b_11100 - 2 u_11100 = 0',
            'u_11100 <= 3',
            'u_11100 >= 0',
            'a_00011 + a_00100 + a_01011 + a_10100 + a_11101 + b_11101 - 2 u_11101 = 0',
            'u_11101 <= 3',
            'u_11101 >= 0',
            'a_00100 + a_00101 + a_01100 + a_10101 + a_11110 + b_11110 - 2 u_11110 = 0',
            'u_11110 <= 3',
            'u_11110 >= 0',
            'a_00101 + a_00110 + a_01101 + a_10110 + a_11111 + b_11111 - 2 u_11111 = 0',
            'u_11111 <= 3',
            'u_11111 >= 0']

        self.mds1_constraints_template = [
            'a_00000 + a_01001 + a_01010 + a_01011 + a_01101 + a_10011 + a_10100 + a_10101 + a_10111 + a_11001 + a_11010 + a_11100 + a_11111 + b_00000 - 2 u_00000 = 0',
            'u_00000 <= 7',
            'u_00000 >= 0',
            'a_00001 + a_01010 + a_01011 + a_01100 + a_01110 + a_10000 + a_11000 + a_11010 + a_11011 + a_11100 + a_11110 + b_00001 - 2 u_00001 = 0',
            'u_00001 <= 6',
            'u_00001 >= 0',
            'a_00010 + a_01011 + a_01100 + a_01101 + a_01111 + a_10001 + a_11001 + a_11011 + a_11100 + a_11101 + a_11111 + b_00010 - 2 u_00010 = 0',
            'u_00010 <= 6',
            'u_00010 >= 0',
            'a_00011 + a_01000 + a_10010 + a_11000 + a_11010 + b_00011 - 2 u_00011 = 0',
            'u_00011 <= 3',
            'u_00011 >= 0',
            'a_00100 + a_01001 + a_10011 + a_11001 + a_11011 + b_00100 - 2 u_00100 = 0',
            'u_00100 <= 3',
            'u_00100 >= 0',
            'a_00101 + a_01010 + a_10100 + a_11010 + a_11100 + b_00101 - 2 u_00101 = 0',
            'u_00101 <= 3',
            'u_00101 >= 0',
            'a_00110 + a_01011 + a_10101 + a_11011 + a_11101 + b_00110 - 2 u_00110 = 0',
            'u_00110 <= 3',
            'u_00110 >= 0',
            'a_00111 + a_01100 + a_10110 + a_11100 + a_11110 + b_00111 - 2 u_00111 = 0',
            'u_00111 <= 3',
            'u_00111 >= 0',
            'a_00001 + a_00010 + a_00011 + a_00101 + a_01000 + a_10001 + a_10010 + a_10100 + a_10111 + a_11011 + a_11100 + a_11101 + a_11111 + b_01000 - 2 u_01000 = 0',
            'u_01000 <= 7',
            'u_01000 >= 0',
            'a_00010 + a_00011 + a_00100 + a_00110 + a_01001 + a_10000 + a_10010 + a_10011 + a_10100 + a_10110 + a_11000 + b_01001 - 2 u_01001 = 0',
            'u_01001 <= 6',
            'u_01001 >= 0',
            'a_00011 + a_00100 + a_00101 + a_00111 + a_01010 + a_10001 + a_10011 + a_10100 + a_10101 + a_10111 + a_11001 + b_01010 - 2 u_01010 = 0',
            'u_01010 <= 6',
            'u_01010 >= 0',
            'a_00000 + a_01011 + a_10000 + a_10010 + a_11010 + b_01011 - 2 u_01011 = 0',
            'u_01011 <= 3',
            'u_01011 >= 0',
            'a_00001 + a_01100 + a_10001 + a_10011 + a_11011 + b_01100 - 2 u_01100 = 0',
            'u_01100 <= 3',
            'u_01100 >= 0',
            'a_00010 + a_01101 + a_10010 + a_10100 + a_11100 + b_01101 - 2 u_01101 = 0',
            'u_01101 <= 3',
            'u_01101 >= 0',
            'a_00011 + a_01110 + a_10011 + a_10101 + a_11101 + b_01110 - 2 u_01110 = 0',
            'u_01110 <= 3',
            'u_01110 >= 0',
            'a_00100 + a_01111 + a_10100 + a_10110 + a_11110 + b_01111 - 2 u_01111 = 0',
            'u_01111 <= 3',
            'u_01111 >= 0',
            'a_00011 + a_00100 + a_00101 + a_00111 + a_01001 + a_01010 + a_01100 + a_01111 + a_10000 + a_11001 + a_11010 + a_11011 + a_11101 + b_10000 - 2 u_10000 = 0',
            'u_10000 <= 7',
            'u_10000 >= 0',
            'a_00000 + a_01000 + a_01010 + a_01011 + a_01100 + a_01110 + a_10001 + a_11010 + a_11011 + a_11100 + a_11110 + b_10001 - 2 u_10001 = 0',
            'u_10001 <= 6',
            'u_10001 >= 0',
            'a_00001 + a_01001 + a_01011 + a_01100 + a_01101 + a_01111 + a_10010 + a_11011 + a_11100 + a_11101 + a_11111 + b_10010 - 2 u_10010 = 0',
            'u_10010 <= 6',
            'u_10010 >= 0',
            'a_00010 + a_01000 + a_01010 + a_10011 + a_11000 + b_10011 - 2 u_10011 = 0',
            'u_10011 <= 3',
            'u_10011 >= 0',
            'a_00011 + a_01001 + a_01011 + a_10100 + a_11001 + b_10100 - 2 u_10100 = 0',
            'u_10100 <= 3',
            'u_10100 >= 0',
            'a_00100 + a_01010 + a_01100 + a_10101 + a_11010 + b_10101 - 2 u_10101 = 0',
            'u_10101 <= 3',
            'u_10101 >= 0',
            'a_00101 + a_01011 + a_01101 + a_10110 + a_11011 + b_10110 - 2 u_10110 = 0',
            'u_10110 <= 3',
            'u_10110 >= 0',
            'a_00110 + a_01100 + a_01110 + a_10111 + a_11100 + b_10111 - 2 u_10111 = 0',
            'u_10111 <= 3',
            'u_10111 >= 0',
            'a_00001 + a_00010 + a_00100 + a_00111 + a_01011 + a_01100 + a_01101 + a_01111 + a_10001 + a_10010 + a_10011 + a_10101 + a_11000 + b_11000 - 2 u_11000 = 0',
            'u_11000 <= 7',
            'u_11000 >= 0',
            'a_00000 + a_00010 + a_00011 + a_00100 + a_00110 + a_01000 + a_10010 + a_10011 + a_10100 + a_10110 + a_11001 + b_11001 - 2 u_11001 = 0',
            'u_11001 <= 6',
            'u_11001 >= 0',
            'a_00001 + a_00011 + a_00100 + a_00101 + a_00111 + a_01001 + a_10011 + a_10100 + a_10101 + a_10111 + a_11010 + b_11010 - 2 u_11010 = 0',
            'u_11010 <= 6',
            'u_11010 >= 0',
            'a_00000 + a_00010 + a_01010 + a_10000 + a_11011 + b_11011 - 2 u_11011 = 0',
            'u_11011 <= 3',
            'u_11011 >= 0',
            'a_00001 + a_00011 + a_01011 + a_10001 + a_11100 + b_11100 - 2 u_11100 = 0',
            'u_11100 <= 3',
            'u_11100 >= 0',
            'a_00010 + a_00100 + a_01100 + a_10010 + a_11101 + b_11101 - 2 u_11101 = 0',
            'u_11101 <= 3',
            'u_11101 >= 0',
            'a_00011 + a_00101 + a_01101 + a_10011 + a_11110 + b_11110 - 2 u_11110 = 0',
            'u_11110 <= 3',
            'u_11110 >= 0',
            'a_00100 + a_00110 + a_01110 + a_10100 + a_11111 + b_11111 - 2 u_11111 = 0',
            'u_11111 <= 3',
            'u_11111 >= 0']            

    @staticmethod
    def ordered_set(seq):
        """
        This method eliminates duplicated elements in a given list,
        and returns a list in which each elements appears only once
        """

        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    @staticmethod
    def flatten_state(s):
        """
        Convert the multidimensional stat array to a flat array
        """

        state_bits = []
        for i in range(len(s)):
            for j in range(len(s[i])):
                for k in range(len(s[i][j])):
                    state_bits += [s[i][j][k]]
        return state_bits

    @staticmethod
    def convert_str_to_binarystatevector(str_hex):
        assert(len(str_hex) == 32)
        state = [0]*128
        for byte in range(16):
            for i in range(7, -1, -1):
                state[i] = ((int(str_hex[2*byte:2*(byte + 1)], base=16) >> i) & 0x1)
        return state

    def apply_swap(self, state):
        """
        Swap two halves of the state array
        """

        temp = [0 for _ in range(4)]
        for i in range(4):
            temp[i] = state[(i + 1) % 4]
        return temp

    def apply_inv_swap(self, state):
        """
        Apply the inverse of swap operation
        """

        temp = [0 for _ in range(4)]
        for i in range(4):
            temp[i] = state[(i + 3) % 4]
        return temp


    def generate_round_x_variables(self, rn):
        """
        Generate the input variables of rn'th round

        :param rn int: round number
        """

        x = [[[f"x_{rn}_{bn}_{byten}_{bitn}" for bitn in range(8)] for byten in range(4)] for bn in range(4)]
        self.binary_variables.extend(self.flatten_state(x))
        return x

    def generate_round_yz_variables(self, rn):
        """
        Generate the variables representing the activeness at the output of MDS matrix

        :param rn int: round number
        """

        y = [[[f"y_{rn}_{bn}_{byten}_{bitn}" for bitn in range(8)] for byten in range(4)] for bn in range(2)]
        z = [[[f"z_{rn}_{bn}_{byten}_{bitn}" for bitn in range(8)] for byten in range(4)] for bn in range(2)]
        self.binary_variables.extend(self.flatten_state(z))
        return y, z

    def generate_round_pr_variables(self, rn):
        """
        Generate the variables encoding the probability of S-boxes
        """
        pr = [[0 for byten in range(4)] for bn in range(2)]
        for bn in range(2):
            for byten in range(4):
                if (bn == 0 and byten % 2 == 0) or (bn == 1 and byten % 2 == 1):
                    pr[bn][byten] =  [f"q_{rn}_{bn}_{byten}_{qn}" for qn in self.s0_correlations.keys()]
                elif (bn == 0 and byten % 2 == 1) or (bn == 1 and byten % 2 == 0):
                    pr[bn][byten] =  [f"q_{rn}_{bn}_{byten}_{qn}" for qn in self.s1_correlations.keys()]
        self.binary_variables.extend(self.flatten_state(pr))
        return pr

    def constraints_by_equality(self, a, b):
        """
        Generate constraints for equality
        a = b
        """
        constraint = f"{a} - {b} = 0\n"
        return constraint

    def constraint_by_byte_equality(self, a, b):
        """
        Generate constraints corresponding
        to equality of two bytes
        """

        constraints = ""
        for bit in range(8):
            constraints += f"{a[bit]} - {b[bit]} = 0\n"
        return constraints

    def constraints_by_xor(self, a, b, c):
        """
        a + b = c
        model:
        - a - b - c >= -2
          a + b - c >= 0
          a - b + c >= 0
        - a + b + c >= 0
        """

        constraints = f"- {a} - {b} - {c} >= -2\n"
        constraints += f"{a} + {b} - {c} >= 0\n"
        constraints += f"{a} - {b} + {c} >= 0\n"
        constraints += f"- {a} + {b} + {c} >= 0\n"
        return constraints

    def constraints_by_byte_xor(self, a, b, c):
        """
        Generate constraints for XOR of bytes
        a + b -> c
        """

        constraints = ""
        for bit in range(8):
            constraints += self.constraints_by_xor(a[bit], b[bit], c[bit])
        return constraints

    def generate_constraints_by_sbox(self, rn, bn, byten, di, do):
        """
        Generate constraints modeling the DDT of S-box

        :return constraints encoding the DDT of S-box:
        :rtype str:
        """

        constraints = ""
        indicator_variable = f"Q_{rn}_{bn}_{byten}"
        pr_indicators = []
        self.binary_variables.append(indicator_variable)
        # Link activeness indicator to input/output
        sum_input = " + ".join(di)
        sum_output = " + ".join(do)
        constraints += f"{sum_input} - {indicator_variable} >= 0\n"
        constraints += f"{sum_output} - {indicator_variable} >= 0\n"
        for i in range(8):
            constraints += f"{indicator_variable} - {di[i]} >= 0\n"
            constraints += f"{indicator_variable} - {do[i]} >= 0\n"

        if (bn == 0 and byten % 2 == 0) or (bn == 1 and byten % 2 == 1):
            for q in self.s0_correlations.keys():
                q_indicator = f"q_{rn}_{bn}_{byten}_{q}"
                self.binary_variables.append(q_indicator)
                pr_indicators.append(q_indicator)
                assert(q_indicator in self.binary_variables)
                for ineq in self.sbox_inequalities[0][q]:
                    for i in range(8):
                        ineq = ineq.replace(f"a{i}", di[i])
                        ineq = ineq.replace(f"b{i}", do[i])
                    ineq = ineq.split(' >= ')
                    lhs = f"{ineq[0]} - {self.big_m} {q_indicator}"
                    rhs = "{}".format(int(ineq[1]) - self.big_m)
                    ineq = f"{lhs} >= {rhs}\n"
                    constraints += ineq
        if (bn == 0 and byten % 2 == 1) or (bn == 1 and byten % 2 == 0):
            for q in self.s1_correlations.keys():
                q_indicator = f"q_{rn}_{bn}_{byten}_{q}"
                self.binary_variables.append(q_indicator)
                pr_indicators.append(q_indicator)
                assert(q_indicator in self.binary_variables)
                for ineq in self.sbox_inequalities[1][q]:
                    for i in range(8):
                        ineq = ineq.replace(f"a{i}", di[i])
                        ineq = ineq.replace(f"b{i}", do[i])
                    ineq = ineq.split(' >= ')
                    lhs = f"{ineq[0]} - {self.big_m} {q_indicator}"
                    rhs = "{}".format(int(ineq[1]) - self.big_m)
                    ineq = f"{lhs} >= {rhs}\n"
                    constraints += ineq
        # Link probability indicators to the activeness indicator
        sum_pr = " + ".join(pr_indicators)
        constraints += f"{sum_pr} - {indicator_variable} = 0\n"
        return constraints

    def generate_constraints_by_mds(self, rn, bn, mi, mo):
        """
        Generate constraints encoding the MDS layers

        :param rn int: round number
        :param bn int: branch number in {0, 1}
        :param mi list: vector of input variables
        :param mo list: vector of output variables
        :rtype: string
        :return: constraints modeling MDS in string format
        """

        constraints = ""
        dummay_variables = [f"mds_{rn}_{bn}_{bitn}" for bitn in range(32)]
        self.integer_variables.extend(dummay_variables)
        mi = [mi[byten][bitn] for byten in range(4) for bitn in range(8)]
        mo = [mo[byten][bitn] for byten in range(4) for bitn in range(8)]
        if bn == 0:
            for ineq in self.mds0_constraints_template:
                for i in range(32):
                    bin_index = bin(i)[2:].zfill(5)
                    ineq = ineq.replace(f"a_{bin_index}", mi[i])
                    ineq = ineq.replace(f"b_{bin_index}", mo[i])
                    ineq = ineq.replace(f"u_{bin_index}", dummay_variables[i])
                constraints += ineq + "\n"
        elif bn == 1:
            for ineq in self.mds1_constraints_template:
                for i in range(32):
                    bin_index = bin(i)[2:].zfill(5)
                    ineq = ineq.replace(f"a_{bin_index}", mi[i])
                    ineq = ineq.replace(f"b_{bin_index}", mo[i])
                    ineq = ineq.replace(f"u_{bin_index}", dummay_variables[i])
                constraints += ineq + "\n"
        return constraints

    def generate_objective_function(self):
        """
        Generate the objective function of MILP model
        The objective is minimizing the summation of variables
        which encode the weight (or probability exponent) the
        linear trail
        """

        objective_function = "minimize\n"
        objective_terms = []
        for rn in range(self.nrounds):
            for bn in range(2):
                for byten in range(4):
                    if (bn == 0 and byten % 2 == 0) or (bn == 1 and byten % 2 == 1):
                        for q in self.s0_correlations.keys():
                            indicator = f"q_{rn}_{bn}_{byten}_{q}"
                            objective_terms += [f"{self.s0_correlations[q]} {indicator}"]
                    elif (bn == 0 and byten % 2 == 1) or (bn == 1 and byten % 2 == 0):
                        for q in self.s1_correlations.keys():
                            indicator = f"q_{rn}_{bn}_{byten}_{q}"
                            objective_terms += [f"{self.s1_correlations[q]} {indicator}"]
        objective_function += " + ".join(objective_terms) + "\n"
        return objective_function

    def generate_constraints(self):
        """
        Generate the constraints describing the propagation
        of linear trails through a reduced-round CLEFIA
        """

        constraints = "subject to\n"
        for rn in range(self.nrounds):
            x_in = self.generate_round_x_variables(rn)
            x_out = self.generate_round_x_variables(rn + 1)
            x_mid = self.apply_inv_swap(x_out)
            y, z = self.generate_round_yz_variables(rn)
            for bn in range(2):
                for byten in range(4):                    
                    constraints += self.generate_constraints_by_sbox(rn, bn, byten, y[bn][byten], z[bn][byten])                    
                    constraints += self.constraint_by_byte_equality(x_in[2*bn + 1][byten], x_mid[2*bn + 1][byten])                    
                    constraints += self.constraints_by_byte_xor(x_in[2*bn][byten], y[bn][byten], x_mid[2*bn][byten])                
                constraints += self.generate_constraints_by_mds(rn, bn, z[bn], x_in[2*bn + 1])
        return constraints

    def declare_binary_vars(self):
        """
        Declare binary variables of MILP model
        """

        self.binary_variables = self.ordered_set(self.binary_variables)
        constraints = "Binary\n"
        constraints += "\n".join(self.binary_variables) + "\n"
        constraints += "General\n"
        constraints += "\n".join(self.integer_variables) + "\n"
        return constraints

    def exclude_trivial_trail(self):
        """
        Exclude all-zero solution from the solution space
        """

        input_diff = self.flatten_state(self.generate_round_x_variables(0))
        input_diff = " + ".join(input_diff)
        constraint = f"{input_diff} >= 1\n"
        return constraint

    def declare_fixed_variables(self):
        lp_contents = ""
        for cond in self.fixed_variables.items():
            var = cond[0]
            val = cond[1]
            var = var.split('_')
            if len(var) == 2:
                assert(var[0] == "x")
                state_vars = self.generate_round_x_variables(var[1])
                state_vars = self.flatten_state(state_vars)
                state_values = list(bin(int(val, 16))[2:].zfill(128))
                for i in range(128):
                    lp_contents += f"{state_vars[i]} = {state_values[i]}\n"
            elif len(var) == 4:
                assert(var[0] == "x")
                state_vars = [f"x_{var[1]}_{var[2]}_{var[3]}_{bit}" for bit in range(8)]
                state_values = list(bin(int(val, 16))[2:].zfill(8))
                for i in range(8):
                    lp_contents += f"{state_vars[i]} = {state_values[i]}\n"
            elif len(var) == 5:
                assert(var[0] == "x")
                state_vars = f"x_{var[1]}_{var[2]}_{var[3]}_{var[4]}"
                state_values = val
                lp_contents += f"{state_vars} = {state_values}\n"
            else:
                pass
        return lp_contents

    def make_model(self):
        """
        Build the MILP model to find the best linear trail
        """

        lp_contents = "\\ Linear attack on {} rounds of CLEFIA\n".format(self.nrounds)
        lp_contents += self.generate_objective_function()
        lp_contents += self.generate_constraints()
        lp_contents += self.exclude_trivial_trail()
        lp_contents += self.declare_fixed_variables()
        lp_contents += self.declare_binary_vars()
        lp_contents += "end"
        with open(self.lp_file_name, "w") as lp_file:
            lp_file.write(lp_contents)

    def exclude_the_previous_sol(self):
        '''
        Let x{S} be the binary variables. Suppose you have a binary
        solution x* in available from the most recent optimization.
        Let N be the subset of S such that x*[n] = 1 for all n in N
        Then, add the following constraint:
        sum{n in N} x[n] - sum{s in S-N} x[s] <= |N|-1
        '''

        all_vars = self.milp_model.getVars()
        nonzero_vars = [v for v in all_vars if v.x == 1]
        zero_vars = [v for v in all_vars if v.x == 0]
        support = len(nonzero_vars)
        first_term = sum(nonzero_vars)
        second_term = sum(zero_vars)
        lhs = first_term - second_term
        self.milp_model.addConstr(lhs <= support - 1)

    def solve(self):
        output = None
        self.milp_model = read(self.lp_file_name)
        os.remove(self.lp_file_name)
        if self.mode == 0:
            output = self.find_characteristic()
        elif self.mode == 1:
            self.find_multiple_characteristics(self.number_of_trails)
        elif self.mode == 2:
            output = self.compute_differential_effect()
            # self.compute_differential_effect_classic_method()
        else:
            print('Enter a number in [0, 1, 2], for the mode parameter please!')
        return output

    def parse_solver_output(self):
        """
        Extract the linear characteristic from the solver output
        """

        get_bit_value = lambda t: str(int(self.milp_model.getVarByName(t).Xn))
        characteristic = dict()
        for r in range(self.nrounds + 1):
            x = self.flatten_state(self.generate_round_x_variables(r))
            x_value = hex(int("0b" + "".join(list(map(lambda t: str(int(self.milp_model.getVarByName(t).Xn)), x))), 2))[2:].zfill(32)
            characteristic[f"x_{r}"] = x_value
        for r in range(self.nrounds):
            round_probability = 0
            for bn in range(2):
                for byten in range(4):
                    if (bn == 0 and byten % 2 == 0) or (bn == 1 and byten % 2 == 1):
                        for q in self.s0_correlations:
                            weight = float(self.s0_correlations[q])*self.milp_model.getVarByName(f"q_{r}_{bn}_{byten}_{q}").Xn
                            round_probability += weight
                    elif (bn == 0 and byten % 2 == 1) or (bn == 1 and byten % 2 == 0):
                        for q in self.s1_correlations:
                            weight = float(self.s1_correlations[q])*self.milp_model.getVarByName(f"q_{r}_{bn}_{byten}_{q}").Xn
                            round_probability += weight
            characteristic[f"pr_{r}"] = f"-{round_probability}"
        characteristic["total_weight"] = "%0.02f" % self.total_weight
        characteristic["nrounds"] = self.nrounds
        return characteristic

    @staticmethod
    def print_trail(trail):
        """
        Print out the discovered linear characteristic
        """

        header = ['x', 'pr']
        # Print everthing
        diff_trail_values = map(str, trail.values())
        col_width = max(len(s) for s in diff_trail_values) + 2
        header_str = "Rounds\t"
        data_str = ""
        current_row = 0
        for entry in header[0:-2]:
            header_str += entry.ljust(col_width)
        header_str += header[-2].ljust(col_width)
        header_str += header[-1].ljust(7)
        for r in range(trail["nrounds"] + 1):
            data_str += str(current_row) + '\t'
            data_str += trail.get(f"x_{r}", 'none').ljust(col_width)
            data_str += trail.get(f"pr_{r}", 'none').ljust(col_width)
            data_str += '\n'
            current_row += 1
        print(header_str)
        print("-"*len(header_str))
        print(data_str)
        total_weight = trail["total_weight"]
        print(f"Weight: -{total_weight}")
        return

    def find_characteristic(self):
        """
        Find the best linear trail for reduced-round CLEFIA
        """
        trail = None
        self.milp_model.Params.OutputFlag = True
        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        obj = self.milp_model.getObjective()
        # Consider the start_weight
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        time_start = time.time()
        #m.setParam(GRB.Param.Threads, 16)
        self.milp_model.optimize()
        # Gurobi syntax: m.Status == 2 represents the model is feasible.
        if (self.milp_model.Status == GRB.OPTIMAL or self.milp_model.Status == GRB.TIME_LIMIT or \
            self.milp_model.Status == GRB.INTERRUPTED):
            self.total_weight = self.milp_model.objVal
            print(f"\nThe squared correlation of the best linear characteristic: 2^-({self.total_weight})")
            print("\nLinear trail:\n")
            trail = self.parse_solver_output()
            self.print_trail(trail=trail)
        # Gurobi syntax: m.Status == 3 represents the model is infeasible. (GRB.Status.INFEASIBLE)
        elif self.milp_model.Status == GRB.INFEASIBLE:
            print("The model is infeasible!")
        else:
            print("Unknown error!")
        elapsed_time = time.time() - time_start
        print("Time used: %0.02f" % elapsed_time)
        return trail

    def find_multiple_characteristics(self, number_of_trails=2):
        """
        Find multiple linear trails for reduced-round of CLEFIA
        """

        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        obj = self.milp_model.getObjective()
        # Consider the start_weight
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        self.milp_model.Params.OutputFlag = False
        self.milp_model.Params.PoolSearchMode = 2
        # Limit number of solutions
        self.milp_model.Params.PoolSolutions = number_of_trails
        time_start = time.time()
        self.milp_model.optimize()
        if (self.milp_model.Status == GRB.OPTIMAL or self.milp_model.Status == GRB.TIME_LIMIT or \
            self.milp_model.Status == GRB.INTERRUPTED):
            # First Method:
            for sol_number in range(number_of_trails):
                if (self.milp_model.Status == GRB.OPTIMAL):
                    self.total_weight = self.milp_model.PoolObjVal
                    trail = self.parse_solver_output()
                    self.print_trail(trail=trail)
                elif (self.milp_model.Status == GRB.TIME_LIMIT or self.milp_model.Status == GRB.INTERRUPTED):
                    self.total_weight = self.milp_model.PoolObjVal
                    trail = self.parse_solver_output()
                    self.print_trail(trail=trail)
                    break
                else:
                    break
                self.exclude_the_previous_sol()
                print("#"*50)
                self.milp_model.optimize()
            # Second Method:
            # number_of_trails = self.milp_model.SolCount
            # for sol_number in range(number_of_trails):
            #     self.milp_model.Params.SolutionNumber = sol_number
            #     # PoolObjVal : This attribute is used to query the objective value of the <span>$</span>k<span>$</span>-th solution stored in the pool of feasible solutions found so far for the problem
            #     self.total_weight = self.milp_model.PoolObjVal
            #     trail = self.parse_solver_output()
            #     self.print_trail(trail=trail)
        # Gurobi syntax: m.Status == 3 represents the model is infeasible. (GRB.INFEASIBLE)
        elif self.milp_model.Status == GRB.INFEASIBLE:
            print("The model is infeasible!")
        else:
            print("Unknown error!")
        elapsed_time = time.time() - time_start
        print("Total time to find %s linear trails: %0.02f" % (number_of_trails, elapsed_time))

    def compute_differential_effect(self):
        """
        Compute the linear effect for a given input/output differences
        Some general information about Gurobi:
        PoolSolutions: It controls the size of the solution pool.
        Changing this parameter won't affect the number of solutions that are found -
        it simply determines how many of those are retained
        You can use the PoolSearchMode parameter to control the approach used to find solutions.
        In its default setting (0), the MIP search simply aims to find one optimal solution.
        Setting the parameter to 2 causes the MIP to do a systematic search for the n best solutions.
        With a setting of 2, it will find the n best solutions,
        where n is determined by the value of the PoolSolutions parameter
        SolCount: Number of solutions found during the most recent optimization.

        Model status:
        LOADED	1	Model is loaded, but no solution information is available.
        OPTIMAL	2	Model was solved to optimality (subject to tolerances), and an optimal solution is available.
        INFEASIBLE	3	Model was proven to be infeasible.
        """

        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        #self.milp_model.Params.PreSolve = 0 # Activating this flag causes the performance to be decreased, but the accuracy will be increased
        self.milp_model.Params.PoolSearchMode = 2
        self.milp_model.Params.PoolSolutions = 1
        self.milp_model.Params.OutputFlag = False

        self.milp_model.printStats()

        # Consider the start_weight
        obj = self.milp_model.getObjective()
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        time_start = time.time()
        self.milp_model.optimize()
        current_probability = 0
        if (self.milp_model.Status == GRB.OPTIMAL):
            self.total_weight = self.milp_model.objVal
            diff_prob = 0
            print('\n')
            while (self.milp_model.Status == GRB.OPTIMAL and self.total_weight <= self.end_weight):
                self.total_weight = self.milp_model.PoolObjVal
                self.milp_model.Params.PoolSolutions = 2000000000 #GRB.MAXINT
                temp_constraint = self.milp_model.addConstr(obj == self.total_weight, name='temp_constraint')
                # self.milp_model.Params.PoolGap = 0
                # self.milp_model.Params.PreSolve = 0
                # self.milp_model.printStats()
                self.milp_model.update()
                self.milp_model.optimize()
                diff_prob += math.pow(2, -self.total_weight) * self.milp_model.SolCount
                print(f"Current weight: {self.total_weight}")
                print(f"Number of trails: {self.milp_model.SolCount}")
                current_probability = math.log(diff_prob, 2)
                print(f"\tCurrent Probability: 2^({current_probability})")
                elapsed_time = time.time() - time_start
                print("Time used = %0.04f seconds\n" % elapsed_time)
                self.milp_model.remove(temp_constraint)
                self.milp_model.Params.PoolSolutions = 1
                self.milp_model.addConstr(obj >= (self.total_weight + self.eps), name='temp_cond')
                #self.milp_model.Params.PreSolve = 0
                self.milp_model.optimize()
        elif (self.milp_model.Status == GRB.INFEASIBLE):
            print("The model is infeasible!")
        else:
            print("Unknown Error!")
        return current_probability

    def compute_differential_effect_classic_method(self):
        """
        Compute linear effect by enumerating all possible linear trails
        """

        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        self.milp_model.Params.OutputFlag = False
        # self.milp_model.printStats()
        # Consider the start_weight
        obj = self.milp_model.getObjective()
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        time_start = time.time()
        self.milp_model.optimize()
        # self.milp_model.Params.Quad = 1
        sol_dict = dict()
        if (self.milp_model.Status == GRB.OPTIMAL):
            self.total_weight = self.milp_model.objVal
            diff_prob = 0
            print('\n')
            while (self.milp_model.Status == GRB.OPTIMAL and self.total_weight <= self.end_weight):
                self.total_weight = self.milp_model.objVal
                diff_prob += math.pow(2, -self.total_weight)
                total_weight_st = 'ntrails_%0.2f' % self.total_weight
                sol_dict[total_weight_st] = sol_dict.get(total_weight_st, 0) + 1
                print('Current weight: %s' % str(self.total_weight))
                print('Number of trails: %d' % sol_dict[total_weight_st])
                print('\tCurrent Probability: 2^(' + str(math.log(diff_prob, 2)) + ')')
                time_end = time.time()
                print('Time used = %0.4f seconds\n' % (time_end - time_start))
                self.exclude_the_previous_sol()
                self.milp_model.optimize()
        elif (self.milp_model.Status == GRB.INFEASIBLE):
            print('The model is infeasible!')
        else:
            print('Unknown Error!')

def loadparameters(args):
        """
        Get parameters from the argument list and inputfile.
        """
        # Load default values
        params = {"nrounds" : 2,
                  "mode" : 0,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : 3600,
                  "numberoftrails" : 1,
                  "fixedVariables" : {}}

        # Check if there is an input file specified
        if args.inputfile:
            with open(args.inputfile[0], 'r') as input_file:
                doc = yaml.load(input_file, Loader=yaml.FullLoader)
                params.update(doc)
                if "fixedVariables" in doc:
                    fixed_vars = {}
                    for variable in doc["fixedVariables"]:
                        fixed_vars = dict(list(fixed_vars.items()) +
                                        list(variable.items()))
                    params["fixedVariables"] = fixed_vars

        # Override parameters if they are set on commandline
        if args.nrounds:
            params["nrounds"] = args.nrounds[0]

        if args.startweight:
            params["startweight"] = args.startweight[0]

        if args.endweight:
            params["endweight"] = args.endweight[0]

        if args.mode:
            params["mode"] = args.mode[0]

        if args.timelimit:
            params["timelimit"] = args.timelimit[0]

        if args.numberoftrails:
            params["numberoftrails"] = args.numberoftrails[0]

        return params

def main():
    """
    Parse the arguments and start the request functionality with the provided
    parameters.
    """

    parser = ArgumentParser(description="This tool finds the best linear"
                                        "trail in a cryptographic primitive"
                                        "using Gurobi",
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('--startweight', nargs=1, type=int,
                        help="Starting weight for the trail search.")
    parser.add_argument('--endweight', nargs=1, type=int,
                        help="Stop search after reaching endweight.")
    parser.add_argument('--nrounds', nargs=1, type=int,
                        help="The number of rounds for the cipher")
    parser.add_argument('--mode', nargs=1, type=int,
                        choices=[0, 1], help=
                        "0 = search characteristic for fixed round\n"
                        "1 = determine the probability of the linear\n")
    parser.add_argument('--timelimit', nargs=1, type=int,
                        help="Set a timelimit for the search in seconds.")
    parser.add_argument('--inputfile', nargs=1, help="Use an yaml input file to"
                                                     "read the parameters.")
    parser.add_argument('--numberoftrails', nargs=1, type=int,
                        help="Number of trails.")

    # Parse command line arguments and construct parameter list.
    args = parser.parse_args()
    params = loadparameters(args)
    clefia = Lin(params)
    clefia.make_model()
    clefia.solve()

if __name__ == "__main__":
    main()