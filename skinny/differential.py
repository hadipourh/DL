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
"""

from argparse import ArgumentParser, RawTextHelpFormatter
import yaml
import time
from gurobipy import read
from gurobipy import GRB
import math
import os
from copy import deepcopy


"""
Modeling the differential analysis of SKINNY by MILP
This tool can:
 - find the best differential trail
 - find multiple differential trails
 - compute the differential probability taking the clustering effect into account

MILP variables:

    x_roundNumber_byteNumber_bitNumber
    x_roundNumber_byteNumber_0: msb
    x_roundNumber_byteNumber_7: lsb

    State before the Subcell after i rounds:
    x_i_0	x_i_1	x_i_2	x_i_3
    x_i_4	x_i_5	x_i_6	x_i_7
    x_i_8	x_i_9	x_i_10	x_i_11
    x_i_12	x_i_13	x_i_14	x_i_15

    State before the AddTweakey:
    y_i_0	y_i_1   y_i_2	y_i_3
    y_i_4	y_i_5	y_i_6	y_i_7
    y_i_8	y_i_9	y_i_10	y_i_11
    y_i_12	y_i_13	y_i_14	y_i_15

    State before the ShiftRows:
    z_i_0	z_i_1	z_i_2	z_i_3
    z_i_4	z_i_5	z_i_6	z_i_7
    z_i_8	z_i_9	z_i_10	z_i_11
    z_i_12	z_i_13	z_i_14	z_i_15

    State before the MixColumns:
    pz_i_0	pz_i_1	pz_i_2	pz_i_3
    pz_i_4	pz_i_5	pz_i_6	pz_i_7
    pz_i_8	pz_i_9	pz_i_10	pz_i_11
    pz_i_12	pz_i_13	pz_i_14	pz_i_15
    Note that pz is only a code variable and it does not appear in the MILP model

    TK1 in round i:
    tk1_i_0	    tk1_i_1	    tk1_i_2	    tk1_i_3
    tk1_i_4	    tk1_i_5	    tk1_i_6	    tk1_i_7
    tk1_i_8	    tk1_i_9	    tk1_i_10	tk1_i_11
    tk1_i_12	tk1_i_13	tk1_i_14	tk1_i_15

    TK2 in round i:
    tk2_i_0	    tk2_i_1	    tk2_i_2	    tk2_i_3
    tk2_i_4	    tk2_i_5	    tk2_i_6	    tk2_i_7
    tk2_i_8	    tk2_i_9	    tk2_i_10	tk2_i_11
    tk2_i_12	tk2_i_13	tk2_i_14	tk2_i_15

    TK3 in round i:
    tk3_i_0	    tk3_i_1	    tk3_i_2	    tk3_i_3
    tk3_i_4	    tk3_i_5	    tk3_i_6	    tk3_i_7
    tk3_i_8	    tk3_i_9	    tk3_i_10	tk3_i_11
    tk3_i_12	tk3_i_13	tk3_i_14	tk3_i_15

    Round tweakey in round i:
    tk_i_0	tk_i_1	tk_i_2	tk_i_3
    tk_i_4	tk_i_5	tk_i_6	tk_i_7
    
    The permuted twakey in round i:
    ptk1_i_0	    ptk1_i_1	    ptk1_i_2	    ptk1_i_3
    ptk1_i_4	    ptk1_i_5	    ptk1_i_6	    ptk1_i_7
    ptk1_i_8	    ptk1_i_9	    ptk1_i_10	    ptk1_i_11
    ptk1_i_12	    ptk1_i_13	    ptk1_i_14	    ptk1_i_15

    ptk2_i_0	    ptk2_i_1	    ptk2_i_2	    ptk2_i_3
    ptk2_i_4	    ptk2_i_5	    ptk2_i_6	    ptk2_i_7
    ptk2_i_8	    ptk2_i_9	    ptk2_i_10	    ptk2_i_11
    ptk2_i_12	    ptk2_i_13	    ptk2_i_14	    ptk2_i_15

    ptk3_i_0	    ptk3_i_1	    ptk3_i_2	    ptk3_i_3
    ptk3_i_4	    ptk3_i_5	    ptk3_i_6	    ptk3_i_7
    ptk3_i_8	    ptk3_i_9	    ptk3_i_10	    ptk3_i_11
    ptk3_i_12	    ptk3_i_13	    ptk3_i_14	    ptk3_i_15
    ptk1, ptk2, and ptk3 are only used in the code, not in the MILP model

    Sbox indicators variables:
    q_roundNumber_byteNumber : shows the activity of an Sbox
    q2_roundNumber_byteNumber
    q2_4150_roundNumber_byteNumber
    q2_6781_roundNumber_byteNumber
    q3_roundNumber_byteNumber
    q3_1926_roundNumber_byteNumber
    q3_4150_roundNumber_byteNumber
    q3_6781_roundNumber_byteNumber
    q4_roundNumber_byteNumber
    q4_4150_roundNumber_byteNumber
    q5_roundNumber_byteNumber
    q5_4150_roundNumber_byteNumber
    q6_roundNumber_byteNumber
    q7_roundNumber_byteNumber
    q_roundNumber_byteNumber = q2_roundNumber_byteNumber + q2_4150_roundNumber_byteNumber + q2_6781_roundNumber_byteNumber
                             + q3_roundNumber_byteNumber + q3_1926_roundNumber_byteNumber + q3_4150_roundNumber_byteNumber
                             + q3_6781_roundNumber_byteNumber + q4_roundNumber_byteNumber + q4_4150_roundNumber_byteNumber
                             + q5_roundNumber_byteNumber + q5_4150_roundNumber_byteNumber + q6_roundNumber_byteNumber + 
                             + q7_roundNumber_byteNumber

    Objective function: Minimize Sum(2*q2_roundNumber_byteNumber + 2.4150*q2_4150_roundNumber_byteNumber + 2.6781*q2_6781_roundNumber_byteNumber
                                   + 3*q3_roundNumber_byteNumber + 3.1926*q3_1926_roundNumber_byteNumber + 3.4150*q3_4150_roundNumber_byteNumber
                                   + 3.6781*q3_6781_roundNumber_byteNumber + 4*q4_roundNumber_byteNumber + 4.4150*q4_4150_roundNumber_byteNumber
                                   + 5*q5_roundNumber_byteNumber + 5.4150*q5_4150_roundNumber_byteNumber + 6*q6_roundNumber_byteNumber + 
                                   + 7*q7_roundNumber_byteNumber)
"""


class Differential:
    '''
    Convert the differential analysis of SKINNY to an MILP problem
    '''

    count = 0
    def __init__(self, param, exact=True):
        self.variant = param['variant']
        if self.variant not in [0, 1, 2, 3, 4]:
            raise ValueError('variant should be in [0, 1, 2, 3, 4]')
        self.cellsize = param['cellsize']
        if self.cellsize not in [4, 8]:
            raise Exception('cellsize must be 4 or 8')
        self.skipsb = param['skipsb']
        self.rounds = param['rounds']
        self.start_round = param['start_round']
        self.end_round = param['end_round']
        self.upperbound1 = param['upperbound1']
        self.upperbound2 = param['upperbound2']
        self.start_weight = param['sweight']
        self.end_weight = param['endweight']
        self.time_limit = param['timelimit']        
        self.mode = param['mode']
        self.fixed_variables = param['fixedVariables']
        self.exact = exact #A Boolean variable indicating whether the model is exact or not. 
        self.accuracy_threshold = 7
        self.total_weight = None
        self.big_m = 2*8         
        self.eps = 1e-2
        self.obj_func = ''
        self.used_variables = []
        self.permuteation = [0x0, 0x1, 0x2, 0x3, 0x7, 0x4, 0x5, 0x6, 0xa, 0xb, 0x8, 0x9, 0xd, 0xe, 0xf, 0xc]
        self.tk_permutation = [0x9, 0xf, 0x8, 0xd, 0xa, 0xe, 0xc, 0xb, 0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7]
        self.model_filename = f"SKINNY-{self.cellsize*16}-{self.cellsize*self.variant*16}-{self.rounds}r.lp"

        # Precomputed constraints for the 8-bit S-box
        # a0, b0: msb
        # a7, b7: lsb
        self.q2_ineqs = ['- a6 - a7 >= -1', 'a6 - b4 >= 0', 'a7 - b2 >= 0', 'a3 - b1 >= 0', '- a4 + b3 >= 0', '- a5 + b7 >= 0', '- b3 - b7 >= -1', 'a2 - b0 >= 0', 'a0 - b6 >= 0', '- a1 + b5 >= 0', '- a3 - a7 >= -1', '- a3 - a6 >= -1', '- a2 - b7 >= -1', '- a0 - b3 >= -1', 'a4 - a6 + b4 >= 0', 'a4 + a5 - a7 + b2 >= 0', '- b5 - b6 >= -1', '- b3 - b5 >= -1', '- a2 - a4 >= -1', '- a0 - a2 >= -1', '- b0 - b1 >= -1', '- b3 - b4 >= -1', 'a1 + a4 + a5 + b0 + b1 + b2 + b4 + b6 >= 1', '- a3 + b6 - b7 >= -1', '- b2 - b3 >= -1', '- b5 - b7 >= -1', '- a0 - b1 >= -1', '- a2 - b5 >= -1', '- a2 - b4 >= -1', '- b2 - b7 >= -1', '- b0 - b2 >= -1', '- a0 - a5 >= -1', '- a0 - b4 >= -1', 'a0 - a3 + b1 + b5 >= 0', '- a0 - b2 >= -1', '- b1 - b5 >= -1', '- b2 - b5 >= -1', '- b4 - b5 >= -1', '- a5 - b4 >= -1', '- a4 - b1 >= -1']
        self.q2_4150_ineqs = ['a3 - b7 >= 0', '- b6 >= 0', '- b5 >= 0', '- b2 >= 0', '- b1 >= 0', '- a7 >= 0', '- a5 >= 0', '- a4 >= 0', '- b3 - b4 >= -1', '- a0 - a1 >= -1', 'a0 - a2 >= 0', 'a1 - a6 >= 0', 'a2 - b0 >= 0', 'b4 - b7 >= 0', '- a3 + b0 + b7 >= 0', 'a3 + b4 >= 1', 'a6 - b4 >= 0']
        self.q2_6781_ineqs = ['- b5 >= 0', '- b3 >= 0', '- b2 >= 0', '- b1 >= 0', '- b0 >= 0', '- a7 >= 0', '- a5 >= 0', '- a4 >= 0', '- a2 >= 0', 'a0 >= 1', 'a1 >= 1', 'a6 >= 1', 'b4 >= 1', 'b6 >= 1', 'a3 - b7 >= 0', '- a3 + b7 >= 0']
        self.q3_ineqs = ['- b0 - b4 >= -1', '- a1 - b2 >= -1', '- b2 - b5 >= -1', '- a4 - b0 >= -1', '- a6 - b0 >= -1', '- b1 - b6 >= -1', '- a0 + b5 + b6 >= 0', 'a4 + b0 + b1 - b3 >= 0', '- a2 + b0 + b1 + b2 >= 0', '- a4 + b1 + b3 >= 0', 'a5 + b4 + b6 - b7 >= 0', '- a5 + b4 + b6 + b7 >= 0', 'a0 + a1 - a3 + b1 >= 0', '- a7 - b6 >= -1', 'a4 - a6 + b2 + b4 >= 0', '- a0 + a1 - b5 >= -1', 'a3 - a7 - b1 >= -1', '- b1 - b7 >= -1', '- a1 + a5 + a7 + b5 >= 0', '- b1 - b4 >= -1', '- a1 - b3 - b7 >= -2', 'b0 - b3 - b6 >= -1', '- a4 + a7 - b5 >= -1', 'a4 + a6 + a7 - b2 >= 0', '- a5 - b0 >= -1', 'a1 + a5 + a6 - b5 >= 0', 'a5 - a6 - a7 - b4 >= -2', 'a2 + a6 + b1 - b3 - b5 >= -1', '- a1 + a6 - a7 + b7 >= -1', 'a4 + a5 + b0 + b2 + b4 + b5 >= 1', '- a6 + a7 - b2 + b7 >= -1', 'a6 + b2 + b3 - b4 >= 0', '- a5 - a6 - b5 >= -2', 'a1 + a6 + b0 + b1 + b3 + b4 + b5 + b7 >= 1', 'a2 + a6 + b1 + b3 + b4 + b6 + b7 >= 1', 'a1 + a3 + b2 + b4 + b5 + b6 + b7 >= 1', '- a1 - a5 - b5 >= -2', '- a0 - a3 + b6 - b7 >= -2', '- a5 - b3 - b4 >= -2', '- a0 - b1 + b3 >= -1', '- a0 - b0 - b6 >= -2', 'a0 + b0 + b5 - b6 >= 0', 'a0 + a1 + a2 + a4 + b1 + b2 + b4 + b5 >= 1', 'a0 - a3 - a6 + b3 - b5 + b7 >= -2', '- a1 - a7 - b5 >= -2', '- a5 - b4 - b6 >= -2', '- a6 - a7 - b1 + b5 >= -2', '- b2 - b6 >= -1', '- b2 - b3 - b7 >= -2', '- a5 - b2 - b4 >= -2', '- a4 - b1 - b2 >= -2', '- b0 - b5 - b6 >= -2', '- b0 - b1 - b5 >= -2', 'a4 + a5 + a7 - b4 + b5 + b6 >= 0', '- a1 - a6 + b3 + b5 >= -1', '- a1 - b4 - b6 >= -2', 'a2 - a6 - b1 >= -1', 'a0 + a7 - b1 - b3 - b5 >= -2', '- a7 - b0 - b1 >= -2', '- a2 - a4 + a6 + b2 >= -1', '- b2 - b3 - b4 >= -2', '- a6 - a7 - b2 - b7 >= -3', '- a1 + b5 - b6 - b7 >= -2', 'a0 - a1 + a3 - b4 - b7 >= -2', '- a0 + a3 - b4 + b6 + b7 >= -1', '- a5 - b5 - b6 - b7 >= -3', '- a2 - b0 - b2 - b3 >= -3']
        self.q3_1926_ineqs = ['- b7 >= 0', '- b5 >= 0', '- b4 >= 0', '- b2 >= 0', '- b1 >= 0', '- a7 >= 0', '- a6 >= 0', '- a5 >= 0', '- a4 >= 0', '- a3 >= 0', '- a1 >= 0', 'a0 >= 1', 'a2 >= 1', 'b0 >= 1', 'b6 >= 1']
        self.q3_4150_ineqs = ['- b5 >= 0', 'a0 - b6 >= 0', '- b2 - b4 >= -1', 'b4 - b7 >= 0', 'a1 - a6 >= 0', '- a1 - b0 >= -1', '- b1 - b4 >= -1', '- a3 - a4 + b7 >= -1', '- a0 + b0 + b6 >= 0', 'a3 - a5 + b7 >= 0', 'a6 - b4 >= 0', 'a0 + a4 + a5 >= 1', '- a4 + b3 >= 0', 'a3 - a7 + b4 >= 0', 'a7 - b2 >= 0', '- a5 + b4 >= 0', 'a4 - a7 - b6 >= -1', 'a7 - b3 - b6 >= -1', '- a4 + a7 >= 0', 'b1 + b2 + b4 >= 1', '- a3 - b1 >= -1', '- a2 - b4 >= -1', 'a3 + b3 - b6 + b7 >= 0', '- a5 - b3 >= -1', '- a3 + b3 - b7 >= -1', 'a3 - b3 - b7 >= -1', 'a2 + b1 + b3 + b4 >= 1', '- a2 - b2 - b3 >= -2']
        self.q3_6781_ineqs = ['- b5 >= 0', '- b4 >= 0', '- b2 >= 0', '- b1 >= 0', '- a7 >= 0', '- a6 >= 0', '- a5 >= 0', '- a4 >= 0', '- a3 >= 0', '- a1 >= 0', 'a0 >= 1', 'a2 >= 1', 'b0 >= 1', 'b6 >= 1', 'b7 >= 1']
        self.q4_ineqs = ['- a6 - b2 - b6 >= -2', '- a0 + b0 + b5 + b6 >= 0', 'a5 + b4 + b6 - b7 >= 0', '- a5 + b4 + b6 + b7 >= 0', 'a4 - a6 + b2 + b4 >= 0', 'a4 + a5 - a7 + b2 >= 0', 'a4 + a6 + b2 - b4 >= 0', 'a1 + a5 + a6 - b5 >= 0', 'a0 + b0 + b5 - b6 >= 0', 'a4 + b0 + b1 - b3 >= 0', '- a2 + b0 + b1 + b2 >= 0', '- a1 + a5 + a6 + b5 >= 0', '- a5 - b1 - b4 >= -2', '- a4 + b0 + b1 + b3 >= 0', 'a2 - b0 + b1 + b2 >= 0', 'a0 + a1 - a3 + b1 >= 0', '- a0 + a1 - b1 - b4 >= -2', '- b0 - b1 - b4 >= -2', '- a5 - a6 - b0 + b2 >= -2', 'a0 + a1 + a3 - b1 >= 0', '- a1 + a7 - b2 - b5 >= -2', '- b0 - b1 - b5 - b6 >= -3', '- a1 - a7 - b5 - b6 >= -3', '- a1 - b2 - b4 + b5 >= -2', '- a1 - a6 - b0 + b4 >= -2', '- a0 + a1 - b2 - b5 >= -2', '- a0 + a1 - b0 - b5 >= -2', '- a4 - b1 + b4 + b5 - b6 >= -2', 'a4 + a6 + a7 - b2 + b4 >= 0', '- a6 - b2 - b4 - b5 >= -3', '- a0 + a5 - a7 - b4 + b6 >= -2', '- a5 - b3 - b4 - b5 >= -3', '- a4 - b2 - b6 >= -2', '- a4 - b0 - b7 >= -2', '- a1 - a5 + a7 - b4 - b6 >= -3', '- a5 - b3 - b4 - b6 >= -3', '- a3 - b0 + b1 + b5 - b6 >= -2', '- a2 - a4 - a5 + a6 + b2 >= -2', '- a0 + a3 - a4 - b0 + b1 >= -2', 'a2 - a4 + a6 - b1 - b4 >= -2', 'a0 - a3 - a5 - b4 - b7 >= -3', '- a1 - a5 - b1 >= -2', '- a0 - a4 - a5 + b6 >= -2', '- a1 - a4 + a5 + a7 - b6 >= -2', 'a7 - b0 - b2 - b4 >= -2', '- a4 - b1 - b3 - b6 >= -3', 'a0 - a1 + a3 - b0 + b6 - b7 >= -2', 'a0 - a1 + a3 - a5 + b7 >= -1', '- a5 - b0 - b1 >= -2', 'a0 - a1 - a3 - b1 + b6 - b7 >= -3', '- a2 - a4 - a6 + b2 - b4 >= -3', '- b1 - b2 + b5 - b6 >= -2', 'a2 + a4 - b0 + b1 + b3 >= 0', '- a2 + a4 + b1 - b2 - b3 >= -2', '- a5 + a6 - a7 + b3 - b4 >= -2', 'a4 + a5 + a6 + a7 + b0 + b6 >= 1', 'a1 + a4 - b1 - b3 - b5 + b6 >= -2', '- a1 + a5 - a6 - a7 - b5 >= -3', 'a0 - a4 + b0 + b3 - b5 + b6 >= -1', '- a0 + a4 + b0 - b1 + b3 + b6 >= -1', '- a5 - b2 - b3 - b4 >= -3', '- a0 - a7 - b5 - b6 >= -3', 'a1 + a5 - a6 - a7 + b2 + b5 >= -1', 'a2 - a6 + b0 - b1 + b2 + b4 >= -1', 'a1 - a4 + a5 + a7 + b2 - b5 >= -1', 'a0 - a1 - a7 + b0 - b2 - b3 >= -3', '- a1 - a5 - b2 - b3 >= -3', 'a4 + a7 - b2 - b5 + b7 >= -1', 'a3 - b0 + b2 + b4 + b5 - b6 + b7 >= -1', 'a0 + a4 - b1 + b3 - b5 - b6 >= -2', '- a4 - b0 - b4 >= -2', '- a0 + a5 + a6 - b0 + b5 - b7 >= -2', '- a3 + a7 - b0 - b2 + b7 >= -2', 'a0 - a3 + a4 - a6 - b0 + b7 >= -2', 'a1 - a5 + a6 - b3 + b5 - b6 + b7 >= -2', '- a7 - b0 + b2 - b6 >= -2', '- a0 - a5 - b1 >= -2', '- a1 - a4 + a7 - b1 + b5 >= -2', '- a4 - a7 - b2 - b5 >= -3', 'a2 - a4 + a6 - b0 + b2 >= -1', '- a1 + a7 + b0 + b1 + b3 + b4 + b7 >= 0', '- a0 - a1 + a5 + b1 - b4 + b5 >= -2', 'a1 - a6 + a7 - b2 + b5 + b7 >= -1', 'a1 + a7 + b4 - b5 - b6 - b7 >= -2', 'a3 + a4 - b0 - b1 + b5 >= -1', '- a1 + a7 + b4 + b5 - b6 - b7 >= -2', '- a1 - a5 - b0 >= -2', 'a0 + a1 + b0 + b2 + b4 + b6 + b7 >= 1', '- a1 - b1 - b2 + b5 >= -2', '- a4 - b2 - b5 - b7 >= -3', '- a0 - a4 - b1 - b3 - b5 >= -4', '- a6 - b0 - b2 - b4 >= -3', '- a0 + a1 - a5 - b4 - b5 >= -3', '- a1 - a5 - b2 - b5 >= -3', '- a0 + a3 - a6 + b1 + b5 + b6 >= -1', '- a1 + a3 + a5 + b0 - b4 + b5 + b7 >= -1', 'a4 - a6 - a7 + b4 - b7 >= -2', '- a4 - b1 - b2 - b7 >= -3', 'a0 - a1 + a4 - a6 - b1 - b3 >= -3', '- a0 + a3 - b0 + b4 + b6 - b7 >= -2', 'a1 - a5 - b2 + b5 - b6 + b7 >= -2', '- a0 + a1 - a6 - a7 + b4 - b7 >= -3', '- a6 - a7 - b1 + b3 - b4 + b6 >= -3', '- a1 - a2 - b4 - b6 >= -3', 'a2 - a5 - b1 - b6 - b7 >= -3', '- a4 - b0 - b1 - b2 >= -3', '- a5 + b0 + b1 + b2 + b3 + b5 + b6 >= 0', 'a5 + b0 + b1 + b2 - b3 + b5 + b6 >= 0', 'a1 + a4 + a5 + b0 + b1 + b2 + b5 >= 1', 'a0 + a3 - b1 - b4 + b7 >= -1', '- a2 - a4 + a6 + b4 - b5 >= -2', '- a0 - a3 - b1 + b5 + b6 >= -2', '- a1 - a5 - b3 - b5 >= -3', 'a6 + a7 + b0 + b1 + b3 + b5 - b6 >= 0', '- a0 - a1 - a3 + b0 + b3 - b4 + b6 + b7 >= -3', '- a0 - a1 + a3 + b0 + b3 - b4 + b6 - b7 >= -3', '- a1 - a3 + a5 + b0 - b4 + b5 - b7 >= -3', 'a3 + a4 - b1 - b3 + b6 - b7 >= -2', 'a1 + b0 + b1 + b4 + b5 + b6 + b7 >= 1', '- a6 - b1 - b2 - b4 >= -3', '- a2 - a6 - b0 - b1 >= -3', 'a2 - a4 - b0 - b2 - b3 >= -3', '- a0 - a1 + a3 + a5 + b1 - b4 + b6 + b7 >= -2', '- a2 - a4 + b1 - b2 + b3 >= -2', '- a0 - b1 - b3 - b5 - b6 >= -4', '- a0 - a1 - a3 + a5 + b1 + b6 - b7 >= -3', '- b2 - b5 - b6 - b7 >= -3', 'a3 + b1 - b2 - b4 - b5 - b7 >= -3', 'a4 + a6 + b1 - b2 + b4 + b5 + b6 + b7 >= 0', '- a3 + a4 + b2 + b4 + b5 + b6 + b7 >= 0', '- b0 - b1 - b2 - b5 >= -3', '- a5 - a6 - a7 + b4 - b5 + b7 >= -3', '- a6 + a7 - b2 - b4 >= -2', 'a0 - a3 + b1 - b3 - b4 + b7 >= -2', 'a0 - a1 + a3 + b1 - b3 - b4 - b7 >= -3', 'a5 + a6 + b0 + b1 + b3 + b5 + b6 >= 1', 'a1 + b0 + b1 + b3 + b4 + b5 + b6 >= 1', 'a4 + b1 + b2 + b4 - b5 + b6 + b7 >= 0', '- a0 - a7 - b0 + b2 + b7 >= -2', 'a0 + a4 + a5 + a7 + b1 + b4 + b5 >= 1', '- a3 + b1 - b2 - b4 + b6 + b7 >= -2', '- a0 + a1 - a5 + a7 - b3 + b7 >= -2', 'a4 - b1 + b3 - b4 - b6 >= -2', 'a4 + a5 + a6 + b0 + b4 + b5 + b6 >= 1', 'a1 + a3 - a4 + b2 + b4 + b5 + b6 - b7 >= -1', '- a1 + a5 + b0 + b1 + b2 + b3 - b5 + b6 >= -1', 'a5 + b0 + b1 + b2 + b3 + b4 >= 1', '- a4 - b1 - b2 - b4 >= -3', 'a0 + a1 + a5 + b0 + b1 + b2 + b3 >= 1', '- a1 - b2 - b6 - b7 >= -3', '- a0 - b0 - b2 - b4 >= -3', 'a0 + a1 + b0 + b1 + b2 + b3 + b4 + b6 >= 1', '- a0 - a3 - b1 - b4 + b6 + b7 >= -3', '- a0 + a2 - b0 + b2 - b6 >= -2', '- a1 - a5 - a6 - a7 + b4 + b7 >= -3', '- a1 + a6 - b3 + b5 - b6 - b7 >= -3', 'a0 + a7 + b0 + b1 + b3 + b4 + b7 >= 1', 'a6 - a7 + b1 - b3 - b5 - b6 - b7 >= -4', '- a2 + a6 - a7 + b2 - b6 + b7 >= -2']
        self.q4_4150_ineqs = ['a1 - b5 >= 0', 'a0 - b6 >= 0', '- a6 + b4 >= 0', 'a1 + b0 >= 1', '- a5 - b5 >= -1', '- a4 + a5 + b2 >= 0', '- a4 - a7 - b2 >= -2', '- a3 - b1 - b4 + b7 >= -2', '- a5 + b4 >= 0', 'a3 + a6 + a7 + b7 >= 1', 'a0 + a5 + a6 >= 1', '- a3 + b2 + b6 - b7 >= -1', '- a7 - b1 - b6 >= -2', 'a6 + b2 + b3 - b4 >= 0', '- a6 - b1 >= -1', 'a1 + a7 - b4 >= 0', '- a2 + a6 + b1 - b3 - b4 >= -2', 'a3 - b1 - b7 >= -1', 'a3 - b0 - b4 + b6 + b7 >= -1', 'a3 - b0 + b5 + b6 - b7 >= -1', 'a5 + a6 + a7 - b2 >= 0', 'a2 - b0 + b1 + b3 >= 0', 'a2 - b0 + b1 + b2 >= 0', '- a6 - b0 - b2 >= -2', '- a0 + b0 + b5 + b6 >= 0', '- a3 + a6 + a7 - b4 - b7 >= -2', '- a3 - a5 - a6 - b7 >= -3', 'a4 + b0 + b1 - b3 >= 0', '- a1 + b4 >= 0', '- a1 - b0 + b5 + b6 >= -1', 'a3 + b2 + b6 + b7 >= 1', 'a3 + a5 - a6 - b7 >= -1', '- a3 + a5 - a6 + b7 >= -1', '- a1 + a6 - b2 - b6 >= -2', '- a3 - b0 - b1 + b6 >= -2', '- a7 + b2 + b3 >= 0', '- a3 + b1 - b5 - b7 >= -2', 'a3 + b1 - b5 + b7 >= 0', 'a3 - a5 - a6 + b7 >= -1', '- a5 + a6 - a7 - b2 >= -2', '- a2 - a5 + b2 >= -1', '- a6 + a7 + b0 + b3 >= 0', 'a3 - b0 + b2 - b7 >= -1', '- a3 + a6 - a7 - b6 + b7 >= -2', '- a4 - a5 - b2 >= -2', '- a7 - b0 + b2 >= -1', '- a3 - a7 - b0 - b6 >= -3', 'b1 + b4 + b7 >= 1', '- b2 + b4 - b7 >= -1', '- b1 + b3 - b5 >= -1', '- a2 - b1 - b6 >= -2', 'a3 + a6 - a7 + b2 - b6 - b7 >= -2']
        self.q5_ineqs = ['- a0 + a3 + a5 + a7 - b0 - b2 - b7 >= -3', '- a4 + b0 + b1 + b3 >= 0', 'a4 - a6 + b2 + b4 >= 0', 'a4 + b0 + b1 - b3 >= 0', '- a4 - a5 - b0 - b4 >= -3', '- a1 + a5 + a6 + b5 >= 0', '- a0 + b0 + b5 + b6 >= 0', 'a0 + b0 + b5 - b6 >= 0', 'a4 + a5 - a7 + b2 >= 0', '- a5 + b4 + b6 + b7 >= 0', 'a2 - b0 + b1 + b2 >= 0', 'a5 + b4 + b6 - b7 >= 0', '- a2 + b0 + b1 + b2 >= 0', 'a1 + a5 + a6 - b5 >= 0', 'a4 + a5 + a7 - b2 >= 0', 'a0 + a1 + a3 - b1 >= 0', 'a0 + a1 - a3 + b1 >= 0', '- a1 - a4 - a5 - b0 >= -3', '- a0 + a1 - a5 - b1 - b4 >= -3', '- a6 - b0 - b1 - b2 - b4 >= -4', '- a0 - a4 - a5 - b2 - b5 >= -4', '- a0 - a6 - b2 - b5 - b6 >= -4', 'a0 + a1 + b0 + b2 + b4 >= 1', 'a4 - a5 + a6 - a7 - b4 >= -2', 'a0 - a1 - a6 - b2 - b6 >= -3', 'a4 + a6 + b2 - b4 >= 0', '- a5 - b0 - b1 - b4 >= -3', '- a0 + a1 + a4 - b0 - b1 - b5 >= -3', 'a0 + a1 + a5 + b0 + b2 >= 1', 'a6 + a7 - b1 - b2 + b5 - b6 >= -2', '- a4 + a6 - b0 - b1 - b2 - b5 >= -4', '- a3 - a4 - a5 - b2 - b4 >= -4', '- a0 - a4 + b3 - b5 - b6 >= -3', 'a1 - a3 - a6 - b0 - b2 - b6 >= -4', 'a4 - a6 + a7 - b2 - b4 >= -2', 'a4 - a5 - a6 - a7 + b4 >= -2', 'a4 + a6 + a7 - b2 + b4 >= 0', '- a2 + a4 + b1 - b2 - b3 >= -2', '- a2 - a4 + b1 - b2 + b3 >= -2', '- a1 - b0 - b1 - b4 - b5 - b6 >= -5', 'a2 - a4 - b0 + b1 - b3 >= -2', 'a2 + a4 - b0 + b1 + b3 >= 0', '- a1 - a5 - b0 - b1 >= -3', 'a1 - a6 - b1 - b2 + b5 - b6 >= -3', '- a1 - a5 + a7 - b0 - b2 - b6 >= -4', '- a1 - a6 - b0 - b1 - b2 >= -4', '- a1 - a4 - b2 - b5 - b6 >= -4', '- a1 - a4 + a5 + a7 + b2 + b5 >= -1', '- a1 + a2 - a6 - b0 + b2 - b6 >= -3', 'a1 + a5 - a6 - a7 + b2 + b5 >= -1', '- a1 - a4 + a5 - a7 - b2 + b5 >= -3', '- a0 - a4 - b0 - b5 - b6 >= -4', '- a1 + a5 + a7 + b0 + b2 - b5 + b6 >= -1', '- a1 + a5 - a6 - a7 + b2 - b5 >= -3', '- a0 + a1 - b2 - b4 - b5 >= -3', '- a0 + a5 + a7 - b2 - b4 + b6 >= -2', '- a4 - b0 - b1 - b2 - b7 >= -4', '- a2 - a6 - b0 - b1 + b2 + b4 >= -3', '- a2 - a4 - a6 + b0 + b2 - b4 >= -3', '- a0 + a4 + b0 - b3 - b5 - b6 >= -3', 'a2 - a4 + a6 - b0 + b2 + b4 >= -1', 'a2 + a6 + b0 - b1 + b2 - b4 >= -1', '- a2 - a4 + a6 + b0 + b2 + b4 >= -1', 'a1 + a5 - a6 + a7 - b2 + b5 >= -1', '- a0 + a1 - a5 - b3 - b4 - b5 >= -4', 'a1 - a4 + a5 + a7 + b2 - b5 >= -1', 'a2 - a6 + b0 - b1 + b2 + b4 >= -1', 'a0 + a1 + a4 + b0 + b4 + b6 >= 1', '- a0 - a4 + b0 - b1 - b3 + b6 >= -3', 'a1 - a4 + a5 - a7 - b2 - b5 >= -3', 'a0 - a4 + b0 - b1 - b3 - b6 >= -3', '- a0 + a4 + b0 - b1 + b3 + b6 >= -1', 'a6 + b0 + b1 + b3 + b4 + b5 >= 1', 'a0 + a4 + b0 - b1 + b3 - b6 >= -1', 'a0 + a4 + b0 - b3 - b5 + b6 >= -1', 'a5 - a6 + a7 + b0 + b3 - b4 + b6 >= -1', 'a0 - a4 + b0 + b3 - b5 + b6 >= -1', '- a1 - a5 - a7 + b1 + b2 + b5 - b6 >= -3', '- a6 - b1 - b2 + b4 + b5 - b6 >= -3', '- a0 + a5 - a7 - b0 + b2 - b4 + b6 >= -3', 'a1 - a3 - a5 - b0 + b1 + b5 - b6 >= -3', '- a1 - a4 - a5 - b3 - b4 - b6 >= -5', '- a4 - a7 - b1 - b2 + b5 - b6 >= -4', 'a1 - a5 + a7 + b2 + b4 + b5 + b7 >= 0', '- a1 - a3 - a5 - b0 - b4 + b5 - b7 >= -5', 'a0 - a1 - a4 - b0 - b6 >= -3', '- a2 + a6 - b0 - b1 + b2 - b4 >= -3', 'a1 - a4 - a5 + a6 + b4 + b5 + b7 >= -1', 'a1 - a5 - a7 - b2 + b4 + b5 + b7 >= -2', '- a0 - a1 + a3 - b0 - b1 + b6 - b7 >= -4', '- a0 - a3 - b0 - b1 - b4 + b6 + b7 >= -4', 'a2 - a4 - a6 - b0 + b2 - b4 >= -3', '- a1 - a5 + a7 + b4 - b5 + b7 >= -2', '- a1 - a5 - b0 - b5 - b6 + b7 >= -4', 'a1 + a7 + b2 + b4 - b5 - b6 - b7 >= -2', 'b0 + b2 + b4 + b6 + b7 >= 1', '- a0 - a1 + a3 - a5 + b1 - b4 + b6 - b7 >= -4', '- a0 - a1 - a3 - a5 + b1 + b6 + b7 >= -3', 'a0 - a1 + a4 - a5 - b1 - b3 - b4 >= -4', 'a1 - a4 + a6 + b4 - b5 - b6 - b7 >= -3', 'a4 + a5 + a6 + b0 + b6 >= 1', 'a1 - a5 - a7 - b2 + b4 - b5 - b6 - b7 >= -5', '- a0 + a1 + a3 - a6 + b1 + b5 + b6 >= -1', 'a0 - a1 - a4 - a5 + b3 - b4 >= -3', '- a1 + a5 - a6 + a7 - b2 - b5 >= -3', '- a0 - a3 + a5 - a7 - b0 - b1 - b7 >= -5', 'a0 + a1 + a4 + b1 + b4 + b6 >= 1', 'b0 + b1 - b2 + b3 + b5 + b6 >= 0', '- a1 - a6 - b0 + b4 - b6 >= -3', 'a1 + a3 - b0 - b1 - b4 - b6 >= -3', 'a0 + a1 + a4 + a6 - a7 + b1 - b2 >= -1', '- a1 - a6 + b0 + b1 - b2 - b4 - b6 >= -4', 'a3 - a5 - b0 - b1 - b6 >= -3', '- a4 - a5 - b2 - b4 - b6 >= -4', 'a4 - a5 - b1 + b3 - b4 - b6 >= -3', '- a1 - a4 - a5 - b1 - b2 >= -4', 'a4 + a5 + b0 + b2 + b4 >= 1', '- a1 + a3 - a5 - b0 - b4 + b5 + b7 >= -3', '- a1 + a3 - b0 - b1 - b4 + b5 + b7 >= -3', '- a3 - a4 - b0 + b1 - b6 + b7 >= -3', '- a2 - a5 - b0 - b1 - b5 - b6 - b7 >= -6', '- a0 - b0 - b1 - b2 - b4 >= -4', '- a0 + a3 - a5 - b1 + b6 + b7 >= -2', '- a0 - a3 - a5 - b1 - b4 + b6 - b7 >= -5', 'a1 + a3 - a4 + a5 + a7 - b0 - b1 + b4 >= -2', '- a0 + a3 - a4 - b0 + b1 - b4 + b6 + b7 >= -3', '- a1 - a2 - a7 + b2 + b4 - b5 - b6 - b7 >= -5', '- a3 - a4 - a7 - b0 - b2 - b6 >= -5', 'a3 - b0 - b1 - b2 - b6 >= -3', '- a0 + a1 - a3 - a4 + a7 - b0 - b1 + b6 >= -4', '- a1 - a4 - b2 + b4 - b5 - b7 >= -4', 'a4 + a5 + a6 + b1 - b5 + b6 >= 0', 'a0 + a3 - a5 - b1 - b4 - b7 >= -3', 'a0 - a1 - a3 - a5 - b1 + b6 + b7 >= -3', '- a0 + a3 + a5 - b1 - b2 + b6 - b7 >= -3', '- a0 + a1 + a3 + b1 + b5 + b6 - b7 >= -1', '- a1 + a3 - a5 + b0 - b2 - b3 - b4 + b5 - b7 >= -5', 'a0 + a1 + b0 + b4 + b6 + b7 >= 1', 'a0 + a3 + a5 - b1 - b4 - b5 + b6 + b7 >= -2', 'a0 - a1 - a3 + a5 - b1 - b5 + b6 - b7 >= -4', 'a5 + b0 + b1 - b2 + b5 + b6 >= 0', '- a1 + a3 + b0 + b1 - b3 - b4 + b5 + b7 >= -2', 'a0 - a1 + a3 - a5 + b1 - b5 + b6 + b7 >= -2', 'a0 - a3 - a5 + b1 - b4 - b5 + b6 - b7 >= -4', 'a0 - a1 + a3 + a5 + b1 - b5 + b6 - b7 >= -2', 'a0 - a3 + a5 + b1 - b4 - b5 + b6 + b7 >= -2', '- a0 + a3 + a4 + a6 + a7 + b5 + b6 + b7 >= 0', '- a0 + a6 + a7 + b0 + b1 + b3 - b5 - b6 >= -2', '- a4 - a5 - b2 - b4 - b5 >= -4', '- a0 + a4 + a5 + b1 + b2 + b5 - b6 >= -1', '- a1 + a6 - b2 + b4 + b5 - b6 - b7 >= -3', '- a1 - a3 - b0 - b1 + b5 - b7 >= -4', '- a0 - a4 - b0 + b4 - b5 - b7 >= -4', 'a4 + a6 - a7 + b1 - b4 + b5 >= -1', '- a0 + a1 - a3 - b0 - b1 + b6 - b7 >= -4', '- a1 + a3 + a5 - b1 - b2 - b4 + b5 + b7 >= -3', '- a1 - a3 + a5 - b1 - b2 + b5 - b7 >= -4', '- a1 + a3 + a5 - b0 - b2 + b5 - b7 >= -3', 'a1 - a5 - a6 - a7 + b2 + b4 - b5 + b7 >= -3', '- a1 - a3 + a5 - a7 - b0 + b1 - b4 + b5 + b7 >= -4', '- a1 + a4 + a5 + b1 + b2 - b5 + b6 >= -1', 'b0 + b1 + b2 + b4 + b5 >= 1', 'a0 + a1 + b2 + b4 + b6 + b7 >= 1', '- a0 + a3 + a5 + b1 - b2 - b4 + b6 + b7 >= -2', 'a0 - a3 - a4 - b0 - b2 - b4 + b7 >= -4', '- a0 - a1 - a3 + a5 + b1 - b2 + b6 - b7 >= -4', '- a1 + a2 + a3 - a5 + b2 - b4 + b5 - b7 >= -3', '- a1 + a2 - a3 - a5 + b2 - b4 + b5 + b7 >= -3', '- a0 + a1 - a5 - b0 - b4 - b5 >= -4', '- a1 - a5 + a6 + b1 - b3 - b5 - b6 + b7 >= -4', 'a0 + a1 + b1 + b4 + b6 + b7 >= 1', '- a1 + b0 + b1 + b2 + b3 - b5 + b6 >= -1', '- a5 + b0 + b1 + b2 - b3 + b5 + b6 >= -1', 'a0 - a1 + a4 + a6 + a7 + b0 + b3 - b6 >= -1', '- a0 + a6 - a7 - b0 + b2 - b4 >= -3', 'a4 + a5 + a6 + b1 + b5 - b6 + b7 >= 0', 'a3 - b0 - b1 + b5 - b6 + b7 >= -2', '- a1 + a5 - a7 + b0 + b3 - b4 + b5 - b6 >= -3', 'a4 - a5 + a6 + b1 - b2 + b4 + b5 - b6 - b7 >= -3', 'a1 + a4 + b1 + b4 + b5 + b6 - b7 >= 0', '- a2 - a5 - b0 - b1 + b5 + b7 >= -3', 'a1 - a4 - a6 + a7 - b2 + b5 - b6 - b7 >= -4', '- a0 + a1 - a3 - a7 - b2 + b4 + b5 + b6 + b7 >= -3', 'a4 + a6 + a7 + b1 + b4 + b5 - b6 - b7 >= -1', '- a1 + a5 - a7 + b0 - b1 - b3 - b4 + b5 + b6 >= -4', '- a1 - a5 - a6 - a7 + b1 + b4 - b6 - b7 >= -5', '- a1 - a5 - b2 - b5 - b6 + b7 >= -4', '- a4 - a5 - a6 + a7 - b2 - b5 + b7 >= -4', '- a0 + a1 + a5 + a7 - b0 - b2 - b4 >= -3', 'a4 + a5 + b2 + b4 + b5 - b6 + b7 >= 0', '- a1 + a7 + b2 + b4 + b5 - b6 - b7 >= -2', '- a1 + a2 + b2 + b4 + b5 - b6 - b7 >= -2', 'a1 - a2 - a7 + b0 + b2 + b4 + b5 - b7 >= -2', '- a1 - a5 - b0 - b4 - b5 - b6 >= -5', '- a1 - a2 + a3 - a5 - a7 + b2 + b5 - b6 + b7 >= -4', '- a3 - a5 + a7 - b1 - b4 + b5 - b6 + b7 >= -4', 'a5 + a6 - a7 - b0 - b1 - b5 - b6 >= -4', '- a0 + a1 - a3 + a4 + a5 + b2 - b4 + b5 + b6 >= -2', 'a4 + a5 + a6 + b2 - b5 + b6 >= 0', 'a4 + a6 - b1 - b2 + b4 + b5 + b6 + b7 >= -1', 'a1 + a5 + b0 + b1 + b2 + b5 >= 1', 'a4 + a5 + b1 + b2 + b4 >= 1', 'a7 + b1 + b2 + b4 + b5 + b6 + b7 >= 1', 'a2 - a7 - b0 + b2 + b4 + b5 - b6 - b7 >= -3', '- a1 + a2 - a5 - b1 - b5 - b6 + b7 >= -4', '- a0 - a3 + a5 - b1 - b2 - b4 + b6 + b7 >= -4', '- a0 + a1 - a2 + a6 - a7 + b2 - b5 + b7 >= -3', 'a0 + a5 + b0 + b2 - b5 + b6 >= 0', 'a0 + a1 + a4 + a5 + b1 + b2 >= 1', 'a4 + a5 + b0 + b1 + b2 >= 1', 'a0 + a1 + b0 + b1 + b2 + b3 >= 1', 'a0 - a1 - b1 + b5 - b6 >= -2', '- a2 - a3 - a5 - a7 + b2 - b4 + b5 - b6 - b7 >= -6', '- a4 - b1 - b2 - b3 - b6 >= -4', 'a0 - a1 - a5 - b4 + b5 - b6 >= -3', 'a3 - a5 + a7 - b1 - b4 - b6 - b7 >= -4', '- a1 - a2 - a5 - a7 + b2 + b4 + b5 + b7 >= -3', 'a2 + a4 + a5 + b2 + b4 + b5 - b7 >= 0', 'a0 - a1 + a5 + a7 + b1 + b2 + b5 + b6 >= 0', 'a1 + a4 + b0 + b1 + b2 + b5 >= 1', '- a3 - a4 + a7 - b0 + b1 + b2 - b6 >= -3', '- a1 + a3 + a5 - a7 - b0 + b1 + b5 - b7 >= -3', '- a5 + a6 + b0 - b1 - b2 - b3 + b4 - b5 + b7 >= -4', 'a1 + b0 + b1 + b4 + b5 + b6 >= 1', 'a0 + a1 + a4 + a7 + b1 + b2 + b4 >= 1', 'a1 + a2 - a5 + b0 + b2 + b4 + b5 + b7 >= 0', '- a1 + a2 + a4 + a5 + b2 - b3 + b5 - b6 >= -2', '- a0 - a2 - a7 + b1 + b2 - b5 - b6 >= -4', '- a1 - a5 - a6 + a7 - b2 - b3 + b7 >= -4', '- a6 - a7 - b0 + b2 + b4 + b5 - b6 - b7 >= -4', '- a4 - a7 + b1 - b2 + b4 + b5 + b6 + b7 >= -2', '- a0 - a1 - a3 - a4 - b0 + b1 + b6 - b7 >= -5', '- a3 + a5 + b0 + b1 + b2 + b5 - b7 >= -1', 'a4 + b0 + b1 + b2 + b4 + b6 >= 1', '- a3 - a5 - b1 - b2 - b4 + b5 - b6 + b7 >= -5', '- a0 + a3 + a5 - a7 - b0 + b1 + b2 - b7 >= -3', 'a0 + a4 + a5 + b2 + b4 + b5 >= 1', '- a0 + a1 + a2 - b1 + b4 - b5 - b6 - b7 >= -4', 'a2 - a5 - a7 - b0 + b2 - b5 + b7 >= -3']
        self.q5_4150_ineqs = ['- a3 - a6 - b7 >= -2', 'a2 + a3 + b3 - b6 + b7 >= 0', 'a3 - a6 + b7 >= 0', 'b2 >= 1', 'b4 >= 1', 'a1 - b5 >= 0', 'a0 - b6 >= 0', 'a4 - a5 >= 0', 'b0 - b1 >= 0', 'a5 + b1 >= 1', '- a1 - b0 + b5 >= -1', '- a4 + b3 >= 0', 'a1 - a4 >= 0', '- a4 - b5 >= -1', '- a6 - b1 >= -1', '- a3 + b5 + b6 - b7 >= -1', 'a7 - b1 >= 0', 'a0 - b1 >= 0', '- a3 - b1 + b6 + b7 >= -1', '- a0 + b1 + b6 >= 0', 'a1 + a3 - b6 >= 0', 'a3 - b5 + b6 - b7 >= -1', 'a3 - a7 + b1 + b7 >= 0', '- a2 - a3 - b3 - b5 - b6 - b7 >= -5', 'a2 + a3 - b3 - b5 - b7 >= -2', '- a2 + a3 - b3 - b5 - b6 + b7 >= -3', '- a1 - a3 - a7 + b5 - b7 >= -3', 'a2 - a3 - b3 - b5 + b7 >= -2', '- a2 + a3 + b3 - b6 - b7 >= -2', 'a2 - a3 + b3 - b5 - b6 - b7 >= -3', '- a2 - a3 + b3 - b5 + b7 >= -2', 'a3 + a6 + a7 - b6 - b7 >= -1', '- a3 + a6 + a7 - b6 + b7 >= -1', 'a3 + b1 + b6 + b7 >= 1']
        self.q6_ineqs = ['a0 + a5 - b0 + b1 + b6 >= 0', 'a4 + b0 + b1 - b3 >= 0', '- a4 + b0 + b1 + b3 >= 0', 'b0 + b1 + b4 + b5 >= 1', 'a4 - a6 + b2 + b4 >= 0', 'a4 + a6 + b2 - b4 >= 0', '- a1 + a5 + a6 + b5 >= 0', 'a1 + a5 + a6 - b5 >= 0', 'a0 + b0 + b5 - b6 >= 0', '- a5 - a6 - a7 + b2 + b4 + b5 + b7 >= -2', 'b0 + b2 + b4 >= 1', 'a0 + a1 + a3 - b1 >= 0', 'a0 + a1 - a3 + b1 >= 0', 'a1 + a5 - a6 + a7 + b5 >= 0', '- a1 + a5 - a6 + a7 - b2 - b5 >= -3', 'a2 + a4 - b0 + b1 + b3 >= 0', '- a1 + a5 - a6 - a7 + b2 - b5 >= -3', 'a2 - a4 - b0 + b1 - b3 >= -2', 'a4 - a6 + a7 - b2 - b4 >= -2', '- a5 + b4 + b6 + b7 >= 0', 'a4 + a5 + a7 - b2 >= 0', '- a0 + b0 + b5 + b6 >= 0', 'a2 - b0 + b1 + b2 >= 0', '- a4 - a5 - b0 - b1 - b2 - b4 >= -5', 'a5 + b4 + b6 - b7 >= 0', 'a4 + a5 - a7 + b2 >= 0', 'a1 + a4 - a5 + b1 + b5 >= 0', 'a5 + b0 + b2 >= 1', '- a4 + a5 + a7 + b2 + b5 >= 0', 'a4 + b0 - b1 + b3 + b6 >= 0', 'a4 - a5 + a6 - a7 - b4 >= -2', '- a0 + a1 - b0 - b1 - b2 - b4 - b5 >= -5', 'a1 - a4 + b1 + b4 + b5 >= 0', 'a4 - a5 - a6 - a7 + b4 >= -2', '- a1 - a4 - b0 - b1 - b2 + b4 - b7 >= -5', '- a1 - a6 - b0 - b1 - b2 + b4 - b6 >= -5', 'a0 + b0 + b3 - b5 + b6 >= 0', '- a0 + a1 - a4 - a5 - b2 - b4 - b5 >= -5', '- a2 - a4 + b1 - b2 + b3 >= -2', 'a4 + a6 - a7 + b0 - b6 >= -1', 'a0 + a1 + b2 + b4 >= 1', '- a2 + a4 + b1 - b2 - b3 >= -2', 'a4 + b1 + b4 + b6 >= 1', 'a0 + a1 + b0 + b2 >= 1', '- a4 + a5 - a7 + b1 - b2 + b6 >= -2', 'b2 + b4 + b6 + b7 >= 1', 'a4 + a6 + a7 - b2 + b4 >= 0', '- a4 + a5 - a7 + b0 + b6 >= -1', '- a1 - a5 - b0 - b1 + b2 - b4 - b5 - b6 >= -6', '- a2 + a6 - b0 - b1 + b2 - b4 >= -3', '- a0 - a4 + b0 + b3 - b5 - b6 >= -3', 'a0 + a5 + a6 + b0 + b6 >= 1', 'a2 - a4 - a6 - b0 + b2 - b4 >= -3', '- a0 + a1 - a5 - b0 - b1 - b4 - b5 >= -5', '- a0 + a4 + b0 - b3 - b5 - b6 >= -3', 'a1 + a3 - a5 - b0 - b1 + b5 - b6 >= -3', '- a1 - a4 - a5 - b0 - b4 - b5 - b6 >= -6', 'a2 - a4 + a6 + b2 + b4 >= 0', '- a0 - a1 - a4 - a5 - b4 + b5 + b6 >= -4', '- a2 - a6 - b1 + b2 + b4 >= -2', '- a4 + a5 - a7 + b1 - b2 + b5 >= -2', 'a0 - a1 - a6 - b0 - b1 - b2 - b6 >= -5', 'a0 + a1 + a5 + b2 >= 1', 'a4 - a5 + a6 + a7 + b1 + b4 >= 0', 'a0 + a1 + a4 + a6 + b4 >= 1', '- a0 - a4 - a5 - b1 - b2 - b3 - b4 >= -6', 'a0 - a4 + b0 - b1 - b3 - b6 >= -3', 'b0 + b4 + b6 + b7 >= 1', 'a0 + a1 + b0 + b4 >= 1', 'a1 - a4 + a5 - a7 - b2 - b5 >= -3', '- a0 + a1 - a3 - a5 - b1 + b5 + b6 >= -3', '- a0 - a1 - a3 - a5 - b1 - b4 + b6 - b7 >= -6', 'b1 + b4 + b6 + b7 >= 1', '- a1 + a5 + a7 + b2 - b5 + b6 >= -1', '- a1 - a5 - a7 - b2 + b4 + b5 - b6 - b7 >= -5', '- a2 - a4 - a6 + b0 + b2 >= -2', '- a0 - a1 + a3 - a5 - b1 + b6 + b7 >= -3', '- a1 - a5 - a7 - b2 + b4 - b5 + b7 >= -4', 'a0 + a4 + b0 - b1 + b3 >= 0', 'a5 + a7 + b0 - b1 - b3 + b6 >= -1', '- a0 - a7 - b0 - b1 - b2 - b4 - b5 - b6 >= -7', 'a1 - a5 - a7 - b2 + b4 - b5 - b6 - b7 >= -5', 'a0 + a1 + b1 + b4 >= 1', 'a0 + a1 + a5 + b1 >= 1', '- a1 - a4 + a6 + b4 + b5 - b6 - b7 >= -3', '- a1 - a4 - a5 + a6 + b4 - b5 + b7 >= -3', 'a1 - a5 - a7 + b4 + b5 + b7 >= -1', 'a1 - a4 + a6 + b4 - b5 - b6 - b7 >= -3', 'a1 - a4 + a5 + a7 + b2 >= 0', '- a1 + a4 - a5 + b1 - b5 + b6 >= -2', '- a0 - a2 + a5 - b1 + b2 + b5 - b6 >= -3', '- a1 - a4 + a5 - a7 - b2 + b5 >= -3', 'a4 + a6 + b1 - b2 + b4 >= 0', '- a1 - a4 - b0 - b1 - b2 + b4 - b6 >= -5', '- a4 + b0 - b1 + b2 - b3 + b6 >= -2', 'a4 + a5 + a6 + b1 - b5 - b6 >= -1', 'a2 + a6 + b0 - b1 + b2 >= 0', 'a0 - a1 - a5 - b1 + b2 - b4 + b5 - b6 >= -4', '- a0 - a1 + a4 - a6 + b0 - b3 - b4 >= -4', '- a0 + a5 + a7 - b0 - b1 - b2 - b4 + b6 >= -4', 'a1 - a5 + b2 + b4 + b5 >= 0', 'a0 + a5 + b2 + b6 >= 1', 'a3 + a5 - b0 - b1 - b2 + b5 - b6 + b7 >= -3', 'a1 + a3 + a4 + b1 + b5 - b6 >= 0', 'a0 - a1 - a4 - a5 - b4 + b5 - b6 >= -4', 'a0 + a1 + a4 + b1 >= 1', '- a5 + a7 + b2 + b4 + b5 - b6 - b7 >= -2', '- a1 - a5 + a7 + b2 + b4 - b5 + b7 >= -2', '- a1 - a5 + b0 + b1 - b2 - b3 - b4 - b6 >= -5', 'a0 - a1 - a5 + b0 - b1 - b2 + b3 - b4 >= -4', '- a0 - a1 - a4 + b0 + b2 + b3 >= -2', 'a1 - a5 + a6 + b4 + b5 + b7 >= 0', 'b0 + b1 + b2 + b3 >= 1', 'a0 + a3 + a4 - a5 - b0 - b4 + b6 - b7 >= -3', 'a0 - a3 - b0 + b1 - b4 + b6 - b7 >= -3', 'a0 - a1 - a3 + a4 - a5 - b0 + b6 + b7 >= -3', 'a0 + a3 - b0 - b1 + b2 - b4 + b6 - b7 >= -3', 'a0 + a3 + a5 - b0 - b4 + b6 + b7 >= -1', 'a0 - a1 + a3 - b0 + b1 + b6 + b7 >= -1', 'a1 - a3 - a5 - b0 + b1 + b5 - b6 >= -3', 'a1 + a5 - a6 + b2 + b5 >= 0', 'a1 + a4 - a5 + a6 + b4 + b5 >= 0', '- a1 - a3 - a5 + a7 - b0 - b1 + b2 - b4 + b5 + b7 >= -5', 'a0 + a1 + a4 + b0 >= 1', '- a1 - a4 - a6 + a7 - b2 - b5 - b6 - b7 >= -6', '- a1 - a3 - a5 - a6 + b0 + b1 + b5 + b7 >= -3', '- a1 - a3 + a5 - b0 - b1 - b2 + b5 - b7 >= -5', 'a1 + a7 + b2 + b4 - b5 - b6 - b7 >= -2', 'a1 + a3 + a4 + b4 + b5 + b6 >= 1', '- a0 - a1 + a4 - a7 - b1 - b4 + b5 + b6 >= -4', '- a1 - a4 - a5 - a6 + a7 - b2 + b4 + b5 + b7 >= -4', 'a4 + a6 - a7 - b1 - b4 + b5 >= -2', '- a0 + a1 + a3 - a5 + b1 + b5 + b6 >= -1', 'a1 + a5 - a6 + b1 + b5 >= 0', '- a0 - a4 + a5 - a6 - b0 + b1 - b4 + b5 - b6 >= -5', '- a0 - a1 + a4 + a7 - b0 - b1 - b2 - b4 + b5 >= -5', '- a0 + b0 - b1 - b3 + b4 + b6 >= -2', 'a4 + a5 + a6 - b0 - b1 - b5 + b6 >= -2', 'a1 - a4 - a5 - a6 + a7 - b2 + b4 - b5 + b7 >= -4', '- a0 + a1 - a4 - a5 - b0 - b4 - b5 >= -5', 'a0 - a1 - a3 - a5 - b1 - b5 + b6 + b7 >= -4', '- a0 - a1 + a3 - a5 + b1 - b4 - b5 + b6 - b7 >= -5', '- a0 + a1 - a3 + a5 - a7 - b2 - b4 + b5 + b6 >= -4', '- a2 + b0 + b1 + b2 >= 0', '- a0 - a1 - a3 - a5 + b1 - b5 + b6 + b7 >= -4', 'a5 + a7 + b1 + b2 + b6 >= 1', '- a1 + a4 + a7 + b0 - b5 + b6 >= -1', 'a0 + a1 + a4 + b4 + b6 >= 1', 'a5 - a6 + a7 + b0 + b3 - b4 - b6 >= -2', 'a0 - a1 - a4 - b0 - b1 - b2 - b6 >= -5', '- a0 + a1 - a4 - a5 - b1 - b2 - b4 >= -5', 'a0 + a1 + b4 + b6 + b7 >= 1', '- a1 + a3 - a5 + a7 - b0 - b1 + b2 - b4 + b5 - b7 >= -5', 'a1 - a4 - a6 + a7 + b4 + b5 - b6 - b7 >= -3', '- a5 + b0 + b1 - b2 - b3 + b5 + b6 >= -2', '- a0 - a3 + a7 - b0 + b1 + b2 + b5 - b6 - b7 >= -4', '- a1 - a3 + a6 + a7 + b0 + b1 + b5 + b7 >= -1', '- a0 + a4 - a6 - b0 + b1 - b2 - b4 + b5 - b6 >= -5', 'a2 + a3 - a5 + a6 - a7 + b2 + b5 - b6 + b7 >= -2', '- a3 + a5 + b0 + b1 + b5 - b7 >= -1', '- a1 - a6 - a7 + b2 + b4 - b5 - b6 - b7 >= -5', 'a0 - a1 - a3 + a5 - b5 + b6 - b7 >= -3', '- a2 + a3 - a5 - b0 - b1 + b2 + b5 - b6 - b7 >= -5', 'a3 + a5 + b0 + b1 + b5 + b7 >= 1', '- a0 - a1 - a3 + a5 - a7 - b0 - b2 - b4 - b5 + b6 + b7 >= -7', 'a4 + a5 + a6 + b2 - b5 >= 0', '- a0 + a5 + a6 + a7 - b0 - b1 - b4 + b5 >= -3', 'a0 + a5 + a7 + b1 + b2 + b5 >= 1', 'a4 + b0 + b4 + b6 >= 1', '- a3 - a4 - a5 - a6 - b0 + b1 - b4 - b6 - b7 >= -7', '- a1 + a5 - a6 + b1 - b5 + b6 >= -2', 'a3 - a4 - a7 - b0 - b1 - b2 + b5 - b6 >= -5', 'a3 - a5 + b0 + b1 + b3 + b5 - b7 >= -1', '- a0 + a3 + a5 - b0 - b1 - b2 - b5 + b6 - b7 >= -5', '- a1 + a3 - a4 - a5 - a6 - b0 + b1 - b4 - b6 + b7 >= -6', '- a2 + a3 + a4 + a7 - b0 - b2 - b3 - b5 - b6 - b7 >= -6', 'a3 - a4 - a5 + b0 - b3 - b4 + b5 + b6 - b7 >= -4', '- a4 - a7 + b4 + b5 + b6 + b7 >= -1', 'a2 + a3 + a4 + a7 - b0 - b2 + b3 - b5 - b6 - b7 >= -4', 'a4 + a5 + b1 + b2 - b5 >= 0', 'a2 + a3 + a4 + a7 - b0 - b1 - b2 - b3 - b5 - b6 + b7 >= -5', 'a1 - a5 - a6 - a7 + b2 + b4 + b7 >= -2', '- a1 + b0 + b2 - b5 + b6 >= -1', '- a0 + a1 - a3 + a5 + a7 - b0 + b4 + b5 - b7 >= -3', '- a1 + a4 - a6 + b0 - b5 + b6 >= -2', '- a3 + a4 + a5 + a6 - b1 + b5 - b6 + b7 >= -2', '- a1 - a2 - a3 - a5 - a6 - b0 - b1 + b5 + b7 >= -6', '- a0 - a1 + a3 + a7 - b0 + b1 + b2 - b4 + b5 - b6 + b7 >= -4', '- a0 - a2 - a3 + a4 + a7 - b0 - b1 - b2 + b3 - b5 - b6 - b7 >= -8', 'a1 + b0 + b1 + b2 + b5 >= 1', '- a1 - a2 - a3 + a7 - b0 - b1 - b2 - b3 - b5 - b6 + b7 >= -8', '- a1 + a2 - a3 + a7 - b0 - b1 - b2 + b3 - b5 - b6 + b7 >= -6', '- a0 + a2 - a3 + a5 + a7 - b0 - b2 - b3 - b7 >= -5', '- a2 + a3 + a7 - b0 - b1 - b2 + b3 - b4 - b5 - b6 + b7 >= -6', '- a0 - a2 + a3 + a5 + a7 - b0 - b1 - b2 - b3 - b5 - b7 >= -7', '- a0 - a2 - a3 + a5 + a7 - b0 - b2 + b3 - b7 >= -5', '- a1 + a2 - a7 + b2 + b4 - b5 - b6 - b7 >= -4', 'a1 + a3 - a6 - b0 - b1 + b5 - b6 >= -3', 'a2 + a3 + a5 + a7 - b0 - b1 - b2 + b3 - b5 - b6 - b7 >= -5', 'a2 + a3 + a5 + a7 - b0 - b2 - b3 - b6 + b7 >= -3', '- a1 - a2 - a5 - b0 - b1 + b4 - b5 + b7 >= -5', 'a3 - a5 - a7 - b0 - b1 - b2 + b5 - b6 - b7 >= -6', '- a0 + a4 - a5 + a7 + b1 + b2 + b5 - b6 >= -2', '- a1 - a3 - a5 - a6 - b0 - b1 - b2 + b5 + b7 >= -6', 'a7 + b0 + b1 + b3 + b5 - b6 >= 0', '- a1 + a3 - a5 + b0 - b1 - b3 - b4 - b5 + b6 - b7 >= -6', '- a1 - a3 - a4 - a5 - b1 - b3 + b5 + b6 + b7 >= -5', '- a1 + a2 - a3 + a6 - a7 + b2 - b4 + b5 - b6 - b7 >= -5', '- a1 + a3 - a4 + a6 + a7 + b1 - b2 + b5 - b6 - b7 >= -4', '- a3 - a4 - a7 - b0 + b1 - b2 + b5 - b6 - b7 >= -6', '- a3 - a4 - a5 + a6 + a7 - b0 - b2 - b4 - b6 + b7 >= -6', '- a0 - a2 - a3 - a6 - a7 - b0 - b1 - b3 + b4 - b5 - b7 >= -9', 'a2 + a3 - a6 - a7 - b0 - b1 - b2 - b3 + b4 - b5 - b7 >= -7', '- a0 + a2 - a3 - a6 - a7 - b0 - b1 - b2 + b3 + b4 - b5 - b7 >= -8', '- a1 + a3 + a6 - a7 + b1 - b2 - b4 + b5 + b7 >= -3', '- a2 + a3 - a6 - a7 - b0 - b1 + b3 + b4 - b5 - b7 >= -6', 'a3 - a4 - a5 - a7 - b1 - b2 - b4 - b6 - b7 >= -7', 'a5 + a6 + b0 + b1 + b5 >= 1', '- a1 - a3 - a4 - a7 - b1 - b2 - b4 + b5 + b7 >= -6', 'a0 - a3 + b1 - b4 - b5 + b6 - b7 >= -3', '- a1 + a3 + a4 + b1 - b5 + b6 - b7 >= -2', 'a0 - a1 + a3 + b1 - b5 + b6 + b7 >= -1', 'a1 - a2 - b1 + b2 + b4 - b6 - b7 >= -3', 'a1 + a2 + a4 - a7 + b2 + b4 + b7 >= 0', '- a0 - a2 - a3 - a6 + a7 - b0 - b1 - b3 + b4 - b5 - b6 + b7 >= -8', '- a2 - b1 + b2 + b4 + b5 - b6 - b7 >= -3', '- a1 + a2 - a3 + a6 + b2 + b5 + b6 + b7 >= -1', '- a0 + a2 - a3 - a6 - a7 - b0 - b1 - b2 - b3 - b5 - b6 + b7 >= -9', '- a0 + a2 - a3 - a6 + a7 - b0 - b1 - b2 + b3 - b5 - b6 + b7 >= -7', '- a2 + a3 - a6 - a7 - b0 - b1 - b3 + b4 - b6 + b7 >= -6', '- a2 + a3 - a6 + a7 - b0 - b1 + b3 + b4 - b6 + b7 >= -4', '- a0 - a2 - a3 - a6 - a7 - b0 - b1 + b3 + b4 - b5 - b6 + b7 >= -8', 'a2 + a6 - a7 + b2 + b4 + b5 + b7 >= 0', 'a2 + a3 - a6 - a7 - b0 - b1 - b2 + b3 - b5 - b6 + b7 >= -6', 'a1 + a5 - a6 + b4 + b5 + b6 >= 0', '- a1 + a2 - a3 - a5 - b0 - b1 - b3 - b4 - b5 - b6 - b7 >= -9', '- a1 - a3 - a4 + a6 + a7 + b0 + b3 + b5 - b6 - b7 >= -4', 'a3 - a4 - a5 + a6 + a7 - b1 - b2 + b3 - b4 + b7 >= -4', 'a3 - a6 + a7 + b0 - b2 + b3 - b4 - b6 - b7 >= -4', '- a3 - a5 - a6 + a7 - b1 - b2 - b4 + b5 - b6 + b7 >= -6', '- a1 + a2 - a3 + a4 - a7 - b0 + b2 - b4 - b6 - b7 >= -6', 'a2 + a3 + a4 - a7 - b0 + b2 - b4 - b6 + b7 >= -3', '- a1 + a3 + a6 - a7 - b0 + b1 + b2 - b4 + b5 - b6 - b7 >= -5', '- a3 - a5 + a6 - a7 - b0 + b1 + b2 - b4 - b6 + b7 >= -5', '- a3 + a4 + a5 + a6 - b1 + b6 - b7 >= -2', '- a0 - a2 - a3 - a4 - a5 - a6 - b0 - b1 - b3 - b5 - b6 - b7 >= -11', '- a0 + a2 - a3 + a6 + a7 - b0 - b2 - b3 + b4 - b5 - b6 >= -6', '- a2 + a3 - a4 - a5 + a7 - b0 - b1 - b2 - b3 + b7 >= -6', '- a0 - a2 - a3 - a4 - a5 + a6 - a7 - b0 - b1 - b2 - b3 + b7 >= -9', '- a0 - a2 - a3 - a4 - a5 + a7 - b2 + b3 - b5 - b6 + b7 >= -7', 'a2 + a3 - a4 - a5 - b0 - b2 - b3 - b6 - b7 >= -6', '- a0 + a2 - a3 - a4 - a5 - b1 - b2 + b3 - b5 - b6 - b7 >= -8', 'a3 + a5 + a6 - a7 + b0 - b3 + b6 + b7 >= -1', '- a2 + a3 - a4 - a5 - b0 - b2 + b3 - b6 - b7 >= -6', 'a2 + a3 - a4 + a6 - a7 - b0 - b3 + b4 - b6 >= -4', 'a2 + a3 - a4 - a5 + a7 - b0 - b1 - b2 + b3 + b7 >= -4', '- a0 + a2 - a3 - a4 + a6 - a7 - b1 + b3 + b4 - b5 - b6 >= -6', '- a2 + a3 - a4 + a6 - a7 - b0 - b2 + b3 - b6 >= -5', 'a3 + a7 + b0 + b1 + b2 + b5 - b7 >= 0', '- a0 + a2 - a3 + a4 + a7 - b1 - b2 - b3 - b5 - b6 - b7 >= -7', '- a1 - a2 - a3 - a5 - b0 - b1 + b3 - b4 - b5 - b6 - b7 >= -9', 'a0 - a1 - a2 - a3 - b1 - b2 - b3 + b5 - b6 - b7 >= -7', '- a0 - a2 - a3 + a6 + a7 - b1 - b2 - b3 + b4 + b6 - b7 >= -6', 'a0 + a2 + a3 - b1 - b2 - b3 + b5 - b6 - b7 >= -4', 'a1 + a2 + a3 - a4 + a6 + a7 - b0 - b3 - b5 - b7 >= -4', 'a0 - a1 + a2 - a3 + a6 - b2 + b3 + b5 - b6 - b7 >= -4', 'a0 - a2 + a3 - a5 - b1 + b3 + b5 - b6 - b7 >= -4', '- a0 + a2 - a3 - a4 - a5 + a6 - a7 - b0 - b3 + b6 - b7 >= -7', '- a0 + a2 - a3 - a4 + a6 + a7 - b0 - b1 + b3 + b6 - b7 >= -5', 'a1 - a2 + a3 - a4 + a6 - a7 - b1 - b2 - b3 - b5 - b7 >= -7', '- a2 + a3 + a6 + a7 - b0 - b2 + b3 + b4 - b5 - b7 >= -4', '- a1 + a2 - a3 + a6 + a7 - b0 - b1 - b3 - b4 + b5 + b7 >= -5', '- a0 - a2 - a3 - a4 - a5 + a6 - a7 - b0 - b2 + b3 + b6 >= -7', 'a0 - a2 + a3 - b1 - b2 - b3 - b4 + b5 - b6 + b7 >= -5', '- a0 + a2 - a3 - a4 - a5 - a6 + a7 - b1 - b2 - b3 + b6 >= -7', 'a1 - a2 + a3 - a4 - a6 + a7 - b1 - b3 - b5 + b6 - b7 >= -6', '- a1 - a2 - a3 + a6 - b0 - b1 + b3 - b4 + b5 + b7 >= -5', '- a0 - a2 - a3 - a4 - a6 + a7 - b0 - b2 + b3 + b6 - b7 >= -7', 'a1 + a2 + a3 - a4 + a6 - a7 - b0 - b1 + b3 - b5 - b7 >= -5', 'a0 + a2 + a3 - b1 - b2 + b3 - b4 + b5 - b6 + b7 >= -3', 'a2 + a3 - a4 - a6 + a7 - b0 - b1 - b2 + b3 - b5 + b6 - b7 >= -6', '- a3 + a4 + a6 - a7 + b3 - b5 + b6 + b7 >= -2', 'a3 - a6 + b0 + b1 + b2 + b5 - b7 >= -1']
        self.q7_ineqs = ['a4 + a5 + a7 >= 1', 'a5 + b4 + b6 >= 1', 'a4 + a6 + b2 >= 1', 'b0 + b1 + b3 >= 1', 'a4 + a6 + b4 >= 1', 'a1 + b1 + b5 >= 1', 'a1 + b2 + b5 >= 1', 'a1 + a5 + a6 - b5 >= 0', 'a2 + b1 + b2 >= 1', 'a0 + a5 + b6 >= 1', 'a0 - b0 + b1 + b6 >= 0', 'b0 - b1 - b3 + b6 >= -1', 'a0 + b2 + b6 >= 1', 'a5 + b0 >= 1', 'b0 + b4 >= 1', 'a0 + a4 + b6 >= 1', 'b0 + b2 >= 1', '- a1 + a5 + a6 + b5 >= 0', 'a4 + b0 >= 1', 'a5 + b1 >= 1', 'b1 + b4 >= 1', 'a1 + a5 - a6 + b5 >= 0', 'a0 + b0 + b5 - b6 >= 0', 'a0 + a1 + b1 >= 1', 'a5 + b2 >= 1', 'b2 + b4 >= 1', 'a0 + a4 + a6 >= 1', 'a0 + a1 + b2 >= 1', 'a4 - a5 - a7 + b4 >= -1', 'a0 + a1 + a4 >= 1', 'a4 + b1 >= 1', 'a6 + b4 + b5 - b6 - b7 >= -1', '- a0 + a1 - a3 - a5 + b5 + b6 >= -2', 'a1 + a2 + a3 - a4 + a6 + a7 - b0 - b1 - b2 + b3 - b5 + b6 >= -4', 'a1 - a4 + a5 - a7 >= -1', 'a1 + b4 >= 1', 'b4 + b6 + b7 >= 1', '- a4 + a5 - a7 + b6 >= -1', '- a1 + b1 - b5 + b6 >= -1', 'a0 + b0 - b1 - b3 >= -1', '- a1 + b2 - b5 + b6 >= -1', '- a4 - a5 - a6 + a7 + b4 + b5 + b7 >= -2', 'a1 - a2 + a3 - a4 + a6 + a7 - b1 - b3 - b5 + b6 >= -4', '- a0 - a1 + a2 + a3 + a6 - a7 - b1 - b2 + b3 + b6 + b7 >= -4', 'a0 + a3 - b0 - b4 + b6 - b7 >= -2', 'a0 - a1 - a3 - b0 + b6 + b7 >= -2', '- a5 + a6 + b4 - b5 + b7 >= -1', 'a1 + a2 + a3 - a4 - a6 - a7 - b0 - b1 + b3 - b5 + b6 >= -5', '- a4 + a5 - a7 + b5 >= -1', '- a0 - b0 + b1 + b5 - b6 >= -2', 'a1 - a2 + a3 - a5 - a6 + a7 - b0 - b2 + b3 - b5 + b6 >= -5', '- a0 + b2 + b5 - b6 >= -1', '- a2 + a6 - b1 + b2 >= -1', 'a0 + a2 + a3 + a6 - b2 - b3 - b4 + b5 - b6 + b7 >= -3', '- a0 - a1 - a2 + a3 + a6 - a7 - b1 - b3 + b6 + b7 >= -5', '- a1 - a2 + a3 - a4 + a6 + a7 - b1 - b3 - b4 + b5 + b6 - b7 >= -6', '- a0 + a1 - a2 - a3 - a5 - a6 + a7 - b1 - b2 - b3 + b6 >= -7', '- a0 + a1 - a2 - a3 - a4 - a6 - a7 - b0 - b2 + b3 + b6 >= -7', '- a5 - a7 + b4 - b5 + b7 >= -2', 'a2 - b0 + b1 - b3 >= -1', 'a1 - a2 + a3 - a4 - a6 - a7 - b1 - b2 - b3 - b5 + b6 >= -7', '- a2 + b1 - b2 + b3 >= -1', 'a0 + a1 + a3 >= 1', 'a2 - a4 - a6 + b2 >= -1', 'a0 - a1 + a2 - a3 + a6 - b2 - b3 + b5 - b6 - b7 >= -5', 'a0 + a2 + a3 - a4 - a6 - a7 - b1 + b3 - b4 + b5 - b6 + b7 >= -5', '- a0 - a1 + a2 + a3 - a4 - a6 - a7 - b3 + b6 + b7 >= -5', '- a0 + a1 + a2 - a3 - a4 + a6 + a7 - b1 - b2 + b3 - b5 - b6 >= -6', '- a0 - a1 - a2 + a3 - a4 - a6 - a7 - b2 + b3 + b6 + b7 >= -6', '- a4 - a6 + a7 + b4 - b5 - b6 - b7 >= -4', 'a1 + a3 - a5 - b0 + b5 - b6 >= -2', '- a0 - a1 - a2 + a3 - a5 - a6 + a7 - b1 - b2 - b3 + b6 + b7 >= -7', 'a0 - a1 + a2 - a3 - a4 - a6 - a7 - b1 + b3 + b5 - b6 - b7 >= -7', '- a0 + a1 + a2 - a3 + a6 - a7 - b0 - b2 - b3 - b5 - b6 >= -7', 'a0 - a1 - a2 - a3 - a5 + a7 - b2 + b3 - b4 + b5 - b6 - b7 >= -7', '- a0 - a2 - a3 - a5 - a6 + a7 - b1 - b2 - b3 - b4 + b6 - b7 >= -9', 'a0 - a2 + a3 - a4 - a6 - a7 - b1 - b2 - b3 - b4 + b5 - b6 + b7 >= -8', '- a0 + a1 - a2 - a3 - a4 + a6 + a7 - b0 - b1 - b3 - b5 - b6 >= -8', '- a1 + a2 - a3 - a5 + a6 - b0 - b1 - b2 + b3 - b4 - b5 - b6 - b7 >= -9', '- a0 + a1 + a2 - a3 - a4 - a6 - a7 - b1 + b3 - b5 - b6 >= -7', 'a1 + a2 + a3 - a4 - a6 - a7 - b0 - b3 - b6 >= -5', 'a0 - a1 - a2 - a3 - a4 - a6 - a7 - b1 - b2 - b3 + b5 - b6 - b7 >= -10', '- a0 + a1 - a2 - a3 - a5 - a6 + a7 - b2 + b3 - b5 - b6 >= -7', 'a1 - a2 + a3 - a5 - a6 + a7 - b0 - b1 - b2 - b3 - b6 >= -7', '- a0 + a1 - a2 - a3 - a4 - a6 - a7 - b0 - b1 - b2 - b3 - b5 - b6 >= -11', '- a1 - a2 - a3 - a4 - a5 - a6 - a7 - b0 - b2 + b3 - b4 - b5 - b6 - b7 >= -12', '- a2 + a3 - a4 - a5 - a6 - a7 - b0 - b2 + b3 - b5 - b6 + b7 >= -8', 'a4 - a6 + a7 - b2 - b4 >= -2', '- a0 - a1 + b0 + b6 >= -1', 'a4 - a5 + a6 - a7 >= -1', '- a0 - a1 + a5 + a7 - b4 - b6 >= -3', '- a0 - a1 + b0 + b3 >= -1', '- a1 + a5 - a6 + a7 - b5 >= -2', '- a0 - a1 - a3 - a5 - b4 - b5 + b6 - b7 >= -6', '- a1 + a4 - a6 - b5 + b6 >= -2', '- a1 + a4 + a6 + a7 - b5 >= -1', '- a0 - a1 + a3 - a5 - b5 + b6 + b7 >= -3', '- a0 + a4 - a6 - b4 + b5 - b6 >= -3', 'a1 + b0 + b3 - b5 - b6 >= -1', 'a1 + b0 + b5 + b6 >= 1', '- a0 - a1 - a3 - a7 - b0 - b4 + b5 - b6 + b7 >= -6', '- a0 + a3 - a6 - b0 - b4 + b5 - b6 - b7 >= -5', 'a0 + b0 - b5 + b6 >= 0', 'a1 + a3 + a4 + b5 + b6 >= 1', 'a1 - a3 + a4 + b5 - b6 >= -1', '- a0 - a1 - a3 - a4 + a6 + a7 - b0 + b5 - b6 - b7 >= -6', '- a5 - a7 + b4 + b5 - b6 - b7 >= -3', '- a0 + a3 - a4 - a5 + a6 + a7 - b0 - b4 + b5 - b6 + b7 >= -5', 'a0 + a1 + a5 >= 1', '- a0 - a1 - a3 - a6 - b0 - b4 + b5 - b6 + b7 >= -6', '- a1 + a2 + a3 - a5 + a6 - b0 - b2 - b3 - b4 - b5 - b6 - b7 >= -8', '- a1 + a2 - a3 - a5 + a6 - b0 - b2 - b3 - b5 - b6 + b7 >= -7', '- a1 - a2 - a3 - a5 + a6 - b0 - b1 - b3 - b4 - b5 - b6 - b7 >= -10', '- a1 - a2 + a3 - a5 + a6 - b0 - b2 + b3 - b4 - b5 - b6 - b7 >= -8', '- a1 - a2 - a3 - a5 + a6 - b0 - b2 + b3 - b5 - b6 + b7 >= -7', '- a1 + a2 + a3 - a5 + a7 - b0 - b1 - b2 + b3 - b4 - b5 - b6 + b7 >= -7', '- a1 - a2 + a3 - a5 + a7 - b0 - b1 - b2 - b3 - b4 - b5 - b6 + b7 >= -9', '- a0 - a1 + a2 - a3 + a6 - a7 - b1 - b2 + b3 - b4 + b6 - b7 >= -7', '- a0 + a3 - a5 - a7 - b0 + b5 - b6 - b7 >= -5', '- a0 - a1 - a2 - a3 + a6 - a7 - b1 - b3 - b4 + b6 - b7 >= -8', 'a2 + a3 - a5 + a6 - a7 - b0 - b1 - b2 + b3 - b5 - b6 + b7 >= -6', '- a0 + a2 - a3 - a4 - a6 - a7 - b1 - b3 - b4 + b6 - b7 >= -8', '- a2 + a3 - a5 + a6 - a7 - b0 - b1 - b3 - b5 - b6 + b7 >= -7', 'a3 - a7 + b0 + b1 + b5 - b7 >= -1', '- a1 + a2 + a3 + a6 - a7 - b1 - b2 + b3 + b5 - b6 - b7 >= -5', '- a0 - a1 - a2 + a3 + a6 - a7 - b2 + b3 - b4 + b5 - b7 >= -6', '- a3 - a7 + b0 + b1 + b5 + b7 >= -1', '- a1 + a2 - a3 - a4 - a5 - a6 - a7 - b0 - b3 - b4 - b5 - b6 - b7 >= -11', '- a1 + a2 + a3 - a4 - a6 + a7 - b0 - b3 - b5 - b6 - b7 >= -7', '- a1 + a2 - a3 - a4 - a6 + a7 - b0 - b1 + b3 - b5 - b6 - b7 >= -8', '- a1 + a2 + a3 - a4 - a5 - a6 - a7 - b0 - b1 + b3 - b4 - b5 - b6 - b7 >= -10', '- a1 + a2 - a3 - a6 + a7 - b0 - b2 - b3 - b4 - b5 - b6 + b7 >= -8', '- a1 - a2 - a3 - a6 + a7 - b0 - b1 - b2 - b3 - b4 - b5 - b6 - b7 >= -11', 'a2 + a3 - a4 - a5 - a6 - a7 - b0 - b3 - b5 - b6 + b7 >= -7', '- a1 + a2 - a3 - a4 - a5 - a6 - a7 - b0 - b1 + b3 - b5 - b6 + b7 >= -9', '- a1 - a2 + a3 - a4 - a5 - a6 - a7 - b0 - b1 - b2 - b3 - b4 - b5 - b6 - b7 >= -13', '- a1 - a2 - a3 - a4 - a5 - a6 - a7 - b0 - b1 - b2 - b3 - b5 - b6 + b7 >= -12', '- a1 - a2 + a3 - a6 + a7 - b0 - b2 + b3 - b4 - b5 - b6 - b7 >= -8', '- a1 - a2 - a3 - a6 + a7 - b0 - b2 + b3 - b4 - b5 - b6 + b7 >= -8', 'a0 - a1 + a2 - a3 + a6 - b1 - b2 + b3 - b4 + b5 - b6 + b7 >= -5', 'a0 + a2 + a3 - a5 + a7 - b1 - b2 + b3 - b4 + b5 - b6 - b7 >= -5', 'a0 - a2 + a3 - a5 + a7 - b1 - b2 - b3 - b4 + b5 - b7 >= -6', 'a0 - a1 - a2 - a3 - a5 + a7 - b1 - b2 - b3 - b4 + b5 + b7 >= -7', 'a0 - a2 + a3 - a5 + a7 - b2 + b3 - b4 + b5 - b6 + b7 >= -4', 'a1 + a2 + a3 - a4 + a6 + a7 - b0 - b2 - b3 - b5 - b6 >= -5', '- a0 + a1 - a2 - a3 + a6 - a7 - b2 + b3 - b5 - b6 >= -6', '- a0 + a1 + a2 - a3 - a5 - a6 + a7 - b0 - b2 - b3 - b5 - b6 >= -8', 'a1 - a2 + a3 - a4 + a6 + a7 - b2 + b3 - b5 - b6 >= -4', 'a1 + a2 + a3 - a5 - a6 + a7 - b1 - b2 + b3 - b5 - b6 >= -5', '- a0 - a1 + a2 - a3 - a4 + a6 + a7 - b2 - b3 - b4 + b5 + b6 - b7 >= -7', '- a1 + a2 + a3 + a6 - a7 - b2 - b3 - b4 + b5 + b6 - b7 >= -5', 'a4 + b4 + b6 >= 1', '- a1 + a2 - a3 + a6 - a7 - b2 - b3 + b5 + b6 + b7 >= -4', '- a0 + a1 + a2 - a3 - a4 + a6 + a7 - b1 - b2 - b3 - b5 + b6 >= -6', 'a1 + a2 + a3 + a6 - a7 - b1 - b2 - b3 - b5 + b6 >= -4', '- a0 + a1 + a2 - a3 + a6 - a7 - b0 - b1 - b2 + b3 - b5 + b6 >= -6', '- a1 + a2 + a3 - a4 + a6 + a7 - b0 - b1 - b2 + b3 - b4 + b5 + b6 - b7 >= -6', 'a3 - a6 + b0 + b1 + b5 - b7 >= -1', '- a0 + a1 - a2 - a3 + a6 - a7 - b1 - b3 - b5 + b6 >= -6', '- a0 - a1 + a2 + a3 - a4 + a6 + a7 - b2 - b3 + b5 + b6 + b7 >= -4', '- a0 - a1 + a2 - a3 - a4 + a6 + a7 - b1 - b2 + b3 + b5 + b6 + b7 >= -5', '- a0 - a1 - a2 - a3 - a4 + a6 + a7 - b2 + b3 - b4 + b5 - b7 >= -7', '- a3 - a6 + b0 + b1 + b5 + b7 >= -1', '- a1 + a2 + a3 - a5 - a6 + a7 - b2 - b3 - b4 + b5 + b6 - b7 >= -6', '- a0 + a2 - a3 - a5 - a6 + a7 - b1 - b2 + b3 - b4 + b5 + b6 - b7 >= -7', '- a1 - a2 - a3 - a4 + a6 + a7 - b1 - b3 + b5 + b6 + b7 >= -5', '- a0 - a1 + a2 + a3 - a4 - a6 - a7 - b1 + b3 - b4 + b5 - b7 >= -7', '- a1 - a2 - a3 + a6 - a7 - b0 - b2 + b3 + b5 + b6 + b7 >= -5', '- a0 + a1 - a2 - a3 - a4 + a6 + a7 - b0 - b2 + b3 - b5 + b6 >= -6', '- a1 + a2 - a3 - a4 - a5 - a6 + a7 - b3 + b5 + b6 + b7 >= -5', 'a1 - a2 + a3 + a6 - a7 - b0 - b2 + b3 - b5 + b6 >= -4', '- a0 + a2 - a3 - a4 - a6 - a7 - b1 + b3 + b5 + b6 + b7 >= -5', 'a1 + a2 + a3 - a5 - a6 + a7 - b1 - b2 - b3 - b5 + b6 >= -5', '- a0 + a1 + a2 - a3 - a5 - a6 + a7 - b0 - b1 - b2 + b3 + b6 >= -6', '- a1 - a2 + a3 - a4 - a6 - a7 - b1 - b2 - b3 - b4 + b5 + b6 - b7 >= -9', '- a0 - a2 - a3 - a4 - a6 - a7 - b2 + b3 - b4 + b5 + b6 - b7 >= -8', '- a1 - a2 - a3 - a4 - a6 - a7 - b1 - b2 - b3 + b5 + b6 + b7 >= -8', '- a0 - a1 - a2 + a3 - a4 + a6 + a7 - b2 + b3 + b5 + b6 + b7 >= -4', '- a0 - a1 + a2 + a3 - a5 - a6 + a7 - b1 - b2 + b3 + b6 + b7 >= -5', '- a0 - a1 - a2 + a3 - a5 - a6 + a7 - b2 + b3 - b4 + b5 - b7 >= -7', '- a0 - a2 - a3 - a5 - a6 + a7 - b2 + b3 + b5 + b6 + b7 >= -5', 'a0 - a2 + a3 + a6 - b1 - b3 + b5 - b6 - b7 >= -4', 'a0 - a1 - a2 - a3 + a6 - b2 + b3 + b5 - b6 - b7 >= -5', 'a0 - a1 + a2 - a3 - a5 + a7 - b2 - b3 - b4 + b5 - b6 - b7 >= -7', 'a0 + a2 + a3 - a4 - a6 - a7 - b3 + b5 - b6 - b7 >= -5', '- a1 - a2 - a3 + a6 - a7 - b0 - b1 - b3 - b4 + b5 - b6 + b7 >= -8', 'a0 - a1 + a2 - a3 - a4 - a6 - a7 - b3 - b4 + b5 + b7 >= -6', 'a0 - a2 + a3 + a6 - b2 + b3 - b4 + b5 - b6 + b7 >= -3', 'a0 + a2 + a3 - a5 + a7 - b2 - b3 - b4 + b5 - b6 + b7 >= -4', '- a1 + a2 - a3 - a4 - a5 - a6 + a7 - b1 + b3 + b5 - b6 + b7 >= -6', '- a2 + a3 - a4 - a6 - a7 - b0 - b2 + b3 + b5 - b6 - b7 >= -7', '- a1 - a2 - a3 - a4 - a6 - a7 - b2 + b3 - b4 + b5 - b6 + b7 >= -8', 'a1 + a2 + a3 + a6 - a7 - b1 - b2 + b3 - b5 - b6 >= -4', 'a1 - a2 + a3 + a6 - a7 - b0 - b1 - b3 - b5 - b6 >= -6', 'a1 - a2 + a3 - a4 - a6 - a7 - b2 + b3 - b5 - b6 >= -6', '- a0 + a1 + a2 - a3 - a4 - a6 - a7 - b1 - b3 + b6 >= -6', '- a3 + a6 + a7 + b0 + b1 + b5 - b6 - b7 >= -2', 'a2 - a3 + a4 + a6 - a7 - b3 - b5 - b7 >= -4', 'a3 + a6 + a7 + b0 + b1 + b5 - b6 + b7 >= 0', 'a2 + a3 + a4 + a6 - a7 - b3 - b5 + b7 >= -2', '- a2 - a3 + a4 + a6 - a7 + b3 - b5 - b7 >= -4', '- a2 + a3 + a4 + a6 - a7 + b3 - b5 + b7 >= -2', 'a2 + a3 + a4 + a6 - a7 + b3 - b5 - b6 - b7 >= -3', 'a2 - a3 + a4 + a6 - a7 + b3 - b6 + b7 >= -2', '- a2 + a3 + a4 + a6 - a7 - b3 - b5 - b6 - b7 >= -5', '- a2 - a3 + a4 + a6 - a7 - b3 - b6 + b7 >= -4', 'a0 + a3 + b1 + b6 - b7 >= 0', 'a0 - a3 + b1 + b6 + b7 >= 0']
        self.possible_probabilities = {8: ['2', '6', '2_4150', '2_6781', '3', '3_1926', '3_4150', '3_6781', '4', '4_4150', '5', '5_4150', '7'], 4: ["2", "3"]}                                          
        # sort weight based on the magnitude: ['2', '6', '2_4150', '2_6781', '3', '3_1926', '3_4150', '3_6781', '4', '4_4150', '5', '5_4150', '7']
        # sort weight based on frequency of appearance in DDT: ['6', '5', '7', '4', '3', '4_4150', '5_4150', '2', '3_4150', '2_4150', '2_6781', '3_6781', '3_1926']
        self.sbox_inequalties_8bit = {'2' : self.q2_ineqs, '2_4150' : self.q2_4150_ineqs, '2_6781' : self.q2_6781_ineqs, '3' : self.q3_ineqs,\
                                      '3_1926' : self.q3_1926_ineqs, '3_4150' : self.q3_4150_ineqs, '3_6781' : self.q3_6781_ineqs, '4' : self.q4_ineqs,\
                                      '4_4150' : self.q4_4150_ineqs, '5' : self.q5_ineqs, '5_4150' : self.q5_4150_ineqs, '6' : self.q6_ineqs, '7' : self.q7_ineqs}
        # Precomputed constraints for the 4-bit S-box
        # a0, b0: msb
        # a3, b4: lsb
        # 3.0000 p0 + 2.0000 p1
        self.sbox_inequalties_4bit = ['- a2 + a3 - b0 - b1 - b2 - b3 >= -4', '- a1 - a2 - a3 - b1 + b2 - b3 >= -4', '- a1 - a2 - a3 - b1 - b2 + b3 >= -4', '- a2 + a3 - b0 - b1 + b2 + b3 >= -2', '- p0 - p1 >= -1', '- a2 - b0 + p0 >= -1', '- b0 - b1 + p0 >= -1', 'a0 + a1 - a3 + b0 >= 0', '- a0 + a1 + a2 + b1 >= 0', '- a1 + a2 + b0 + b2 >= 0', '- a2 + b0 + b1 + b3 >= 0', 'a0 + a1 + a3 - p0 >= 0', 'a2 - b1 - b2 + p0 >= -1', '- a1 + b0 - b3 + p0 >= -1', '- a2 + b2 + b3 + p0 >= 0', 'a1 + b0 - b1 + p1 >= 0', 'a2 + b1 - b2 + p1 >= 0', 'a1 + b0 - b3 + p1 >= 0', 'a1 - a2 - a3 + b0 - b1 >= -2', 'a0 - a1 + a2 + b1 + b2 >= 0', 'a1 + b0 - b1 - b2 - b3 >= -2', '- a1 + b0 - b1 + b2 - b3 >= -2', '- a1 - a2 + b0 - b2 + b3 >= -2', '- a0 + a3 + b1 - b3 + p0 >= -1', 'a0 + a1 + a2 + a3 - p1 >= 0', 'a0 + a1 + a2 - b0 + p1 >= 0', 'a1 - a2 + a3 - b1 + p1 >= -1', '- a0 + a1 - a3 + b1 + p1 >= -1', 'a0 + a2 - b0 - b2 + p1 >= -1', '- a0 - a1 + a2 + b2 + p1 >= -1', '- a1 + a2 - a3 + b1 - b2 - b3 >= -3', '- a0 + a2 - a3 - b0 + b1 + b3 >= -2', 'a2 + a3 - b0 + b1 - b2 + b3 >= -1', '- a0 - a3 - b0 + b1 - b3 + p1 >= -3', 'a0 + a3 - b0 + b1 - b3 + p1 >= -1', 'a0 - a1 - a3 + b1 + b3 + p1 >= -1', '- a0 - a1 + a3 + b1 + b3 + p1 >= -1']
        # tight lower bound for the number of active S-boxes for 1 <= rounds <= 22
        self.lower_bound = [[1, 2, 5, 8, 12, 16, 26, 36, 41, 46, 51, 55, 58, 61, 66, 75, 82, 88, 92, 96, 102, 108],
                            [0, 0, 1, 2, 3 , 6 , 10, 13, 16, 23, 32, 38, 41, 45, 49, 54, 59, 62, 66, 70, 75 , 79],
                            [0, 0, 0, 0, 1 , 2 , 3 , 6 , 9 , 12, 16, 21, 25, 31, 35, 40, 43, 47, 52, 57, 59 , 64],
                            [0, 0, 0, 0, 0 , 0 , 1 , 2 , 3 , 6 , 10, 13, 16, 19, 24, 27, 31, 35, 43, 45, 48 , 51],
                            [0, 0, 0, 0, 0 , 0 , 0 , 0 , 1 , 2 , 3 , 6 , 9 , 12, 16, 19, 21, 24, 30, 35, 39 , 41]]
    
    def create_objective_function(self, start_round=0, end_round=None):
        '''
        Create the objective function
        '''
        
        if end_round == None:
            end_round = self.rounds
        lp_contents = ""
        minus_log2_p = []
        if self.skipsb == 1:
            start_round = start_round + 1
        if self.cellsize == 8:
            for r in range(start_round, end_round):
                for cell_number in range(16):                
                    for pr in self.possible_probabilities[self.cellsize][0:self.accuracy_threshold]:
                        minus_log2_p.append(pr.replace('_', '.') + ' ' + 'q' + pr + '_' + str(r) + '_' + str(cell_number))
                    if self.exact == True:
                        for pr in self.possible_probabilities[self.cellsize][self.accuracy_threshold:]:
                            minus_log2_p.append(pr.replace('_', '.') + ' ' + 'q' + pr + '_' + str(r) + '_' + str(cell_number))
        elif self.cellsize == 4:
            for r in range(start_round, end_round):
                for cell_number in range(16):
                    for pr in ["2", "3"]:
                        minus_log2_p.append(f"{pr} q{pr}_{r}_{cell_number}")
        else:
            raise Exception("Cell size should be 4 or 8")
        lp_contents += ' + '.join(minus_log2_p)
        return lp_contents

    def create_state_variables(self, r, s):
        '''
        Generate the state variables
        '''

        array = [['' for _ in range(0, self.cellsize)] for _ in range(0, 16)]
        for i in range(0, 16):
            for j in range(0, self.cellsize):                
                array[i][j] = f"{s}_{r}_{i}_{j}"
                self.used_variables.append(array[i][j])
        return array

    def create_half_state_variables(self, r, s):
        '''
        Generate the variables to denote the first half of the state matrix
        '''

        array = [['' for _ in range(0, self.cellsize)] for _ in range(0, 8)]
        for i in range(0, 8):
            for j in range(0, self.cellsize):
                array[i][j] = f"{s}_{r}_{i}_{j}"
                self.used_variables.append(array[i][j])
        return array

    def flatten(self, state_array):
        '''
        Get a state array and output a flatten list
        '''

        flat_list = []
        for cell_number in range(len(state_array)):
            for bit_number in range(len(state_array[0])):
                flat_list.append(state_array[cell_number][bit_number])
        return flat_list

    def xor(self, a, b, c):
        '''
        Generate the constraints of a binary XOR
        a xor b = c can be modeled with 4 inequalities (without definition of dummy variable) by removing all impossible vectors (a, b, c)
        '''

        lp_contents = ""
        lp_contents += f"{a} + {b} - {c} >= 0\n"                
        lp_contents += f"{a} - {b} + {c} >= 0\n"
        lp_contents += f"-1 {a} + {b} + {c} >= 0\n"        
        lp_contents += f"-1 {a} - {b} - {c} >= -2\n"        
        return lp_contents

    def xor3(self, b, a2, a1, a0):
        '''
        Generate the constraints of a three-input XOR  (b = a0 xor a1 xor a2)    
        b - a2 - a1 - a0 >= -2
        - b + a2 - a1 - a0 >= -2
        - b - a2 + a1 - a0 >= -2
        b + a2 + a1 - a0 >= 0
        - b - a2 - a1 + a0 >= -2
        b + a2 - a1 + a0 >= 0
        b - a2 + a1 + a0 >= 0
        - b + a2 + a1 + a0 >= 0
        The above inequalities are derived with QuineMcCluskey algorithm
        '''

        lp_contents = ""
        lp_contents += f"{b} - {a2} - {a1} - {a0} >= -2\n"
        lp_contents += f"-1 {b} + {a2} - {a1} - {a0} >= -2\n"
        lp_contents += f"-1 {b} - {a2} + {a1} - {a0} >= -2\n"
        lp_contents += f"{b} + {a2} + {a1} - {a0} >= 0\n"
        lp_contents += f"-1 {b} - {a2} - {a1} + {a0} >= -2\n"
        lp_contents += f"{b} + {a2} - {a1} + {a0} >= 0\n"
        lp_contents += f"{b} - {a2} + {a1} + {a0} >= 0\n"
        lp_contents += f"-1 {b} + {a2} + {a1} + {a0} >= 0\n"
        return lp_contents
    
    def equality(self, x, y):
        '''
        Generate the MILP constraints modeling the equality of two bits
        '''

        lp_contents = f"{x} - {y} = 0\n"
        return lp_contents

    def mix_columns(self, x, y):
        '''
        Model the MixColumns
        '''

        lp_contents = ""
        for cell_number in range(4):            
            for bit_number in range(self.cellsize):
                lp_contents += self.xor3(y[cell_number][bit_number], x[cell_number][bit_number], x[cell_number + 8][bit_number], x[cell_number + 12][bit_number])
                lp_contents += self.equality(y[cell_number + 4][bit_number], x[cell_number][bit_number])
                lp_contents += self.xor(x[cell_number + 4][bit_number], x[cell_number + 8][bit_number], y[cell_number + 8][bit_number])
                lp_contents += self.xor(x[cell_number][bit_number], x[cell_number + 8][bit_number], y[cell_number + 12][bit_number])
        return lp_contents

    def add_tweakey(self, tk, x, y):
        '''
        Model the AddTweaey
        tk xor x = y
        '''

        lp_contents = ""
        for cell_number in range(8):
            for bit_number in range(self.cellsize):
                lp_contents += self.xor(tk[cell_number][bit_number], x[cell_number][bit_number], y[cell_number][bit_number])
        for cell_number in range(8, 16):
            for bit_number in range(self.cellsize):
                lp_contents += self.equality(x[cell_number][bit_number], y[cell_number][bit_number])
        return lp_contents

    def shift_rows(self, state):
        '''
        ShiftRows
        This function does not generate MILP constraints
        '''

        temp = [0]*16
        for i in range(16):
            temp[i] = state[self.permuteation[i]]
        return temp
    
    def permute_tweakey(self, state):
        '''
        Tweakey permutation (Q)
        This function does not generate MILP constraints
        '''

        temp = [0]*16
        for i in range(16):
            temp[i] = state[self.tk_permutation[i]]
        return temp
    
    def subcells_8bit(self, x, y, r):
        '''
        Model the 8-bit S-box
        x: input state
        y: output state
        y = S(x)
        r: indicates the round number
        '''

        lp_contents = ""
        for cell_number in range(16):
            q_byte_variables = []           
            for pr in self.possible_probabilities[self.cellsize][0:7]:
                q = 'q' + pr + '_' + str(r) + '_' + str(cell_number)
                q_byte_variables.append(q)
                for ineq in self.sbox_inequalties_8bit[pr]:
                    for i in range(8):
                        ineq = ineq.replace('a' + str(i), x[cell_number][i])
                        ineq = ineq.replace('b' + str(i), y[cell_number][i])
                    # Big-M Method:
                    ineq = ineq.split(' >= ')
                    ineq = ineq[0] + ' - ' + str(self.big_m) + ' ' + q + ' >= ' + str(int(ineq[1]) - self.big_m) + '\n'
                    # Indicator constraint:
                    # ineq  = q + ' = 1'+ ' -> ' + ineq + '\n'                    
                    lp_contents += ineq
            if self.exact == True:                            
                for pr in self.possible_probabilities[self.cellsize][self.accuracy_threshold:]:
                    q = 'q' + pr + '_' + str(r) + '_' + str(cell_number)
                    q_byte_variables.append(q)
                    for ineq in self.sbox_inequalties_8bit[pr]:
                        for i in range(8):
                            ineq = ineq.replace('a' + str(i), x[cell_number][i])
                            ineq = ineq.replace('b' + str(i), y[cell_number][i])
                        # Big-M method:
                        ineq = ineq.split(' >= ')
                        ineq = ineq[0] + ' - ' + str(self.big_m) + ' ' + q + ' >= ' + str(int(ineq[1]) - self.big_m) + '\n'
                        # Indicator constraint:
                        # ineq  = q + ' = 1'+ ' -> ' + ineq + '\n'
                        lp_contents += ineq                        
            # q = sum(qi)
            self.used_variables.extend(q_byte_variables)
            q_byte = 'q_' + str(r) + '_' + str(cell_number)
            self.used_variables.append(q_byte)
            ineq = q_byte + ' - ' + ' - '.join(q_byte_variables) + ' = 0\n'            
            lp_contents += ineq
            # q[cell_number] >= x[cell_number][i] for i in range(8)
            for i in range(8):
                ineq = q_byte + ' - ' + x[cell_number][i] + ' >= 0\n'                
                lp_contents += ineq
            # sum(x[cell_number][i] for i in range(8)) >= q[cell_number]
            ineq = ' + '.join(x[cell_number]) + ' - ' + q_byte + ' >= 0\n'            
            lp_contents += ineq
            # Model the bijectivity of Sbox (x = 0 => y = 0, y = 0 => x = 0)
            ineq = '8 ' + ' + 8 '.join(y[cell_number]) + ' - ' + ' - '.join(x[cell_number]) + ' >= 0\n'            
            lp_contents += ineq
            ineq = '8 ' + ' + 8 '.join(x[cell_number]) + ' - ' + ' - '.join(y[cell_number]) + ' >= 0\n'            
            lp_contents += ineq
        return lp_contents
    
    def subcells_4bit(self, x, y, r):
        '''
        Model the 4-bit S-box
        '''

        lp_contents = ""
        for cell_number in range(16):
            for ineq in self.sbox_inequalties_4bit:
                for i in range(4):
                    ineq = ineq.replace(f"a{i}", x[cell_number][i])
                    ineq = ineq.replace(f"b{i}", y[cell_number][i])
                ineq = ineq.replace("p0", f"q3_{r}_{cell_number}")
                ineq = ineq.replace("p1", f"q2_{r}_{cell_number}")
                lp_contents += ineq + "\n"
            self.used_variables.append(f"q3_{r}_{cell_number}")
            self.used_variables.append(f"q2_{r}_{cell_number}")            
            lp_contents += '4 ' + ' + 4 '.join(y[cell_number]) + ' - ' + ' - '.join(x[cell_number]) + ' >= 0\n'
            lp_contents += '4 ' + ' + 4 '.join(x[cell_number]) + ' - ' + ' - '.join(y[cell_number]) + ' >= 0\n'            
        return lp_contents
    
    def subcells(self, x, y, r, cellsize):
        '''
        decide between 4-bit or 8-bit S-box
        '''

        if cellsize == 4:
            return self.subcells_4bit(x, y, r)
        elif cellsize == 8:
            return self.subcells_8bit(x, y, r)
    
    def lfsr_tk2(self, a, b):
        '''
        Generate the MILP constraints corresponding to the 8-bit LFSR used through the tweakey path
        It is supposed that a, and b, are 8-bit vectors
        Seminal paper: (x7||x6||x5||x4||x3||x2||x1||x0) -> (x6||x5||x4||x3||x2||x1||x0||x7 xor x5) where x0 is the LSB
        a[0]: msb, a[7]: lsb
        (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]) -> (a[1], a[2], a[3], a[4], a[5], a[6], a[0] xor a[2]) = (b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7])
        '''

        lp_contents = ""
        for i in range(self.cellsize - 1):
            lp_contents += self.equality(a[i + 1], b[i])
        if self.cellsize == 8:
            lp_contents += self.xor(a[0], a[2], b[7])
        elif self.cellsize == 4:
            lp_contents += self.xor(a[0], a[1], b[3])
        return lp_contents

    def lfsr_tk3(self, a, b):
        '''
        Generate the MILP constraints corresponding to the 8-bit LFSR used through the tweakey path
        It is supposed that a, and b, are 8-bit vectors
        Seminal paper: (x7||x6||x5||x4||x3||x2||x1||x0) -> (x0 xor x6||x7||x6||x5||x4||x3||x2||x1) where x0 is the LSB
        a[0]: msb, a[7]: lsb
        (a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]) -> (a[7] xor a[1], a[0], a[1], a[2], a[3], a[4], a[5], a[6]) = (b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7])
        '''
        
        lp_contents = ""        
        for i in range(0, self.cellsize - 1):
            lp_contents += self.equality(a[i], b[i + 1])
        if self.cellsize == 8:
            lp_contents += self.xor(a[7], a[1], b[0])
        elif self.cellsize == 4:
            lp_contents += self.xor(a[3], a[0], b[0])
        return lp_contents
    
    def lfsr_tk4(self, a, b):
        '''
        Generate the MILP constraints for the linear map employed in the fourth tweakey line of SKINNYe-v2
        Note that this linear map is not an LFSR!
        Reference: https://ia.cr/2020/542
        (a[0], a[1], a[2], a[3]) -> (a[2], a[3], a[0] xor a[1], a[1] xor a[2]) = (b[0], b[1], b[2], b[3])
        '''

        lp_contents = ""
        lp_contents += self.xor(a[0], a[1], b[2])
        lp_contents += self.xor(a[1], a[2], b[3])
        lp_contents += self.equality(a[2], b[0])
        lp_contents += self.equality(a[3], b[1])
        return lp_contents

    def tweakey_schedule(self):
        '''
        Model the difference propagation through the tweakey schedule
        '''

        lp_contents = ""
        if self.variant >= 1:
            tk1 = self.create_state_variables(0, 'tk1')
        if self.variant >= 2:
            tk2 = self.create_state_variables(0, 'tk2')
        if self.variant >= 3:
            tk3 = self.create_state_variables(0, 'tk3')
        if self.variant >= 4:
            tk4 = self.create_state_variables(0, 'tk4')
            tk34 = self.create_half_state_variables(0, 'tk34')
        tk = self.create_half_state_variables(0, 'tk')
        # model the round tweakey generation in the first round: TK = FirstHalf(TK1) xor FirstHalf(TK2) xor FirstHalf(TK3)
        for cell_number in range(8):
            for bit_number in range(self.cellsize):
                if self.variant == 0:
                    lp_contents += self.equality(tk[cell_number][bit_number], 0) # single-tweakey differential analysis
                elif self.variant == 1:
                    lp_contents += self.equality(tk[cell_number][bit_number], tk1[cell_number][bit_number])
                elif self.variant == 2:
                    lp_contents += self.xor(tk[cell_number][bit_number], tk1[cell_number][bit_number], tk2[cell_number][bit_number])
                elif self.variant == 3:
                    lp_contents += self.xor3(tk[cell_number][bit_number], tk1[cell_number][bit_number], tk2[cell_number][bit_number], tk3[cell_number][bit_number])
                elif self.variant == 4:
                    lp_contents += self.xor(tk34[cell_number][bit_number], tk3[cell_number][bit_number], tk4[cell_number][bit_number])
                    lp_contents += self.xor3(tk[cell_number][bit_number], tk1[cell_number][bit_number], tk2[cell_number][bit_number], tk34[cell_number][bit_number])
        for r in range(1, self.rounds):
            if self.variant >= 1:
                ptk1 = self.permute_tweakey(tk1)
                tk1 = self.create_state_variables(r, 'tk1')            
                # tk1 = ptk1
                for cell_number in range(16):
                    for bit_number in range(self.cellsize):
                        lp_contents += self.equality(tk1[cell_number][bit_number], ptk1[cell_number][bit_number])
            if self.variant >= 2:
                # Apply LFSR to the the first half of ptk2
                ptk2 = self.permute_tweakey(tk2)
                tk2 = self.create_state_variables(r, 'tk2')
                for cell_number in range(8):
                    lp_contents += self.lfsr_tk2(ptk2[cell_number], tk2[cell_number])
                for cell_number in range(8, 16):
                    for bit_number in range(self.cellsize):
                        lp_contents += self.equality(ptk2[cell_number][bit_number], tk2[cell_number][bit_number])
            if self.variant >= 3:
                # Apply LFSR to the first half of ptk3
                ptk3 = self.permute_tweakey(tk3)
                tk3 = self.create_state_variables(r, 'tk3')
                for cell_number in range(8):
                    lp_contents += self.lfsr_tk3(ptk3[cell_number], tk3[cell_number])
                for cell_number in range(8, 16):
                    for bit_number in range(self.cellsize):
                        lp_contents += self.equality(ptk3[cell_number][bit_number], tk3[cell_number][bit_number])
            if self.variant >= 4:
                # Apply the linear map to the first half of ptk4
                ptk4 = self.permute_tweakey(tk4)
                tk4 = self.create_state_variables(r, 'tk4')
                tk34 = self.create_half_state_variables(r, 'tk34')
                for cell_number in range(8):
                    lp_contents += self.lfsr_tk4(ptk4[cell_number], tk4[cell_number])
                for cell_number in range(8, 16):
                    for bit_number in range(self.cellsize):
                        lp_contents += self.equality(ptk4[cell_number][bit_number], tk4[cell_number][bit_number])
            tk = self.create_half_state_variables(r, 'tk')
            # model the round tweakey generation: TK = FirstHalf(TK1) xor FirstHalf(TK2) xor FirstHalf(TK3) xor FirstHalf(TK4)
            for cell_number in range(8):
                for bit_number in range(self.cellsize):
                    if self.variant == 0:
                        lp_contents += self.equality(tk[cell_number][bit_number], 0)
                    elif self.variant == 1:
                        lp_contents += self.equality(tk[cell_number][bit_number], tk1[cell_number][bit_number])
                    elif self.variant == 2:
                        lp_contents += self.xor(tk[cell_number][bit_number], tk1[cell_number][bit_number], tk2[cell_number][bit_number])
                    elif self.variant == 3:
                        lp_contents += self.xor3(tk[cell_number][bit_number], tk1[cell_number][bit_number], tk2[cell_number][bit_number], tk3[cell_number][bit_number])
                    elif self.variant == 4:
                        lp_contents += self.xor(tk34[cell_number][bit_number], tk3[cell_number][bit_number], tk4[cell_number][bit_number])
                        lp_contents += self.xor3(tk[cell_number][bit_number], tk1[cell_number][bit_number], tk2[cell_number][bit_number], tk34[cell_number][bit_number])
        return lp_contents
    
    def encryption(self):
        '''
        Generate the MILP constraints modeling the propagation of differences through the encryption data path
        '''
        
        lp_contents = ""
        for r in range(self.rounds):
            x = self.create_state_variables(r, 'x')
            y = self.create_state_variables(r, 'y')
            z = self.create_state_variables(r, 'z')
            tk = self.create_half_state_variables(r, 'tk')
            x_next = self.create_state_variables(r + 1, 'x')
            lp_contents += self.subcells(x, y, r, self.cellsize)
            lp_contents += self.add_tweakey(tk, y, z)
            pz = self.shift_rows(z)
            lp_contents += self.mix_columns(pz, x_next)
        return lp_contents
    
    def exclude_trivial_trail(self):
        lp_contents = ""
        x = self.create_state_variables(0, 'x')
        if self.variant == 0:
            lp_contents += " + ".join(self.flatten(x)) + " >= 1\n"
        elif self.variant == 1:
            tk1 = self.create_state_variables(0, 'tk1')
            temp = self.flatten(x) + self.flatten(tk1)
            lp_contents += " + ".join(temp) + " >= 1\n"
        elif self.variant == 2:
            tk1 = self.create_state_variables(0, 'tk1')
            tk2 = self.create_state_variables(0, 'tk2')
            temp = self.flatten(x) + self.flatten(tk1) + self.flatten(tk2)
            lp_contents += " + ".join(temp) + " >= 1\n"
        elif self.variant == 3:
            tk1 = self.create_state_variables(0, 'tk1')
            tk2 = self.create_state_variables(0, 'tk2')
            tk3 = self.create_state_variables(0, 'tk3')
            temp = self.flatten(x) + self.flatten(tk1) + self.flatten(tk2) + self.flatten(tk3)
            lp_contents += " + ".join(temp) + " >= 1\n"
        elif self.variant == 4:
            tk1 = self.create_state_variables(0, 'tk1')
            tk2 = self.create_state_variables(0, 'tk2')
            tk3 = self.create_state_variables(0, 'tk3')
            tk4 = self.create_state_variables(0, 'tk4')
            temp = self.flatten(x) + self.flatten(tk1) + self.flatten(tk2) + self.flatten(tk3) + self.flatten(tk4)
            lp_contents += " + ".join(temp) + " >= 1\n"
        return lp_contents
    
    def declare_fixed_variables(self):
        lp_contents = ""
        for cond in self.fixed_variables.items():            
            var = cond[0]
            val = cond[1]
            var = var.split('_')
            if len(var) == 2:                
                state_vars = self.create_state_variables(var[1], var[0])
                state_vars = self.flatten(state_vars)
                if "X" not in val:
                    state_values = list(bin(int(val, 16))[2:].zfill(self.cellsize*16))
                    for i in range(self.cellsize*16):
                        lp_contents += f"{state_vars[i]} = {state_values[i]}\n"
                else:
                    fixed_positions = [i for i in range(len(val)) if val[i] != "X"]
                    for i in fixed_positions:
                        cell_value = list(bin(int(val[i], 16))[2:].zfill(self.cellsize))
                        for j in range(self.cellsize):
                            lp_contents += f"{state_vars[i*self.cellsize + j]} = {cell_value[j]}\n"
                    
            elif len(var) == 3:
                state_vars = [f"{var[0]}_{var[1]}_{var[2]}_{i}" for i in range(self.cellsize)]
                if val != "Y":
                    state_values = list(bin(int(val, 16))[2:].zfill(self.cellsize))
                    for i in range(self.cellsize):
                        lp_contents += f"{state_vars[i]} = {state_values[i]}\n"
                elif val == "Y":
                    lp_contents += " + ".join(state_vars) + " >= 1\n"
            elif len(var) == 4:
                lp_contents += f"{cond[0]} = {cond[1]}\n"
        return lp_contents

    def declare_variables_type(self):
        '''
        Specifying variables' type in the LP file
        '''
        
        lp_contents = 'binary\n'
        self.used_variables = list(set(self.used_variables))
        for var in self.used_variables:
            lp_contents += var + '\n'            
        lp_contents += "end\n"
        return lp_contents

    def make_model(self):
        '''
        Generate the MILP model of Skinny-128-256 for differential cryptanalysis
        '''
        
        lp_contents = ""
        print('Generating the MILP model ...')
        lp_contents += "minimize\n"
        self.obj_func = self.create_objective_function()
        lp_contents += self.obj_func
        lp_contents += "\nsubject to\n"
        if self.rounds <= 22:
            if self.skipsb:
                if self.rounds >= 2:
                    lp_contents += f"{self.obj_func} >= {2*self.lower_bound[self.variant][self.rounds - 2]}\n"                
            else:
                lp_contents += f"{self.obj_func} >= {2*self.lower_bound[self.variant][self.rounds - 1]}\n"
        if self.upperbound1 != None:
            lp_contents += f"{self.obj_func} <= {self.upperbound1}\n"
        if self.start_round != None and self.end_round != None and self.upperbound2 != None:
            halfway_weight = self.create_objective_function(start_round=self.start_round, end_round = self.end_round)
            lp_contents += f"{halfway_weight} <= {self.upperbound2}\n"
        lp_contents += self.exclude_trivial_trail()
        lp_contents += self.tweakey_schedule()
        lp_contents += self.encryption()
        lp_contents += self.declare_fixed_variables()
        lp_contents += self.declare_variables_type() 
        if os.path.exists(self.model_filename):
            os.remove(self.model_filename)
        with open(self.model_filename, 'w') as fileobj:
            fileobj.write(lp_contents)
        print(f"MILP model was written into {self.model_filename}\n")
    
    def generate_tweakey(self, total_rounds, fixed_round_tweakey, target_round):
        '''
        Generate the tweakey schedule for the given number of rounds
        '''
        
        lp_contents = "subject to\n"
        self.used_variables = []
        self.fixed_variables = fixed_round_tweakey
        self.rounds = total_rounds
        lp_contents += self.tweakey_schedule()
        lp_contents += self.declare_fixed_variables()
        lp_contents += self.declare_variables_type()
        with open("temp.lp", "w") as fileobj:
            fileobj.write(lp_contents)
        milp_model = read("temp.lp")
        os.remove("temp.lp")
        milp_model.Params.OutputFlag = False
        milp_model.optimize()
        characteristic = ""
        if self.variant >= 1:
            tk1 = self.create_state_variables(target_round, 'tk1')
            tk1_flat = self.flatten(tk1)
            tk1_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(milp_model.getVarByName(t).X)), tk1_flat))), 2))[2:].zfill(self.cellsize*4)
            characteristic += tk1_value
        if self.variant >= 2:
            tk2 = self.flatten(self.create_state_variables(target_round, 'tk2'))
            tk2_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(milp_model.getVarByName(t).X)), tk2))), 2))[2:].zfill(self.cellsize*4)
            characteristic += tk2_value
        if self.variant >= 3:
            tk3 = self.flatten(self.create_state_variables(target_round, 'tk3'))
            tk3_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(milp_model.getVarByName(t).X)), tk3))), 2))[2:].zfill(self.cellsize*4)
            characteristic += tk3_value 
        if self.variant >= 4:
            tk4 = self.flatten(self.create_state_variables(target_round, 'tk4'))
            tk4_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(milp_model.getVarByName(t).X)), tk4))), 2))[2:].zfill(self.cellsize*4)
            characteristic += tk4_value 
        return characteristic  

    
    def find_characteristic(self):
        '''
        Find the best differential trail under the given constraints, e.g., satisfying an activeness pattern
        '''
        
        status = False
        if self.time_limit != -1:
            self.model.Params.TIME_LIMIT = self.time_limit
        obj = self.model.getObjective()
        # Consider the start_weight
        if self.start_weight != None:
            self.model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        time_start = time.time()
        #self.model.Params.Threads = 16
        #self.model.Params.PreSolve = 0
        self.model.Params.OutputFlag = True
        self.model.optimize()
        if (self.model.Status in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.INTERRUPTED, GRB.SOLUTION_LIMIT]):
            # obj = self.model.getObjective()
            # objVal = obj.getValue()
            self.total_weight = self.model.objVal
            print("\nThe probability of the best differential characteristic: 2^-(%s)" % self.total_weight)
            print("\nDifferential trail:\n")
            diff_trail = self.parse_solver_output()
            self.print_trail(diff_trail)
            status = True
        elif self.model.Status == GRB.INFEASIBLE:
            print("The model is infeasible!")
        else:
            print("Unknown error!")
        time_end = time.time()
        print("Time used = {:0.02f}".format(time_end - time_start))
        return status

    def find_multiple_characteristics(self):
        '''
        Find multiple differential trails for the given number of rounds (and the given fixed input/output differences)
        '''
        
        status = False
        if self.time_limit != -1:
            self.model.Params.TIME_LIMIT = self.time_limit
        #m.setParam(GRB.Param.Threads, 16)
        obj = self.model.getObjective()
        # Consider the start_weight
        if self.start_weight != None:
            self.model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        self.model.Params.OutputFlag = False 
        # You can use the PoolSearchMode parameter to control the approach used to find solutions. In its default setting (0), the MIP search simply aims to find one optimal solution. Setting the parameter
        # to 1 causes the MIP search to expend additional effort to find more solutions, but in a non-systematic way. 
        # You will get more solutions, but not necessarily the best solutions. Setting the parameter to 2 causes the MIP to do a systematic search for the n best solutions. For both non-default settings, 
        # the PoolSolutions parameter sets the target for the number of solutions to find.
        self.model.Params.PoolSearchMode = 2
        self.model.Params.PoolSolutions = 10
        time_start = time.time()
        self.model.optimize()
        if (self.model.Status == GRB.OPTIMAL or self.model.Status == GRB.TIME_LIMIT or self.model.Status == GRB.INTERRUPTED):
            status = True
            # First Method:
            number_of_trails = 10
            for sol_number in range(number_of_trails):
                if (self.model.Status == GRB.OPTIMAL):
                    self.total_weight = self.model.PoolObjVal
                    diff_trail = self.parse_solver_output()
                    self.print_trail(diff_trail)
                elif (self.model.Status == GRB.TIME_LIMIT or self.model.Status == GRB.INTERRUPTED):
                    self.total_weight = self.model.PoolObjVal
                    diff_trail = self.parse_solver_output()
                    self.print_trail(diff_trail)
                    break
                else:
                    break
                self.exclude_the_previous_sol()
                self.model.optimize()
            # Second Method:
            # number_of_trails = self.model.SolCount
            # for sol_number in range(number_of_trails):
            #     self.model.Params.SolutionNumber = sol_number
            #     # PoolObjVal : This attribute is used to query the objective value of the <span>$</span>k<span>$</span>-th solution stored in the pool of feasible solutions found so far for the problem
            #     self.total_weight = self.model.PoolObjVal                
            #     diff_trail = self.parse_solver_output()
            #     self.print_trail(diff_trail)
        elif self.model.Status == GRB.INFEASIBLE:
            print("The model is infeasible!")
        else:
            print("Unknown error!")
        time_end = time.time()
        print("Time used = {:0.02f}".format(time_end - time_start))
        return status
    
    def compute_differential_effect(self, log=1):
        '''
        Compute the differential effect for a given input/output differences

        Some general information about Gurobi:

        PoolSolutions: It controls the size of the solution pool. Changing this parameter won't affect the number of solutions that are found - 
        it simply determines how many of those are retained

        You can use the PoolSearchMode parameter to control the approach used to find solutions. In its default setting (0), the MIP search simply aims to find one optimal solution. 
        Setting the parameter to 2 causes the MIP to do a systematic search for the n best solutions. With a setting of 2, it will find the n best solutions, 
        where n is determined by the value of the PoolSolutions parameter        

        SolCount: Number of solutions found during the most recent optimization.
        
        Model status:
        LOADED	1	Model is loaded, but no solution information is available.
        OPTIMAL	2	Model was solved to optimality (subject to tolerances), and an optimal solution is available.
        INFEASIBLE	3	Model was proven to be infeasible.
        '''
        
        status = False
        if self.time_limit != -1:
            self.model.Params.TIME_LIMIT = self.time_limit
        #self.model.Params.PreSolve = 0 # Activating this flag causes the performance to be decreased        
        self.model.Params.PoolSearchMode = 2
        self.model.Params.PoolSolutions = 1
        self.model.Params.OutputFlag = False                
        obj = self.model.getObjective()
        # Consider the start_weight
        if self.start_weight != None:            
            self.model.addConstr(obj >= self.start_weight, 'start_weight_constraint')       
        time_start = time.time()
        self.model.optimize()
        if (self.model.Status == GRB.OPTIMAL):
            status = True
            self.total_weight = self.model.objVal
            diff_prob = 0
            print('\n')
            while (self.model.Status == GRB.OPTIMAL and self.total_weight <= self.end_weight):
                self.total_weight = self.model.objVal
                self.model.Params.PoolSolutions = 2000000000 #GRB.MAXIN, Default value for PoolSolutions: 10                
                temp_constraint = self.model.addConstr(obj == self.total_weight, name='temp_constraint')
                self.model.update()
                #self.model.Params.PreSolve = 1
                self.model.optimize()
                diff_prob += math.pow(2, -self.total_weight) * self.model.SolCount
                time_end = time.time()
                if log == 1:
                    print('Current weight: %s' % str(self.total_weight))
                    print('Number of trails: %s' % str(self.model.SolCount))
                    print('\tCurrent Probability: 2^(' + str(math.log(diff_prob, 2)) + ')')
                    print('Time used = %0.4f seconds\n' % (time_end - time_start))
                self.model.remove(temp_constraint)
                self.model.Params.PoolSolutions = 1                
                self.model.addConstr(obj >= (self.total_weight + self.eps))
                #self.model.Params.PreSolve = 0
                self.model.optimize()
        elif (self.model.Status == GRB.INFEASIBLE):
            print('The model is infeasible!')
            return status
        else: 
            print('Unknown Error!')
            return status
        print("Total weight = {:0.02f}".format(math.log(diff_prob, 2)))
        return math.log(diff_prob, 2)
            
    def compute_differential_effect_classic_method(self):
        status = False
        if self.time_limit != -1:
            self.model.Params.TIME_LIMIT = self.time_limit
        self.model.Params.OutputFlag = False
        # self.model.printStats()
        # Consider the start_weight
        obj = self.model.getObjective()
        if self.start_weight != None:            
            self.model.addConstr(obj >= self.start_weight, 'start_weight_constraint')       
        time_start = time.time()
        self.model.optimize()
        self.model.Params.Quad = 1
        sol_dict = dict()
        if (self.model.Status == GRB.OPTIMAL):
            status = True
            self.total_weight = self.model.objVal
            diff_prob = 0
            print('\n')
            while (self.model.Status == GRB.OPTIMAL and self.total_weight <= self.end_weight):  
                self.total_weight = self.model.objVal
                diff_prob += math.pow(2, -self.total_weight)
                total_weight_st = 'ntrails_%0.2f' % self.total_weight
                sol_dict[total_weight_st] = sol_dict.get(total_weight_st, 0) + 1
                print('Current weight: %s' % str(self.total_weight))
                print('Number of trails: %d' % sol_dict[total_weight_st])
                print('\tCurrent Probability: 2^(' + str(math.log(diff_prob, 2)) + ')')
                time_end = time.time()
                print('Time used = %0.4f seconds\n' % (time_end - time_start))           
                self.exclude_the_previous_sol()
                self.model.optimize()
        elif (self.model.Status == GRB.INFEASIBLE):
            print('The model is infeasible!')
        else: 
            print('Unknown Error!')
        return status

    def exclude_the_previous_sol(self):
        '''
        Let x{S} be the binary variables. Suppose you have a binary solution x* in available from the most recent optimization. 
        Let N be the subset of S such that x*[n] = 1 for all n in N
        Then, add the following constraint:
        sum{n in N} x[n] - sum{s in S-N} x[s] <= |N|-1
        '''

        all_vars = self.model.getVars()
        nonzero_vars = [v for v in all_vars if v.x == 1]
        zero_vars = [v for v in all_vars if v.x == 0]
        support = len(nonzero_vars)
        first_term = sum(nonzero_vars)
        second_term = sum(zero_vars)
        lhs = first_term - second_term
        self.model.addConstr(lhs <= support - 1)

    def solve(self, log=1, solution_limit=None, mip_focus=None):        
        self.model = read(self.model_filename)
        if solution_limit != None:
            self.model.Params.SolutionLimit = solution_limit
        if mip_focus != None:
            self.model.Params.MIPFocus = mip_focus
        status = False
        if self.mode == 0:
            status = self.find_characteristic()
        elif self.mode == 1:                     
            status = self.find_multiple_characteristics()
        elif self.mode == 2:
            status = self.compute_differential_effect(log)
            #self.compute_differential_effect_classic_method()
        else:
            print("mode should be in [0, 1, 2]")
        os.remove(self.model_filename)
        return status

    def parse_solver_output(self):
        '''
        Extract the differential characteristic from the solver output
        '''

        characteristic = dict()
        for r in range(self.rounds + 1):
            x = self.flatten(self.create_state_variables(r, 'x'))            
            x_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(self.model.getVarByName(t).X)), x))), 2))[2:].zfill(self.cellsize*4)
            characteristic['x_' + str(r)] = x_value
        tk1 = self.create_state_variables(0, 'tk1')
        for r in range(self.rounds):
            y = self.flatten(self.create_state_variables(r, 'y'))
            y_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(self.model.getVarByName(t).X)), y))), 2))[2:].zfill(self.cellsize*4)
            characteristic['y_' + str(r)] = y_value
            z = self.flatten(self.create_state_variables(r, 'z'))
            z_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(self.model.getVarByName(t).X)), z))), 2))[2:].zfill(self.cellsize*4)        
            characteristic['z_' + str(r)] = z_value
            if self.variant >= 1:
                tk1_flat = self.flatten(tk1)
                tk1_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(self.model.getVarByName(t).X)), tk1_flat))), 2))[2:].zfill(self.cellsize*4)
                characteristic['tk1_' + str(r)] = tk1_value
                tk1 = self.permute_tweakey(tk1)
            if self.variant >= 2:
                tk2 = self.flatten(self.create_state_variables(r, 'tk2'))
                tk2_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(self.model.getVarByName(t).X)), tk2))), 2))[2:].zfill(self.cellsize*4)
                characteristic['tk2_' + str(r)] = tk2_value
            if self.variant >= 3:
                tk3 = self.flatten(self.create_state_variables(r, 'tk3'))
                tk3_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(self.model.getVarByName(t).X)), tk3))), 2))[2:].zfill(self.cellsize*4)
                characteristic['tk3_' + str(r)] = tk3_value 
            if self.variant >= 4:
                tk4 = self.flatten(self.create_state_variables(r, 'tk4'))
                tk4_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(self.model.getVarByName(t).X)), tk4))), 2))[2:].zfill(self.cellsize*4)
                characteristic['tk4_' + str(r)] = tk4_value                      
            tk = self.flatten(self.create_half_state_variables(r, 'tk'))[0:64]
            tk_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(self.model.getVarByName(t).X)), tk))), 2))[2:].zfill(self.cellsize*2)
            characteristic['tk_' + str(r)] = tk_value        
            round_probability = 0
            for cell_number in range(16):
                if self.cellsize == 8 and self.exact == False:
                    for pr in self.possible_probabilities[self.cellsize][0:self.accuracy_threshold]:
                        q = 'q' + pr + '_' + str(r) + '_' + str(cell_number)
                        round_probability += float(pr.replace('_', '.')) * int(self.model.getVarByName(q).X)
                else:
                    for pr in self.possible_probabilities[self.cellsize]:
                        q = 'q' + pr + '_' + str(r) + '_' + str(cell_number)
                        round_probability += float(pr.replace('_', '.')) * int(self.model.getVarByName(q).X)
            characteristic['pr_' +str(r)] = '-' + str(round_probability)
        
        return characteristic
    
    def print_trail(self, diff_trail):
        '''
        Print out the obtained differential trail
        '''
        
        if self.variant == 0:
            header = ['x', 'y', 'z', 'tk', 'pr']
        elif self.variant == 1:
            header = ['x', 'y', 'z', 'tk1', 'tk', 'pr']
        elif self.variant == 2:
            header = ['x', 'y', 'z', 'tk1', 'tk2', 'tk', 'pr']
        elif self.variant == 3:
            header = ['x', 'y', 'z', 'tk1', 'tk2', 'tk3', 'tk', 'pr']
        elif self.variant == 4:
            header = ['x', 'y', 'z', 'tk1', 'tk2', 'tk3', 'tk4', 'tk', 'pr']
        # Print everything        
        col_width = max(len(s) for s in diff_trail.values()) + 2
        header_str = "Rounds\t"
        data_str = ""
        current_row = 0
        for entry in header[0:-2]:
            header_str += entry.ljust(col_width)
        header_str += header[-2].ljust(col_width//2)
        header_str += header[-1].ljust(7)
        for r in range(self.rounds + 1):
            data_str += str(current_row) + '\t'            
            data_str += diff_trail.get('x_' + str(r), 'none').ljust(col_width)
            data_str += diff_trail.get('y_' + str(r), 'none').ljust(col_width)
            data_str += diff_trail.get('z_' + str(r), 'none').ljust(col_width)
            if self.variant >= 1:
                data_str += diff_trail.get('tk1_' + str(r), 'none').ljust(col_width)
            if self.variant >= 2:
                data_str += diff_trail.get('tk2_' + str(r), 'none').ljust(col_width)
            if self.variant >= 3:
                data_str += diff_trail.get('tk3_' + str(r), 'none').ljust(col_width)
            if self.variant >= 4:
                data_str += diff_trail.get('tk4_' + str(r), 'none').ljust(col_width)
            data_str += diff_trail.get('tk_' + str(r), 'none').ljust(col_width//2)
            data_str += diff_trail.get('pr_' + str(r), 'none').ljust(7)
            data_str += '\n'
            current_row += 1
        print(header_str)
        print("-"*len(header_str))
        print(data_str)        
        print("Weight: " + '-' + str(self.total_weight))
        return

def loadparameters(args):
    '''
    Extract parameters from the argument list and input file
    '''

    # Load default values
    params = {"rounds" : 9,
              "variant" : 4,
              "upperbound1": None,
              "upperbound2": None,
              "start_round": 0,
              "end_round": None,
              "cellsize" : 4,
              "skipsb": 0,
              "mode" : 0,
              "sweight" : 0,
              "endweight" : 128,              
              "timelimit" : -1,
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

    # Override parameters if they are set on command line    
    if args.rounds:
        params["rounds"] = args.rounds[0]    
    
    if args.variant:
        params["variant"] = args.variant[0]

    if args.upperbound1:
        params["upperbound1"] = args.upperbound1[0]
    
    if args.upperbound2:
        params["upperbound2"] = args.upperbound2[0]
    
    if args.start_round:
        params["start_round"] = args.start_round[0]
    
    if args.end_round:
        params["end_round"] = args.end_round[0]
    
    if args.cellsize:
        params["cellsize"] = args.cellsize[0]

    if args.skipsb:
        params["skipsb"] = args.skipsb[0]
    
    if args.mode:
        params["mode"] = args.mode[0]

    if args.sweight:
        params["sweight"] = args.sweight[0]
    
    if args.endweight:
        params["endweight"] = args.endweight[0]
    
    if args.timelimit:
        params["timelimit"] = args.timelimit[0]

    return params

def main():
    '''
    Parse the arguments and start the request functionality with the provided
    parameters.
    '''
    parser = ArgumentParser(description="This tool finds the best differential"
                                        "trail in a cryptographic primitive"
                                        "using Gurobi",
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-i", "--inputfile", nargs=1, type=str,
                        help="Input file with parameters")
    parser.add_argument("-r", "--rounds", nargs=1, type=int,
                        help="Number of rounds")
    parser.add_argument("-v", "--variant", nargs=1, type=int,
                        help="variant of SKINNY (0, 1, 2, 3, 4)")
    parser.add_argument("-u1", "--upperbound1", nargs=1, type=int,
                        help="Upper bound on the weight of differential trail", default=None)
    parser.add_argument("-u2", "--upperbound2", nargs=1, type=int,
                        help="Upper bound on the weight of halfway differential trail", default=None)
    parser.add_argument("-str", "--start_round", nargs=1, type=int,
                        help="Start round", default=0)
    parser.add_argument("-enr", "--end_round", nargs=1, type=int,
                        help="End round", default=None)
    parser.add_argument("-c", "--cellsize", nargs=1, type=int,
                        help="cell size")
    parser.add_argument("-ssb", "--skipsb", nargs=1, type=int,
                        help="skip the 1st S-box layer", default=0)
    parser.add_argument('--mode', nargs=1, type=int, 
                        choices=[0, 1, 2], help=
                        "0 = search for the best differential characteristic\n"                        
                        "1 = search for multiple differential characteristics\n"
                        "2 = compute the differential effect")
    parser.add_argument("-sw", "--sweight", nargs=1, type=int,
                        help="starting weight for the trail search")
    parser.add_argument("-ew", "--endweight", nargs=1, type=int,
                        help="ending weight for the trail search")
    parser.add_argument("-t", "--timelimit", nargs=1, type=int,
                        help="time limit for the search")       

    # Parse command line arguments and construct parameter list
    args = parser.parse_args()
    params = loadparameters(args)
    skinny = Differential(params, True)
    skinny.make_model()
    skinny.solve()

if __name__ == "__main__":
    main()