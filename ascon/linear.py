#!/usr/env/bin python3
#-*- coding: UTF-8 -*-

"""
An Automatic Tool for Linear Analysis of Ascon
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

Disclaimer: We acknowledge that the Ascon permutation doesn't adhere to statistical assumptions 
in differential analysis, such as the random sub-key assumption (essentially, there's no sub-key) 
or Markov cipher assumption. The tool's primary function is to identify and verify impossible-differential, 
zero-correlation, and integral attacks on Ascon, not finding some differential/linear characteristics.
"""

from argparse import ArgumentParser, RawTextHelpFormatter
import yaml
import time
from gurobipy import read
from gurobipy import GRB
import math
import os


"""
Modeling the linear analysis of Ascon
This tool can:
 - find the best linear trail
 - find multiple linear trails
 - compute the squared correlation taking the clustering effect into account

MILP variables:

    x_RoundNumber_RowNumber_ColumnNumber
    y_RoundNumber_RowNumber_ColumnNumber

    x_round_row_0:  msb (left most bit)
    x_round_row_63: lsb (right most bit)
                                
We have used sboxanalyzer to model the differential/linear behavior of Ascon's S-box:
https://github.com//sboxanalyzer

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.04 seconds
Number of constraints: 96
Input:	a0||a1||a2||a3||a4; a0: msb
Output:	b0||b1||b2||b3||b4; b0: msb
Weight: 4.0000 p0 + 2.0000 p1
sage: milp
['- a0 - a3 - p1 >= -2',
 '- a1 - a3 - p1 >= -2',
 '- a2 + p0 + p1 >= 0',
 '- b2 + p0 + p1 >= 0',
 '- b3 + p0 + p1 >= 0',
 '- b4 + p0 + p1 >= 0',
 'a1 + a2 + a3 - p0 >= 0',
 'a0 + a1 + a4 - p0 >= 0',
 'a0 + a1 + b2 - p0 >= 0',
 'a2 - a3 + b2 + p0 >= 0',
 '- a0 + a4 + b3 + p0 >= 0',
 'a0 - a1 - a4 - p1 >= -2',
 '- a0 - a1 + a4 - p1 >= -2',
 '- a2 + a3 - b0 - p1 >= -2',
 'a0 - a4 - b2 - p1 >= -2',
 'a0 - a1 + a3 + a4 + b0 >= 0',
 'a0 + a1 - a2 - a3 - b2 >= -2',
 '- a0 - a1 - a4 - b1 - b2 >= -4',
 'a1 - a2 + a3 + b0 + b2 >= 0',
 '- a0 - a1 - a4 + b1 + b2 >= -2',
 '- a0 + a3 - a4 + b0 - p0 >= -2',
 'a0 + a3 - a4 + b3 + p0 >= 0',
 '- a0 + a1 + b0 - b2 - p1 >= -2',
 'a0 + a1 + b1 - b2 - p1 >= -1',
 '- a1 + a2 - a4 - b4 - p1 >= -3',
 '- a1 + a3 - a4 - b0 - b2 - b4 >= -4',
 'a1 + a3 + a4 - b0 - b2 - b4 >= -2',
 'a2 - b0 + b1 + b2 - b3 - b4 >= -2',
 '- a0 + a2 + a3 - a4 - b0 + b4 >= -2',
 'a0 + a2 + a3 + a4 - b0 + b4 >= 0',
 '- a1 + b0 + b1 + b2 + b3 + b4 >= 0',
 '- a1 + a2 + a3 - b0 - b4 - p0 >= -3',
 'a1 + a2 + b2 + b3 + b4 - p0 >= 0',
 'a1 + a3 - b0 - b3 - b4 + p0 >= -2',
 'a0 + a1 + a2 + a3 + a4 - p1 >= 0',
 'a1 + a2 + a4 - b0 + b2 - p1 >= -1',
 '- a0 - a4 + b0 + b2 - b4 - p1 >= -3',
 'a1 + a2 + b2 - b3 - b4 + p1 >= -1',
 '- a1 + a2 + a3 + b0 + b4 + p1 >= 0',
 'a3 + b0 - b1 + b2 + b4 + p1 >= 0',
 '- a1 - a2 - a3 - a4 - b1 - b2 + b3 >= -5',
 'a0 - a2 - a3 - a4 + b1 + b2 + b3 >= -2',
 'a0 + a2 + a4 + b1 - b2 + b3 - b4 >= -1',
 'a0 + a4 + b0 + b1 - b2 + b3 - b4 >= -1',
 'a0 + a2 - b0 - b1 + b2 + b3 - b4 >= -2',
 '- a0 - a1 - a2 + a3 - a4 + b1 + b4 >= -3',
 'a1 + a2 - a3 + a4 - b2 - b3 + b4 >= -2',
 'a0 + a2 + a4 + b1 - b2 - b3 + b4 >= -1',
 'a0 + a4 + b0 + b1 - b2 - b3 + b4 >= -1',
 'a1 + a2 + a3 + b0 + b4 + p0 - p1 >= 0',
 '- a1 - a2 + a3 + b0 - b2 - b4 + p1 >= -3',
 '- a1 - a2 + a3 - a4 + b2 + b4 + p1 >= -2',
 'a0 - a1 + a2 - a3 - a4 + b0 - b1 - b3 >= -4',
 'a0 - a1 + a2 - a3 - a4 - b0 + b1 - b3 >= -4',
 '- a0 - a1 + a2 - a3 + a4 + b0 - b1 - b4 >= -4',
 '- a0 - a1 - a3 + a4 + b0 - b1 + b2 - b4 >= -4',
 'a0 - a1 + a2 + a4 - b1 - b2 - b3 - b4 >= -4',
 'a0 - a3 + a4 + b0 - b1 + b2 + b3 - b4 >= -2',
 '- a0 - a1 + a2 - a3 - b0 - b1 - b2 + b4 >= -5',
 '- a0 - a2 + a4 - b0 + b1 - b2 - b3 + b4 >= -4',
 'a0 + a2 + a4 + b1 + b2 + b3 + b4 - p0 >= 0',
 'a0 + a4 + b0 - b1 - b2 - b3 - b4 + p1 >= -3',
 'a0 + a2 + a4 - b1 - b2 + b3 + b4 + p1 >= -1',
 'a0 + a4 + b0 - b1 - b2 + b3 + b4 + p1 >= -1',
 '- a1 + a2 - a3 + a4 - b0 - b1 + b2 - b3 + b4 >= -4',
 'a0 - a2 - a3 + a4 - b0 - b1 + b2 + b3 + b4 >= -3',
 '- a0 - a2 + a3 + a4 - b0 - b2 + b4 - p0 + p1 >= -4',
 'a0 - a1 - a2 - a3 - a4 + b1 - b2 - b3 - p0 + p1 >= -6',
 'a0 - a1 - a2 - a3 - a4 - b1 + b2 - b3 - p0 + p1 >= -6',
 'a0 - a1 + a2 - a3 - a4 - b0 - b1 + b3 - p0 + p1 >= -5',
 'a0 - a1 + a2 - a3 - a4 + b0 + b1 + b3 - p0 + p1 >= -3',
 '- a0 - a1 - a3 + a4 - b0 - b1 + b2 + b4 - p0 + p1 >= -5',
 '- a0 - a1 - a2 - a3 + a4 + b0 + b1 - b2 - b4 - p0 + p1 >= -6',
 '- a1 - a2 - a3 + a4 - b0 + b1 - b2 + b3 + b4 - p0 + p1 >= -5',
 'a0 - a2 - a3 + a4 - b0 - b1 - b2 - b3 + b4 >= -5',
 'a0 - a2 - a3 + a4 - b0 + b1 + b2 - b3 + b4 >= -3',
 '- a0 - a1 - a2 - a3 - b0 - b1 - b2 - b4 >= -7',
 '- a0 - a1 - a2 - a3 + b0 - b1 - b2 + b4 >= -5',
 '- a0 - a1 + a2 + a4 - b0 + b1 - b4 >= -3',
 'a1 - a2 - a3 + b0 + b3 + b4 >= -1',
 '- a2 - a3 - b0 + b1 + b2 + b3 - b4 >= -3',
 '- a0 - a1 + b0 + b1 + b2 + b4 >= -1',
 'a0 - a2 - a3 - b0 + b1 - b2 - b3 - b4 >= -5',
 'a0 - a2 - a3 - b0 - b1 + b2 - b3 - b4 >= -5',
 'a0 + a4 + b0 + b1 + b2 - b3 - b4 >= -1',
 'a1 - a2 - a3 - b0 - b3 + b4 >= -3',
 'a1 - a2 - a3 + b0 - b3 - b4 >= -3',
 '- a2 - a3 - b0 - b1 - b2 + b3 - b4 >= -5',
 '- a0 + a1 + a2 - a4 + b3 + b4 >= -1',
 '- a0 - a1 + a2 + a4 + b0 + b1 + b4 >= -1',
 '- a0 + a1 + a2 - a4 - b3 - b4 >= -3',
 'a0 + b0 - b1 + b2 - b3 + b4 >= -1',
 'a1 - a2 - a3 - b0 + b3 - b4 >= -3',
 '- a0 - a1 + a3 - b0 + b2 - b4 >= -3',
 '- a0 - a1 - b0 + b1 + b2 - b4 >= -3',
 'a1 + a2 + a4 - b2 + b3 - b4 >= -1']
"""


class Differential:
    '''
    Convert the differential analysis of SKINNY to an MILP problem
    '''

    count = 0
    def __init__(self, param, exact=True):
        self.booster = False
        self.no_rounds = param['rounds']        
        self.time_limit = param['timelimit']            
        self.fixed_variables = param['fixedVariables']
        self.mode = param['mode']
        self.start_weight = param['sweight']
        self.end_weight = param['endweight']
        self.eps = 1e-3
        self.exact = exact #A Boolean variable indicating whether we model DDT or *-DDT        
        self.total_weight = None                
        self.obj_func = ''
        self.used_variables = []
        self.rotation = [[0, 19, 28], 
                         [0, 61, 39], 
                         [0, 1, 6],
                         [0, 10, 17], 
                         [0, 7, 41]]
        self.pr_weights = [4, 2]
        self.model_filename = f"Ascon-{self.no_rounds}r.lp"

        #######################################################################################################
        #######################################################################################################
        #######################################################################################################
        #  __  __             _        _   _    _             _         _   _____ 
        # |  \/  |  ___    __| |  ___ | | | |_ | |__    ___  | |       / \ |_   _|
        # | |\/| | / _ \  / _` | / _ \| | | __|| '_ \  / _ \ | |      / _ \  | |  
        # | |  | || (_) || (_| ||  __/| | | |_ | | | ||  __/ | |___  / ___ \ | |  
        # |_|  |_| \___/  \__,_| \___||_|  \__||_| |_| \___| |_____|/_/   \_\|_|  
                                                                            
        self.sbox_exact_model = ['- a0 - a3 - p1 >= -2',
                                '- a1 - a3 - p1 >= -2',
                                '- a2 + p0 + p1 >= 0',
                                '- b2 + p0 + p1 >= 0',
                                '- b3 + p0 + p1 >= 0',
                                '- b4 + p0 + p1 >= 0',
                                'a1 + a2 + a3 - p0 >= 0',
                                'a0 + a1 + a4 - p0 >= 0',
                                'a0 + a1 + b2 - p0 >= 0',
                                'a2 - a3 + b2 + p0 >= 0',
                                '- a0 + a4 + b3 + p0 >= 0',
                                'a0 - a1 - a4 - p1 >= -2',
                                '- a0 - a1 + a4 - p1 >= -2',
                                '- a2 + a3 - b0 - p1 >= -2',
                                'a0 - a4 - b2 - p1 >= -2',
                                'a0 - a1 + a3 + a4 + b0 >= 0',
                                'a0 + a1 - a2 - a3 - b2 >= -2',
                                '- a0 - a1 - a4 - b1 - b2 >= -4',
                                'a1 - a2 + a3 + b0 + b2 >= 0',
                                '- a0 - a1 - a4 + b1 + b2 >= -2',
                                '- a0 + a3 - a4 + b0 - p0 >= -2',
                                'a0 + a3 - a4 + b3 + p0 >= 0',
                                '- a0 + a1 + b0 - b2 - p1 >= -2',
                                'a0 + a1 + b1 - b2 - p1 >= -1',
                                '- a1 + a2 - a4 - b4 - p1 >= -3',
                                '- a1 + a3 - a4 - b0 - b2 - b4 >= -4',
                                'a1 + a3 + a4 - b0 - b2 - b4 >= -2',
                                'a2 - b0 + b1 + b2 - b3 - b4 >= -2',
                                '- a0 + a2 + a3 - a4 - b0 + b4 >= -2',
                                'a0 + a2 + a3 + a4 - b0 + b4 >= 0',
                                '- a1 + b0 + b1 + b2 + b3 + b4 >= 0',
                                '- a1 + a2 + a3 - b0 - b4 - p0 >= -3',
                                'a1 + a2 + b2 + b3 + b4 - p0 >= 0',
                                'a1 + a3 - b0 - b3 - b4 + p0 >= -2',
                                'a0 + a1 + a2 + a3 + a4 - p1 >= 0',
                                'a1 + a2 + a4 - b0 + b2 - p1 >= -1',
                                '- a0 - a4 + b0 + b2 - b4 - p1 >= -3',
                                'a1 + a2 + b2 - b3 - b4 + p1 >= -1',
                                '- a1 + a2 + a3 + b0 + b4 + p1 >= 0',
                                'a3 + b0 - b1 + b2 + b4 + p1 >= 0',
                                '- a1 - a2 - a3 - a4 - b1 - b2 + b3 >= -5',
                                'a0 - a2 - a3 - a4 + b1 + b2 + b3 >= -2',
                                'a0 + a2 + a4 + b1 - b2 + b3 - b4 >= -1',
                                'a0 + a4 + b0 + b1 - b2 + b3 - b4 >= -1',
                                'a0 + a2 - b0 - b1 + b2 + b3 - b4 >= -2',
                                '- a0 - a1 - a2 + a3 - a4 + b1 + b4 >= -3',
                                'a1 + a2 - a3 + a4 - b2 - b3 + b4 >= -2',
                                'a0 + a2 + a4 + b1 - b2 - b3 + b4 >= -1',
                                'a0 + a4 + b0 + b1 - b2 - b3 + b4 >= -1',
                                'a1 + a2 + a3 + b0 + b4 + p0 - p1 >= 0',
                                '- a1 - a2 + a3 + b0 - b2 - b4 + p1 >= -3',
                                '- a1 - a2 + a3 - a4 + b2 + b4 + p1 >= -2',
                                'a0 - a1 + a2 - a3 - a4 + b0 - b1 - b3 >= -4',
                                'a0 - a1 + a2 - a3 - a4 - b0 + b1 - b3 >= -4',
                                '- a0 - a1 + a2 - a3 + a4 + b0 - b1 - b4 >= -4',
                                '- a0 - a1 - a3 + a4 + b0 - b1 + b2 - b4 >= -4',
                                'a0 - a1 + a2 + a4 - b1 - b2 - b3 - b4 >= -4',
                                'a0 - a3 + a4 + b0 - b1 + b2 + b3 - b4 >= -2',
                                '- a0 - a1 + a2 - a3 - b0 - b1 - b2 + b4 >= -5',
                                '- a0 - a2 + a4 - b0 + b1 - b2 - b3 + b4 >= -4',
                                'a0 + a2 + a4 + b1 + b2 + b3 + b4 - p0 >= 0',
                                'a0 + a4 + b0 - b1 - b2 - b3 - b4 + p1 >= -3',
                                'a0 + a2 + a4 - b1 - b2 + b3 + b4 + p1 >= -1',
                                'a0 + a4 + b0 - b1 - b2 + b3 + b4 + p1 >= -1',
                                '- a1 + a2 - a3 + a4 - b0 - b1 + b2 - b3 + b4 >= -4',
                                'a0 - a2 - a3 + a4 - b0 - b1 + b2 + b3 + b4 >= -3',
                                '- a0 - a2 + a3 + a4 - b0 - b2 + b4 - p0 + p1 >= -4',
                                'a0 - a1 - a2 - a3 - a4 + b1 - b2 - b3 - p0 + p1 >= -6',
                                'a0 - a1 - a2 - a3 - a4 - b1 + b2 - b3 - p0 + p1 >= -6',
                                'a0 - a1 + a2 - a3 - a4 - b0 - b1 + b3 - p0 + p1 >= -5',
                                'a0 - a1 + a2 - a3 - a4 + b0 + b1 + b3 - p0 + p1 >= -3',
                                '- a0 - a1 - a3 + a4 - b0 - b1 + b2 + b4 - p0 + p1 >= -5',
                                '- a0 - a1 - a2 - a3 + a4 + b0 + b1 - b2 - b4 - p0 + p1 >= -6',
                                '- a1 - a2 - a3 + a4 - b0 + b1 - b2 + b3 + b4 - p0 + p1 >= -5',
                                'a0 - a2 - a3 + a4 - b0 - b1 - b2 - b3 + b4 >= -5',
                                'a0 - a2 - a3 + a4 - b0 + b1 + b2 - b3 + b4 >= -3',
                                '- a0 - a1 - a2 - a3 - b0 - b1 - b2 - b4 >= -7',
                                '- a0 - a1 - a2 - a3 + b0 - b1 - b2 + b4 >= -5',
                                '- a0 - a1 + a2 + a4 - b0 + b1 - b4 >= -3',
                                'a1 - a2 - a3 + b0 + b3 + b4 >= -1',
                                '- a2 - a3 - b0 + b1 + b2 + b3 - b4 >= -3',
                                '- a0 - a1 + b0 + b1 + b2 + b4 >= -1',
                                'a0 - a2 - a3 - b0 + b1 - b2 - b3 - b4 >= -5',
                                'a0 - a2 - a3 - b0 - b1 + b2 - b3 - b4 >= -5',
                                'a0 + a4 + b0 + b1 + b2 - b3 - b4 >= -1',
                                'a1 - a2 - a3 - b0 - b3 + b4 >= -3',
                                'a1 - a2 - a3 + b0 - b3 - b4 >= -3',
                                '- a2 - a3 - b0 - b1 - b2 + b3 - b4 >= -5',
                                '- a0 + a1 + a2 - a4 + b3 + b4 >= -1',
                                '- a0 - a1 + a2 + a4 + b0 + b1 + b4 >= -1',
                                '- a0 + a1 + a2 - a4 - b3 - b4 >= -3',
                                'a0 + b0 - b1 + b2 - b3 + b4 >= -1',
                                'a1 - a2 - a3 - b0 + b3 - b4 >= -3',
                                '- a0 - a1 + a3 - b0 + b2 - b4 >= -3',
                                '- a0 - a1 - b0 + b1 + b2 - b4 >= -3',
                                'a1 + a2 + a4 - b2 + b3 - b4 >= -1']
    
        #######################################################################################################
        #######################################################################################################
        #######################################################################################################
        #  __  __             _        _   _    _             ____   _                _     ____   ____  _____ 
        # |  \/  |  ___    __| |  ___ | | | |_ | |__    ___  / ___| | |_  __ _  _ __ | |_  |  _ \ |  _ \|_   _|
        # | |\/| | / _ \  / _` | / _ \| | | __|| '_ \  / _ \ \___ \ | __|/ _` || '__|| __| | | | || | | | | |  
        # | |  | || (_) || (_| ||  __/| | | |_ | | | ||  __/  ___) || |_| (_| || |   | |_  | |_| || |_| | | |  
        # |_|  |_| \___/  \__,_| \___||_|  \__||_| |_| \___| |____/  \__|\__,_||_|    \__| |____/ |____/  |_|  
                                        
        self.sbox_star_model = ['- a0 - a1 + a3 - a4 - b0 >= -3',
                                'a0 + a1 - a2 - a3 - b2 >= -2',
                                '- a0 - a1 - a4 - b1 - b2 >= -4',
                                'a0 + a1 + a4 + b1 - b2 >= 0',
                                'a0 + a1 + a2 - a3 + b2 >= 0',
                                'a0 + a1 - a2 + a3 + b2 >= 0',
                                'a1 - a2 + a3 + b0 + b2 >= 0',
                                '- a0 - a1 - a4 + b1 + b2 >= -2',
                                'a1 - a2 + a3 + a4 - b0 - b2 >= -2',
                                'a0 + a1 + a3 + a4 + b2 - b3 >= 0',
                                'a0 + a2 + a3 - a4 - b0 - b4 >= -2',
                                '- a0 + a2 + a3 + a4 - b0 - b4 >= -2',
                                '- a0 + a2 + a3 - a4 + b0 - b4 >= -2',
                                '- a0 + a3 - a4 + b0 + b2 - b4 >= -2',
                                'a1 + a2 + a3 - b2 - b3 - b4 >= -2',
                                'a1 + a2 + a3 + a4 + b3 - b4 >= 0',
                                'a1 + a2 + a3 + b0 + b3 - b4 >= 0',
                                'a0 + a2 + a3 + a4 - b0 + b4 >= 0',
                                'a0 + a2 + a3 - a4 + b0 + b4 >= 0',
                                'a1 + a2 + a3 - a4 - b2 + b4 >= -1',
                                'a2 + a3 + a4 + b0 - b2 + b4 >= 0',
                                '- a0 + a1 + a2 + a3 + b2 + b4 >= 0',
                                'a0 - a2 + a3 - a4 + b2 + b4 >= -1',
                                'a3 + a4 + b0 - b1 + b2 + b4 >= 0',
                                'a1 + a2 + a3 - b0 + b3 + b4 >= 0',
                                '- a0 + a1 + a2 + b2 + b3 + b4 >= 0',
                                '- a3 + b0 + b1 + b2 + b3 + b4 >= 0',
                                '- a1 - a2 - a3 - a4 - b1 - b2 + b3 >= -5',
                                'a0 - a2 - a3 - a4 + b1 + b2 + b3 >= -2',
                                'a0 - a1 - a2 + a3 - a4 - b2 - b4 >= -4',
                                '- a1 - a2 + a3 + b0 - b1 - b2 - b4 >= -4',
                                'a0 + a2 + a4 + b1 + b2 - b3 - b4 >= -1',
                                'a2 + a4 - b0 + b1 - b2 + b3 - b4 >= -2',
                                'a0 + a2 + a4 - b1 + b2 + b3 - b4 >= -1',
                                'a0 + a4 + b0 - b1 + b2 + b3 - b4 >= -1',
                                '- a0 - a2 + a3 + a4 - b0 - b2 + b4 >= -3',
                                '- a0 - a2 + a3 - a4 + b0 - b2 + b4 >= -3',
                                'a1 + a2 - a3 + a4 - b2 - b3 + b4 >= -2',
                                'a0 + a2 + a4 + b1 - b2 - b3 + b4 >= -1',
                                'a0 + a4 + b0 + b1 - b2 - b3 + b4 >= -1',
                                'a0 - a1 + a2 - a3 - a4 + b0 - b1 - b3 >= -4',
                                'a0 - a1 + a2 - a3 - a4 - b0 + b1 - b3 >= -4',
                                'a0 - a1 - a2 - a3 - a4 + b1 - b2 - b3 >= -5',
                                'a0 - a1 - a2 - a3 - a4 - b1 + b2 - b3 >= -5',
                                'a0 - a1 + a2 - a3 - a4 - b0 - b1 + b3 >= -4',
                                'a0 - a1 + a2 - a3 - a4 + b0 + b1 + b3 >= -2',
                                '- a0 - a1 + a2 - a3 + b0 - b1 - b2 - b4 >= -5',
                                '- a0 - a1 - a2 + a4 + b0 + b1 - b2 - b4 >= -4',
                                '- a0 - a1 - a3 + a4 + b0 - b1 + b2 - b4 >= -4',
                                '- a0 + a1 + a2 - a3 + a4 + b2 - b3 - b4 >= -3',
                                '- a0 - a1 + a2 - a3 - b0 - b1 - b2 + b4 >= -5',
                                '- a0 - a1 - a3 + a4 - b0 - b1 + b2 + b4 >= -4',
                                '- a0 - a2 + a4 - b0 + b1 - b2 - b3 + b4 >= -4',
                                '- a1 + a2 - a3 - b0 - b1 - b2 + b3 + b4 >= -4',
                                'a0 + a2 + a4 - b0 + b1 + b2 + b3 + b4 >= 0',
                                'a0 - a2 - a3 + a4 + b0 - b1 - b2 - b3 - b4 >= -5',
                                'a0 - a1 - a3 + a4 + b0 + b1 - b2 + b3 - b4 >= -3',
                                '- a1 + a2 - a3 + a4 - b0 - b1 + b2 - b3 + b4 >= -4',
                                'a0 - a1 - a3 + a4 + b0 - b1 - b2 + b3 + b4 >= -3',
                                '- a1 - a2 - a3 + a4 - b0 + b1 - b2 + b3 + b4 >= -4',
                                'a0 - a2 - a3 + a4 - b0 - b1 + b2 + b3 + b4 >= -3',
                                'a0 - a2 - a3 + a4 - b0 + b1 + b2 - b3 + b4 >= -3',
                                'a0 - a2 - a3 + a4 - b0 - b1 - b2 - b3 + b4 >= -5',
                                'a0 - a1 + a2 + a4 - b1 - b2 - b3 - b4 >= -4',
                                '- a0 - a1 - a2 - a3 - b0 - b1 - b2 - b4 >= -7',
                                '- a0 - a1 - a2 - a3 + b0 - b1 - b2 + b4 >= -5',
                                '- a0 - a1 + b0 + b1 + b2 + b4 >= -1',
                                '- a0 - a1 + a2 + a4 - b0 + b1 - b4 >= -3',
                                'a1 - a2 - a3 + b0 + b3 + b4 >= -1',
                                '- a2 - a3 - b0 + b1 + b2 + b3 - b4 >= -3',
                                'a0 - a2 - a3 - b0 + b1 - b2 - b3 - b4 >= -5',
                                'a0 - a2 - a3 - b0 - b1 + b2 - b3 - b4 >= -5',
                                'a0 + a4 + b0 + b1 + b2 - b3 - b4 >= -1',
                                'a1 - a2 - a3 - b0 - b3 + b4 >= -3',
                                '- a2 - a3 - b0 - b1 - b2 + b3 - b4 >= -5',
                                '- a0 + a1 + a2 - a4 + b3 + b4 >= -1',
                                'a1 - a2 - a3 + b0 - b3 - b4 >= -3',
                                '- a0 - a1 + a2 + a4 + b0 + b1 + b4 >= -1',
                                '- a0 + a1 + a2 - a4 - b3 - b4 >= -3',
                                'a0 + b0 - b1 + b2 - b3 + b4 >= -1',
                                'a1 - a2 - a3 - b0 + b3 - b4 >= -3',
                                '- a0 - a1 + a3 - b0 + b2 - b4 >= -3',
                                '- a0 - a1 - b0 + b1 + b2 - b4 >= -3',
                                '- a0 + a1 + a3 - a4 + b0 >= -1',
                                'a1 + a2 + a4 - b2 + b3 - b4 >= -1',
                                'a0 - a1 + a3 + a4 + b0 >= 0']
    
    def create_objective_function(self):
        '''
        Create the objective function
        '''
        
        minus_log2_p = []
        if self.exact:
            for round in range(self.no_rounds):
                for column in range(64):
                    minus_log2_p += [f"4 pr4_{round}_{column} + 2 pr2_{round}_{column}"]
            lp_contents = ' + '.join(minus_log2_p)
        else:
            lp_contents = "0"        
        return lp_contents

    def create_state_variables(self, r, s):
        '''
        Generate the state variables
        '''

        array = [['' for _ in range(64)] for _ in range(5)]
        for i in range(0, 5):
            for j in range(0, 64):
                array[i][j] = f"{s}_{r}_{i}_{j}"
                self.used_variables.append(array[i][j])
        return array

    def flatten(self, state_array):
        '''
        Get a state array and output a flatten list
        '''

        flat_list = []
        for frame in range(len(state_array)):
            for bit_number in range(len(state_array[0])):
                flat_list.append(state_array[frame][bit_number])
        return flat_list
    
    def create_probability_variables(self, r):
        '''
        Generate the variables corresponding to differential probabilities
        '''
        
        array = [['' for _ in range(3)] for _ in range(64)]
        for col in range(64):
            array[col] = [f"pr4_{r}_{col}", f"pr2_{r}_{col}"]
            if not self.booster:
                self.used_variables += array[col]
        return array

    def xor3(self, b, a0, a1, a2):
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

    def linear_layer(self, x, y):
        '''
        Model the inverse of MixColumns in linear analysis
        '''

        lp_contents = ""
        for row in range(5):
            for col in range(64):
                lp_contents += self.xor3(y[row][col], 
                                         x[row][(col + self.rotation[row][0])%64], 
                                         x[row][(col + self.rotation[row][1])%64], 
                                         x[row][(col + self.rotation[row][2])%64])
        return lp_contents
    
    def subcells_exact(self, x, y, p):
        '''
        Model the 5-bit S-box of Ascon (exact)
        '''

        lp_contents = ""      
        
        for col in range(64):
            for ineq in self.sbox_exact_model:
                for row in range(5):
                    ineq = ineq.replace(f"a{row}", x[row][col])
                    ineq = ineq.replace(f"b{row}", y[row][col])
                for i in range(2):
                    ineq = ineq.replace(f"p{i}", p[col][i])
                lp_contents += ineq + "\n"            
        return lp_contents
    
    def subcells_star(self, x, y):
        '''
        Model the 5-bit S-box of Ascon (*-DDT)
        '''

        lp_contents = ""        
        for col in range(64):
            for ineq in self.sbox_star_model:
                for row in range(5):
                    ineq = ineq.replace(f"a{row}", x[row][col])
                    ineq = ineq.replace(f"b{row}", y[row][col])
                lp_contents += ineq + "\n"            
        return lp_contents

    def ascon_permutation(self):
        '''
        Generate the MILP constraints modeling the propagation of differences through the permutation
        '''
        
        lp_contents = ""
        for r in range(self.no_rounds):
            x = self.create_state_variables(r, 'x')
            y = self.create_state_variables(r, 'y')
            lp_contents += self.linear_layer(x, y)
            x_next = self.create_state_variables(r + 1, 'x')
            if self.exact:
                p = self.create_probability_variables(r)
                lp_contents += self.subcells_exact(y, x_next, p)
            else:
                lp_contents += self.subcells_star(y, x_next)
        return lp_contents
    
    def exclude_trivial_trail(self):
        lp_contents = ""
        x = self.create_state_variables(0, 'x')
        temp = self.flatten(x)
        lp_contents += " + ".join(temp) + " >= 1\n"
        return lp_contents
    
    def declare_fixed_variables(self):
        lp_contents = ""
        for cond in self.fixed_variables.items():            
            var = cond[0]
            val = cond[1]
            var = var.split('_')
            assert(len(var) == 3)
            state_vars = [f"{var[0]}_{var[1]}_{var[2]}_{i}" for i in range(64)]
            for i in range(64):
                if val[i] != "?":
                    lp_contents += f"{state_vars[i]} = {val[i]}\n"
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
        lp_contents += self.exclude_trivial_trail()
        lp_contents += self.ascon_permutation()
        lp_contents += self.declare_fixed_variables()
        lp_contents += self.declare_variables_type()
        if os.path.exists(self.model_filename):
            os.remove(self.model_filename)
        with open(self.model_filename, 'w') as fileobj:
            fileobj.write(lp_contents)
        print(f"MILP model was written into {self.model_filename}\n")  

    
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
            # self.compute_differential_effect_classic_method()
        else:
            print("mode should be in [0, 1, 2]")
        os.remove(self.model_filename)
        return status

    def parse_solver_output(self):
        '''
        Extract the linear characteristic from the solver output
        '''

        characteristic = dict()
        for r in range(self.no_rounds + 1):
            x = self.create_state_variables(r, 'x')
            for row in range(0, 5):                
                x_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(self.model.getVarByName(t).X)), x[row]))), 2))[2:].zfill(16)
                characteristic[f"x_{r}_{row}"] = x_value
        for r in range(self.no_rounds):
            y = self.create_state_variables(r, 'y')
            for row in range(0, 5):
                y_value = hex(int('0b' + ''.join(list(map(lambda t: str(int(self.model.getVarByName(t).X)), y[row]))), 2))[2:].zfill(16)
                characteristic[f"y_{r}_{row}"] = y_value      
            round_probability = 0
            if self.exact:
                p = self.create_probability_variables(r)
                for col in range(64):
                    for i in range(2):
                        round_probability += float(self.pr_weights[i]) * int(self.model.getVarByName(p[col][i]).X)
                characteristic['pr_' +str(r)] = '-' + str(round_probability)
        return characteristic
    
    def print_trail(self, diff_trail):
        '''
        Print out the obtained linear trail
        '''
        
        if self.exact:
            header = ['x', 'y', 'x', 'pr']
        else:
            header = ['x', 'y']        
        # Print everything
        col_width = max(len(s) for s in diff_trail.values()) + 2
        header_str = "Rounds\t"
        data_str = ""
        current_row = 0
        for entry in header[0:-2]:
            header_str += entry.ljust(col_width)
        header_str += header[-2].ljust(col_width)
        header_str += header[-1].ljust(7)
        for r in range(self.no_rounds + 1):
            for row in range(5):
                data_str += str(current_row) + '\t'            
                data_str += diff_trail.get(f"x_{r}_{row}", 'none').ljust(col_width)
                data_str += diff_trail.get(f"y_{r}_{row}", 'none').ljust(col_width)
                data_str += diff_trail.get(f"x_{r + 1}_{row}", 'none').ljust(col_width)
                if row == 4 and self.exact:
                    data_str += diff_trail.get('pr_' + str(r), 'none').ljust(7)                
                data_str += "\n"
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
    params = {"rounds" : 1,
              "mode" : 0,
              "sweight" : 0,
              "endweight" : 1000,
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