#!/usr/env/bin python3
#-*- coding: UTF-8 -*-

"""
MIT License

Copyright (c) 2024 

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

import logging
from pathlib import Path
from random import randint
logging.basicConfig(filename="minizinc-python.log", level=logging.DEBUG)
import time
import minizinc
import datetime
from argparse import ArgumentParser, RawTextHelpFormatter
from diff import Diff
from lin import Lin
from draw import *
import subprocess
# Check if "OR Tools" appears in the output of "minizinc --solvers" command 
try:
    output = subprocess.run(['minizinc', '--solvers'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if "com.google.ortools.sat" in output.stdout.decode("utf-8"):
        ortools_available = True
        print("OR Tools is available")
    else:
        ortools_available = False
        print("OR Tools is not available")
except FileNotFoundError:
    ortools_available = False
    print("OR Tools is not available")


class DiffLin:
    DL_counter = 0

    def __init__(self, param) -> None:
        DiffLin.DL_counter += 1
        self.id = DiffLin.DL_counter
        self.name = "DiffLin" + str(self.id)
        self.type = "DiffLin"
        self.blocksize = param["blocksize"]
        self.halfblocksize = self.blocksize // 2
        self.RU = param["RU"]
        self.RM = param["RM"]
        self.RL = param["RL"]
        self.RMU = param["RMU"]
        self.RML = param["RML"]
        self.WU = param["WU"]
        self.WM = param["WM"]
        self.WL = param["WL"]
        self.RD = self.RU + self.RM + self.RL
        self.cp_solver_name = param["solver"]
        ##################################################
        if ortools_available:
            if self.cp_solver_name == "ortools":
                self.cp_solver_name = "com.google.ortools.sat"
        #################################################
        self.cp_solver = minizinc.Solver.lookup(self.cp_solver_name)
        self.time_limit = param["timelimit"]
        self.num_of_threads = param["np"]
        self.mzn_file_name = None
        self.output_file_name = param["output"]
        self.mzn_file_name = "attack.mzn"
    
    #############################################################################################################################################
    #############################################################################################################################################    
    #  ____                           _        __                        ____   _       _    _                       _       _                 
    # / ___|   ___   __ _  _ __  ___ | |__    / _|  ___   _ __    __ _  |  _ \ (_) ___ | |_ (_) _ __    __ _  _   _ (_) ___ | |__    ___  _ __ 
    # \___ \  / _ \ / _` || '__|/ __|| '_ \  | |_  / _ \ | '__|  / _` | | | | || |/ __|| __|| || '_ \  / _` || | | || |/ __|| '_ \  / _ \| '__|
    #  ___) ||  __/| (_| || |  | (__ | | | | |  _|| (_) || |    | (_| | | |_| || |\__ \| |_ | || | | || (_| || |_| || |\__ \| | | ||  __/| |   
    # |____/  \___| \__,_||_|   \___||_| |_| |_|   \___/ |_|     \__,_| |____/ |_||___/ \__||_||_| |_| \__, | \__,_||_||___/|_| |_| \___||_|   
    #                                                                                                  |___/                                   
    # Search for a distinguisher using MiniZinc

    def search(self):
        """
        Search for a distinguisher
        """

        if self.time_limit != -1:
            time_limit = datetime.timedelta(seconds=self.time_limit)
        else:
            time_limit = None
    
        start_time = time.time()
        #############################################################################################################################################
        print(f"Searching a distinguisher for {self.RD} rounds of SIMECK-{self.blocksize} ...")
        self.cp_model = minizinc.Model()
        self.cp_model.add_file(self.mzn_file_name)
        self.cp_inst = minizinc.Instance(solver=self.cp_solver, model=self.cp_model)
        self.cp_inst["blocksize"] = self.blocksize
        self.cp_inst["RU"] = self.RU
        self.cp_inst["RM"] = self.RM
        self.cp_inst["RL"] = self.RL
        self.cp_inst["RMU"] = self.RMU
        self.cp_inst["RML"] = self.RML
        self.cp_inst["WU"] = self.WU
        self.cp_inst["WM"] = self.WM
        self.cp_inst["WL"] = self.WL
        self.result = self.cp_inst.solve(timeout=time_limit, 
                                         processes=self.num_of_threads, 
                                         verbose=False, 
                                         debug_output=Path("./debug_output.txt",
                                         intermediate_solutions=True),
                                         random_seed=randint(0, 100),
                                         optimisation_level=2)
        #############################################################################################################################################
        elapsed_time = time.time() - start_time
        print(f"Solver status: {self.result.status}")
        print("Time used to find a distinguisher: {:0.02f} seconds".format(elapsed_time))
        if minizinc.Status.has_solution(self.result.status) or self.result.status == minizinc.Status.ERROR:
            self.attack_summary, self.upper_trail, self.lower_trail = self.parse_solution()
            self.attack_summary += "Time used to find a distinguisher: {:0.02f} seconds\n".format(elapsed_time)
            self.attack_summary += f"Solver status: {self.result.status}\n"
            self.diff_effect = self.compute_differential_effect(self.upper_trail)
            self.attack_summary += "#"*55 + "\n"
            self.attack_summary += "-Log2(PU)               = \t{:02d}\n".format(self.result["PU"])
            self.attack_summary += "-Log2(CM)              ~= \t{:02d}\n".format(self.result["CM"])
            self.attack_summary += "-Log2(CL)               = \t{:02d}\n".format(self.result["CL"])
            self.attack_summary += f"Differential effect = 2^({self.diff_effect})"
            print(self.attack_summary)
            draw = Draw(self, output_file_name=self.output_file_name, attack_summary=self.attack_summary)
            draw.generate_distinguisher_shape()
        elif self.result.status == minizinc.Status.UNSATISFIABLE:
            print("Model is unsatisfiable") 
        elif self.result.status == minizinc.Status.UNKNOWN:
            print("Unknown error!")
        else:
            print("Solving process was interrupted")

    #############################################################################################################################################
    #############################################################################################################################################
    #  ____                           _    _             ____          _         _    _               
    # |  _ \  __ _  _ __  ___   ___  | |_ | |__    ___  / ___|   ___  | | _   _ | |_ (_)  ___   _ __  
    # | |_) |/ _` || '__|/ __| / _ \ | __|| '_ \  / _ \ \___ \  / _ \ | || | | || __|| | / _ \ | '_ \ 
    # |  __/| (_| || |   \__ \|  __/ | |_ | | | ||  __/  ___) || (_) || || |_| || |_ | || (_) || | | |
    # |_|    \__,_||_|   |___/ \___|  \__||_| |_| \___| |____/  \___/ |_| \__,_| \__||_| \___/ |_| |_|
    # Parse the solution and print the distinguisher's specifications

    def parse_solution(self):
        """
        Parse the solution and print the distinguisher's specifications
        """
        
        upper_trail = {"xl": [[0 for _ in range(self.halfblocksize)] for _ in range(self.RU + self.RM + 1)],
                       "xr": [[0 for _ in range(self.halfblocksize)] for _ in range(self.RU + self.RM + 1)],}
        for rn in range(self.RU):
            upper_trail["xl"][rn] = self.result["xu_left"][rn]
            upper_trail["xr"][rn] = self.result["xu_right"][rn]

        for rn in range(self.RU, self.RU + self.RM + 1):
            upper_trail["xl"][rn] = self.result["xmu_left"][rn - self.RU]
            upper_trail["xr"][rn] = self.result["xmu_right"][rn - self.RU]
        
        lower_trail = {"xl": [[0 for _ in range(self.halfblocksize)] for _ in range(self.RM + self.RL + 1)],
                       "xr": [[0 for _ in range(self.halfblocksize)] for _ in range(self.RM + self.RL + 1)],}
        for rn in range(self.RM):
                lower_trail["xl"][rn] = self.result["xml_left"][rn]
                lower_trail["xr"][rn] = self.result["xml_right"][rn]
        for rn in range(self.RM, self.RM + self.RL + 1):
            lower_trail["xl"][rn] = self.result["xl_left"][rn - self.RM]
            lower_trail["xr"][rn] = self.result["xl_right"][rn - self.RM]

        input_diff = f"const char *DP_STR = \"" + "".join(map(str, upper_trail["xl"][0] + upper_trail["xr"][0])) + "\";\n"
        input_diff_middle = f"const char *DP_STR = \"" + "".join(map(str, upper_trail["xl"][self.RU] + upper_trail["xr"][self.RU])) + "\";\n"
        output_mask_middle = f"const char *LC_STR = \"" + "".join(map(str, lower_trail["xl"][self.RM] + lower_trail["xr"][self.RM])) + "\";\n"
        output_mask = f"const char *LC_STR = \"" + "".join(map(str, lower_trail["xl"][self.RM + self.RL] + lower_trail["xr"][self.RM + self.RL])) + "\";\n"
        
        attack_summary = f"Attack summary:\n"
        attack_summary += f"Setting: RU: {self.RU}, RM: {self.RM}, RL: {self.RL}, RMU: {self.RMU}, RML: {self.RML}\n"
        attack_summary += f"Weight factors: WU: {self.WU}, WM: {self.WM}, WL: {self.WL}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff.: \n{input_diff}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff. middle: \n{input_diff_middle}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask middle: \n{output_mask_middle}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask: \n{output_mask}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"PU:  {self.result['PU']}\n"
        attack_summary += f"CM:  {self.result['CM']}\n"
        attack_summary += f"CL:  {self.result['CL']}\n"
        attack_summary += f"Number of effective bit-positions in the middle: {self.result['CM']}\n"
        attack_summary += "#"*50 + "\n"
        # print the upper trail
        attack_summary += "Upper trail:\n"
        for rn in range(self.RU + self.RM + 1):
            attack_summary += "Round {:02d}: ".format(rn)
            attack_summary += "".join(map(str, upper_trail["xl"][rn] + upper_trail["xr"][rn])).replace("-1", "*") + "\n"
        attack_summary += "#"*50 + "\n"
        # print the lower trail
        attack_summary += "Lower trail: \n"
        for rn in range(self.RM + self.RL + 1):
            attack_summary += "Round {:02d}: ".format(rn)
            attack_summary += "".join(map(str, lower_trail["xl"][rn] + lower_trail["xr"][rn])).replace("-1", "*") + "\n"
        return attack_summary, upper_trail, lower_trail

    #############################################################################################################################################
    #############################################################################################################################################
    #############################################################################################################################################
    #  _____ _           _    ____                          _         ____  _  __  __                     _   _       _   _____          _ _     
    # |  ___(_)_ __   __| |  / ___|___  _ __   ___ _ __ ___| |_ ___  |  _ \(_)/ _|/ _| ___ _ __ ___ _ __ | |_(_) __ _| | |_   _| __ __ _(_) |___ 
    # | |_  | | '_ \ / _` | | |   / _ \| '_ \ / __| '__/ _ \ __/ _ \ | | | | | |_| |_ / _ \ '__/ _ \ '_ \| __| |/ _` | |   | || '__/ _` | | / __|
    # |  _| | | | | | (_| | | |__| (_) | | | | (__| | |  __/ ||  __/ | |_| | |  _|  _|  __/ | |  __/ | | | |_| | (_| | |   | || | | (_| | | \__ \
    # |_|   |_|_| |_|\__,_|  \____\___/|_| |_|\___|_|  \___|\__\___| |____/|_|_| |_|  \___|_|  \___|_| |_|\__|_|\__,_|_|   |_||_|  \__,_|_|_|___/
                    
    def compute_differential_effect(self, upper_trail):
        """
        Compute differential effect of the upper trail
        """
        
        params = {"nrounds" : self.RU,
                  "blocksize" : self.blocksize,
                  "mode" : 2,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : 3600,
                  "numberoftrails" : 1,
                  "fixedVariables" : {}}
        for bit in range(self.halfblocksize):
            params["fixedVariables"]["xl_0_" + str(bit)] = upper_trail["xl"][0][bit]
            params["fixedVariables"]["xr_0_" + str(bit)] = upper_trail["xr"][0][bit]
        for bit in range(self.halfblocksize):
            params["fixedVariables"]["xl_" + str(self.RU) + "_" + str(bit)] = upper_trail["xl"][self.RU][bit]
            params["fixedVariables"]["xr_" + str(self.RU) + "_" + str(bit)] = upper_trail["xr"][self.RU][bit]
        UDiff = Diff(params)
        UDiff.make_model()
        output = UDiff.solve()
        return output


#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#  _   _                    ___         _                __                   
# | | | | ___   ___  _ __  |_ _| _ __  | |_  ___  _ __  / _|  __ _   ___  ___ 
# | | | |/ __| / _ \| '__|  | | | '_ \ | __|/ _ \| '__|| |_  / _` | / __|/ _ \
# | |_| |\__ \|  __/| |     | | | | | || |_|  __/| |   |  _|| (_| || (__|  __/
#  \___/ |___/ \___||_|    |___||_| |_| \__|\___||_|   |_|   \__,_| \___|\___|
                                                                            
def loadparameters(args):
    '''
    Extract parameters from the argument list and input file
    '''

    # Load default values
    params = {
            "blocksize": 32,
            "RU": 0,
            "RM": 4,
            "RL": 0,
            "RMU": 0,
            "RML": 0,
            "WU": 20,
            "WM": 4,
            "WL": 20,
            "np" : 8,
            "tl"  : -1,
            "solver"  : "ortools",
            "output"  : "output.tex"}

    # Override parameters if they are set on command line
    if args.blocksize is not None:
        params["blocksize"] = args.blocksize
    if args.RU is not None:
        params["RU"] = args.RU
    if args.RM is not None:
        params["RM"] = args.RM
    if args.RL is not None:
        params["RL"] = args.RL
    if args.RMU is not None:
        params["RMU"] = args.RMU
    if args.RML is not None:
        params["RML"] = args.RML
    if args.WU is not None:
        params["WU"] = args.WU
    if args.WM is not None:
        params["WM"] = args.WM
    if args.WL is not None:
        params["WL"] = args.WL
    if args.np is not None:
        params["np"] = args.np
    if args.timelimit is not None:
        params["timelimit"] = args.timelimit
    if args.solver is not None:
        params["solver"] = args.solver
    if args.output is not None:
        params["output"] = args.output

    return params

def main():
    '''
    Parse the arguments and start the request functionality with the provided
    parameters.
    '''
    
    parser = ArgumentParser(description="This tool finds a nearly optimum differential-linear"
                                        "distinguisher for SIMECK family of block ciphers.",
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-bs", "--blocksize", type=int, default=32, help="Block size in bits")
    parser.add_argument("-RU", type=int, default=2, help="Number of rounds for EU")
    parser.add_argument("-RM", type=int, default=6, help="Number of rounds in the middle")
    parser.add_argument("-RL", type=int, default=2, help="Number of rounds for EL")

    parser.add_argument("-RMU", type=int, default=0, help="Number of rounds passed probabilistically at the beginning of EM")
    parser.add_argument("-RML", type=int, default=0, help="Number of rounds passed probabilistically at the end of EM")
    parser.add_argument("-WU", type=int, default=20, help="Weight factor for EU")
    parser.add_argument("-WM", type=int, default=2, help="Weight factor for EM")
    parser.add_argument("-WL", type=int, default=20, help="Weight factor for EL")

    parser.add_argument("-np", type=int, default=8, help="Number of parallel threads")
    parser.add_argument("-tl", "--timelimit", type=int, default=2.5*36000, help="Time limit in seconds")
    parser.add_argument("-sl", "--solver", default="ortools", type=str,
                        choices=['gecode', 'chuffed', 'coin-bc', 'gurobi', 'picat', 'scip', 'choco', 'ortools', 'cplex', 'cbc'],
                        help="choose a cp solver\n")    
    parser.add_argument("-o", "--output", default="output.tex", type=str, help="Output file name")

    # Parse command line arguments and construct parameter list
    args = parser.parse_args()
    params = loadparameters(args)
    dld = DiffLin(params)
    dld.search()

if __name__ == "__main__":
    main()

# Golden configurations
# bocksize, RU, RM, RL, RMU, RML = 32, 3, 8, 3, 1, 1
