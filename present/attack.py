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

import logging
from pathlib import Path
from random import randint
logging.basicConfig(filename="minizinc-python.log", level=logging.DEBUG)
import time
import minizinc
import datetime
from argparse import ArgumentParser, RawTextHelpFormatter
from drawdistinguisher import *
from diff import Diff
import subprocess

# Check if "OR Tools" appears in the output of "minizinc --solvers" command 
try:
    output = subprocess.run(['minizinc', '--solvers'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if "cp-sat" in output.stdout.decode("utf-8"):
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
        self.RU = param["RU"]
        self.RM = param["RM"]
        self.RL = param["RL"]
        self.WU = param["WU"]
        self.WM = param["WM"]
        self.WL = param["WL"]
        self.RMU = param["RMU"]
        self.RML = param["RML"]
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
        print(f"Searching a distinguisher for {self.RD} rounds of PRESENT ...")
        self.cp_model = minizinc.Model()
        self.cp_model.add_file(self.mzn_file_name)
        self.cp_inst = minizinc.Instance(solver=self.cp_solver, model=self.cp_model)
        self.cp_inst["RU"] = self.RU
        self.cp_inst["RM"] = self.RM
        self.cp_inst["RL"] = self.RL
        self.cp_inst["RMU"] = self.RMU
        self.cp_inst["RML"] = self.RML
        self.cp_inst["WU"] = self.WU
        self.cp_inst["WM"] = self.WM
        self.cp_inst["WL"] = self.WL
        self.cp_inst["offset"] = 0
        self.result = self.cp_inst.solve(timeout=time_limit, 
                                         processes=self.num_of_threads, 
                                         verbose=False, 
                                         debug_output=Path("./debug_output.txt",
                                         intermediate_solutions=True),
                                         random_seed=randint(0, 100),
                                         optimisation_level=2)
        #############################################################################################################################################
        elapsed_time = time.time() - start_time
        print("Time used to find a distinguisher: {:0.02f} seconds".format(elapsed_time))
        print(f"Solver status: {self.result.status}")
        if minizinc.Status.has_solution(self.result.status) or self.result.status == minizinc.Status.ERROR:
            self.attack_summary, self.upper_trail, self.lower_trail = self.parse_solution()            
            diff_effect = self.compute_diff_effect()
            self.attack_summary = self.attack_summary + f"Diff. effect: 2^({diff_effect})\n"
            print(self.attack_summary)
            draw = DrawDL(self, output_file_name=self.output_file_name)
            draw.generate_distinguisher_shape()
        elif self.result.status == minizinc.Status.UNSATISFIABLE:
            print("Model is unsatisfiable") 
        elif self.result.status == minizinc.Status.UNKNOWN:
            print("Unknown error!")
        else:
            print("Solving process was interrupted")

    #############################################################################################################################################
    #############################################################################################################################################
    #   ____                                 _          ____   _   __   __                          _    _         _   _____   __   __              _   
    #  / ___| ___   _ __ ___   _ __   _   _ | |_  ___  |  _ \ (_) / _| / _|  ___  _ __  ___  _ __  | |_ (_)  __ _ | | | ____| / _| / _|  ___   ___ | |_ 
    # | |    / _ \ | '_ ` _ \ | '_ \ | | | || __|/ _ \ | | | || || |_ | |_  / _ \| '__|/ _ \| '_ \ | __|| | / _` || | |  _|  | |_ | |_  / _ \ / __|| __|
    # | |___| (_) || | | | | || |_) || |_| || |_|  __/ | |_| || ||  _||  _||  __/| |  |  __/| | | || |_ | || (_| || | | |___ |  _||  _||  __/| (__ | |_ 
    #  \____|\___/ |_| |_| |_|| .__/  \__,_| \__|\___| |____/ |_||_|  |_|   \___||_|   \___||_| |_| \__||_| \__,_||_| |_____||_|  |_|   \___| \___| \__|
    #                         |_|                                                                                                                       
    # Compute differential effect
    
    def compute_diff_effect(self):
        """
        Compute the differential effect of the differential characteristic
        considering the clustering effect
        """
        time_limit = 10000
        params = {"nrounds" : self.RU,
                  "mode" : 2,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : time_limit,
                  "numberoftrails" : 1,
                  "fixedVariables" : {}}

        for bit in range(64):
            params["fixedVariables"][f"x_0_{bit}"] = self.upper_trail["x"][0][bit]
        for bit in range(64):
            params["fixedVariables"][f"x_{self.RU}_{bit}"] = self.upper_trail["x"][self.RU][bit]
        diff = Diff(params)
        diff.make_model()
        diff_effect_upper = diff.solve()
        return diff_effect_upper

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
        
        upper_trail = {"x": [[0 for _ in range(64)] for _ in range(self.RU + self.RM + 1)],
                       "y": [[0 for _ in range(64)] for _ in range(self.RU + self.RM)]}
        for r in range(self.RU):           
            upper_trail["x"][r] = self.result["xu"][r]
            upper_trail["y"][r] = self.result["yu"][r]
        for r in range(self.RU, self.RU + self.RM + 1):
            upper_trail["x"][r] = self.result["xmu"][r - self.RU]
            if r < self.RU + self.RM:
                upper_trail["y"][r] = self.result["ymu"][r - self.RU]
        lower_trail = {"x": [[0 for _ in range(64)] for _ in range(self.RM + self.RL + 1)],
                       "y": [[0 for _ in range(64)] for _ in range(self.RM + self.RL)]}
        for r in range(self.RM):
            lower_trail["x"][r] = self.result["xml"][r]
            lower_trail["y"][r] = self.result["yml"][r]
        for r in range(self.RM, self.RM + self.RL + 1):
            lower_trail["x"][r] = self.result["xl"][r - self.RM]
            if r < self.RM + self.RL:
                lower_trail["y"][r] = self.result["yl"][r - self.RM]
        input_diff = f"uint64_t inputdiff = 0x" + hex(int("".join(list(map(str, self.result["xu"][0]))), 2))[2:].zfill(16) + ";\n"        
        input_diff_middle = f"uint64_t inputdiff = 0x" + hex(int("".join(list(map(str, self.result["xmu"][0]))), 2))[2:].zfill(16) + ";\n"                
        output_mask_middle = f"uint64_t outputmask = 0x" + hex(int("".join(list(map(str, self.result["xml"][self.RM]))), 2))[2:].zfill(16) + ";\n"        
        output_mask = f"uint64_t outputmask = 0x" + hex(int("".join(list(map(str, self.result["xl"][self.RL]))), 2))[2:].zfill(16) + ";\n"
        
        attack_summary = f"Attack summary:\n"
        attack_summary += f"Setting: RU: {self.RU}, RM: {self.RM}, RL: {self.RL}, RMU: {self.RMU}, RML: {self.RML}, WU: {self.WU}, WM: {self.WM}, WL: {self.WL}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff.: \n{input_diff}"
        for r in range(self.RU + 1):
            attack_summary += hex(int("".join(list(map(str, upper_trail["x"][r]))), 2))[2:].zfill(16) + "\n"        
        attack_summary += "#"*50 + "\n"
        attack_summary += f"input diff. middle: \n{input_diff_middle}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask middle: \n{output_mask_middle}"
        for r in range(self.RM, self.RL + self.RM + 1):
            attack_summary += hex(int("".join(list(map(str, lower_trail["x"][r]))), 2))[2:].zfill(16) + "\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"output mask: \n{output_mask}"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"PU:  {self.result['PU']}\n"
        attack_summary += f"CMB:  {self.result['CMB']}\n"
        attack_summary += f"CMW:  {self.result['CMW']}\n"
        attack_summary += f"Q^2: {self.result['QL']}\n"
        attack_summary += f"Number of effective S-boxes in the middle:       {self.result['CMW']}\n"
        attack_summary += f"Number of effective bit-positions in the middle: {self.result['CMB']}\n"
        attack_summary += "#"*50 + "\n"
        return attack_summary, upper_trail, lower_trail

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
    params = {"RU": 1,
            "RM": 2,
            "RL": 1,
            "WU": 1,
            "WM": 1,
            "WL": 1,
            "RMU": 0,
            "RML": 0,
            "np" : 8,
            "tl"  : -1,
            "solver"  : "ortools",
            "output"  : "output.tex"}

    # Override parameters if they are set on command line
    if args.RU is not None:
        params["RU"] = args.RU
    if args.RM is not None:
        params["RM"] = args.RM
    if args.RL is not None:
        params["RL"] = args.RL
    if args.WU is not None:
        params["WU"] = args.WU
    if args.WM is not None:
        params["WM"] = args.WM
    if args.WL is not None:
        params["WL"] = args.WL
    if args.RMU is not None:
        params["RMU"] = args.RMU
    if args.RML is not None:
        params["RML"] = args.RML
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
                                        "distinguisher for PRESENT block ciphers.",
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-RU", type=int, default=1, help="Number of rounds for EU")
    parser.add_argument("-RM", type=int, default=9, help="Number of rounds in the middle")
    parser.add_argument("-RL", type=int, default=1, help="Number of rounds for EL")
    
    parser.add_argument("-WU", type=int, default=1, help="Weight of EU")
    parser.add_argument("-WM", type=int, default=2, help="Weight of middle")
    parser.add_argument("-WL", type=int, default=1, help="Weight of EL")
    parser.add_argument("-RMU", type=int, default=0, help="Number of rounds passed probabilistically at the beginning of EM")
    parser.add_argument("-RML", type=int, default=0, help="Number of rounds passed probabilistically at the end of EM")
    parser.add_argument("-np", type=int, default=8, help="Number of parallel threads")
    parser.add_argument("-tl", "--timelimit", type=int, default=1000, help="Time limit in seconds")
    # Fetch available solvers from MiniZinc
    available_solvers = [solver_name for solver_name in minizinc.default_driver.available_solvers().keys()]
    parser.add_argument("-sl", "--solver", default="cp-sat", type=str,
                        choices=available_solvers,
                        help="Choose a CP solver")     
    parser.add_argument("-o", "--output", default="output.tex", type=str, help="Output file name")

    # Parse command line arguments and construct parameter list
    args = parser.parse_args()
    params = loadparameters(args)
    dld = DiffLin(params)
    dld.search()

if __name__ == "__main__":
    main()
