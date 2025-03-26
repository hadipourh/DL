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
from diff import Diff
from lin import Lin
import io
from contextlib import redirect_stdout
from draw import DrawDL

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
        print(f"Searching a distinguisher for {self.RD} rounds of AES ...")
        self.cp_model = minizinc.Model()
        self.cp_model.add_file(self.mzn_file_name)
        self.cp_inst = minizinc.Instance(solver=self.cp_solver, model=self.cp_model)
        self.cp_inst["RU"] = self.RU
        self.cp_inst["RM"] = self.RM
        self.cp_inst["RL"] = self.RL
        self.cp_inst["WU"] = self.WU
        self.cp_inst["WM"] = self.WM
        self.cp_inst["WL"] = self.WL
        self.cp_inst["RMU"] = self.RMU
        self.cp_inst["RML"] = self.RML
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
            self.attack_summary, self.upper_trail, self.lower_trail = self.parse_minizinc_solution()
            print(self.attack_summary)
            output_buffer = io.StringIO()
            output_buffer.write("#"*50 + "\n")
            output_buffer.write(f"Upper trail:\n")
            if self.RU + self.RMU >= 1:                                
                self.diff, self.diff_trail, self.diff_effect_upper = self.find_differential_trail()
                with redirect_stdout(output_buffer):
                    self.diff.print_trail(self.diff_trail)
            else:                
                self.diff = None
                self.diff_trail = None
                self.diff_effect_upper = 0
            if self.RML + self.RL >= 1:                
                self.lin, self.lin_trail, self.lin_trail_weight = self.find_linear_trail()                
            else:
                self.lin_trail_weight = 0
            output_buffer.write("#"*50 + "\n")               
            output_buffer.write(f"Sandwich {self.RM} rounds in the middle\n")
            if self.RU + self.RMU >= 1:
                input_diff_middle = "".join([self.diff_trail[f"x_{self.RU}"][(8*row + 2*column):(8*row + 2*column + 2)] for column in range(4) for row in range(4)])
            else:
                input_diff_middle = "".join([str(self.upper_trail["x"][self.RU][row][column]).zfill(2) for column in range(4) for row in range(4)])            
            output_buffer.write(f"char DP_STR[] = \"{input_diff_middle}\";\n")            
            if self.RML + self.RL >= 1:
                output_mask_middle = "".join([self.lin_trail[f"x_{self.RML}"][(8*row + 2*column):(8*row + 2*column + 2)] for column in range(4) for row in range(4)])
            else:
                output_mask_middle = "".join([str(self.lower_trail["x"][self.RM][row][column]).zfill(2) for column in range(4) for row in range(4)])            
            output_buffer.write(f"char LC_STR[] = \"{output_mask_middle}\";\n")            
            output_buffer.write("#"*50 + "\n")                        
            output_buffer.write(f"Lower trail:\n")
            if self.RML + self.RL >= 1:
                with redirect_stdout(output_buffer):
                    self.lin.print_trail(self.lin_trail)                
            output_buffer.write("#"*50 + "\n")            
            output_buffer.write(f"Differential effect of the upper differential trail: 2^({self.diff_effect_upper})\n")            
            output_buffer.write(f"Number of common active S-boxes in the middle      : {self.result['CM']}\n")
            output_buffer.write(f"Squared correaltion of the lower linear trail      : 2^({self.lin_trail_weight})\n")
            self.attack_summary += output_buffer.getvalue()
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
    #  _____  _             _    ____                                _          _____             _  _      
    # |  ___|(_) _ __    __| |  / ___| ___   _ __    ___  _ __  ___ | |_  ___  |_   _|_ __  __ _ (_)| | ___ 
    # | |_   | || '_ \  / _` | | |    / _ \ | '_ \  / __|| '__|/ _ \| __|/ _ \   | | | '__|/ _` || || |/ __|
    # |  _|  | || | | || (_| | | |___| (_) || | | || (__ | |  |  __/| |_|  __/   | | | |  | (_| || || |\__ \
    # |_|    |_||_| |_| \__,_|  \____|\___/ |_| |_| \___||_|   \___| \__|\___|   |_| |_|   \__,_||_||_||___/
    # Find concrete trails

    def find_differential_trail(self):
        """
        Find concrete differential trail
        """
        params = {"nrounds" : self.RU + self.RMU,
                "variant": 1,
                "is_related_key": 0,
                "mode" : 0,
                "startweight" : 0,
                "endweight" : 128,
                "timelimit" : 60,
                "numberoftrails" : 1,
                "fixedVariables" : {}}
        for r in range(self.RU + self.RMU + 1):
            for row in range(4):
                for col in range(4):
                    if self.upper_trail["x"][r][row][col] == 0:
                        params["fixedVariables"][f"x_{r}_{row}_{col}"] = "0"
        diff = Diff(params)
        diff.make_model()
        diff_trail = diff.solve()
        if self.RU > 0:
            params = {"nrounds" : self.RU,
                    "variant": 1,
                    "is_related_key": 0,
                    "mode" : 2,
                    "startweight" : 0,
                    "endweight" : 128,
                    "timelimit" : 60,
                    "numberoftrails" : 1,
                    "fixedVariables" : {}}
            params["fixedVariables"] = {"x_0": diff_trail["x_0"], f"x_{self.RU}": diff_trail[f"x_{self.RU}"]}        
            diff = Diff(params)
            diff.make_model()
            diff_effect_upper = diff.solve()
        else: 
            diff_effect_upper = 0
        return diff, diff_trail, diff_effect_upper
    
    def find_linear_trail(self):
        """
        Find concrete linear trail
        """
        
        params = {"nrounds" : self.RML + self.RL,
                "mode" : 0,
                "startweight" : 0,
                "endweight" : 128,
                "timelimit" : 60,
                "numberoftrails" : 1,
                "fixedVariables" : {}}
        for r in range(self.RML + self.RL + 1):
            for row in range(4):
                for col in range(4):
                    if self.lower_trail["x"][self.RM - self.RML + r][row][col] == 0:
                        params["fixedVariables"][f"x_{r}_{row}_{col}"] = "0"
        lin = Lin(params)
        lin.make_model()
        lin_trail = lin.solve()
        if self.RL > 0:
            params = {"nrounds" : self.RL,
                    "mode" : 2,
                    "startweight" : 0,
                    "endweight" : 128,
                    "timelimit" : 60,
                    "numberoftrails" : 1,
                    "fixedVariables" : {}}
            params["fixedVariables"] = {"x_0": lin_trail[f"x_{self.RML}"], f"x_{self.RL}": lin_trail[f"x_{self.RML + self.RL}"]}        
            lin = Lin(params)
            lin.make_model()
            lin_effect_lower = lin.solve()
        else:
            lin_effect_lower = 0
        return lin, lin_trail, lin_effect_lower

    #############################################################################################################################################
    #############################################################################################################################################
    #  ____                           _    _             ____          _         _    _               
    # |  _ \  __ _  _ __  ___   ___  | |_ | |__    ___  / ___|   ___  | | _   _ | |_ (_)  ___   _ __  
    # | |_) |/ _` || '__|/ __| / _ \ | __|| '_ \  / _ \ \___ \  / _ \ | || | | || __|| | / _ \ | '_ \ 
    # |  __/| (_| || |   \__ \|  __/ | |_ | | | ||  __/  ___) || (_) || || |_| || |_ | || (_) || | | |
    # |_|    \__,_||_|   |___/ \___|  \__||_| |_| \___| |____/  \___/ |_| \__,_| \__||_| \___/ |_| |_|
    # Parse the solution and print the distinguisher's specifications

    def parse_minizinc_solution(self):
        """
        Parse the solution and print the distinguisher's specifications
        """
        
        upper_trail = {"x": [[[0 for _ in range(4)] for _ in range(4)] for _ in range(self.RU + self.RM + 1)],
                       "y": [[[0 for _ in range(4)] for _ in range(4)] for _ in range(self.RU + self.RM)]}
        for r in range(self.RU):
            for row in range(4):
                for col in range(4):
                    upper_trail["x"][r][row][col] = self.result["xu"][r][row][col]
                    upper_trail["y"][r][row][col] = self.result["yu"][r][row][col]
        for r in range(self.RU, self.RU + self.RM + 1):
            for row in range(4):
                for col in range(4):
                    upper_trail["x"][r][row][col] = self.result["xmu"][r - self.RU][row][col]
                    if r < self.RU + self.RM:
                        upper_trail["y"][r][row][col] = self.result["ymu"][r - self.RU][row][col]        
        lower_trail = {"x": [[[0 for _ in range(4)] for _ in range(4)] for _ in range(self.RM + self.RL + 1)],
                       "y": [[[0 for _ in range(4)] for _ in range(4)] for _ in range(self.RM + self.RL)]}
        for r in range(self.RM):
            for row in range(4):
                for col in range(4):
                    lower_trail["x"][r][row][col] = self.result["xml"][r][row][col]
                    lower_trail["y"][r][row][col] = self.result["yml"][r][row][col]
        for r in range(self.RM, self.RM + self.RL + 1):
            for row in range(4):
                for col in range(4):
                    lower_trail["x"][r][row][col] = self.result["xl"][r - self.RM][row][col]
                    if r < self.RM + self.RL:
                        lower_trail["y"][r][row][col] = self.result["yl"][r - self.RM][row][col]     
        
        attack_summary = f"Attack summary:\n"
        attack_summary += f"Setting: RU: {self.RU}, RM: {self.RM}, RL: {self.RL}, RMU: {self.RMU}, RML: {self.RML}, WU: {self.WU}, WM: {self.WM}, WL: {self.WL}\n"
        attack_summary += "#"*50 + "\n"
        attack_summary += f"PU:  {self.result['PU']}\n"
        attack_summary += f"CM:  {self.result['CM']}\n"
        attack_summary += f"Q^2: {self.result['QL']}\n"        
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
    params = {
        "RU": 1,
        "RM": 3,
        "RL": 0,
        "RMU": 0,
        "RML": 0,
        "WU": 1,
        "WM": 1,
        "WL": 1,
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
                                        "distinguisher for AES in the single-key setting.",
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-RU", type=int, default=1, help="Number of rounds for EU")
    parser.add_argument("-RM", type=int, default=3, help="Number of rounds in the middle")
    parser.add_argument("-RL", type=int, default=1, help="Number of rounds for EL")

    parser.add_argument("-RMU", type=int, default=0, help="Number of rounds passed probabilistically at the beginning of EM")
    parser.add_argument("-RML", type=int, default=1, help="Number of rounds passed probabilistically at the end of EM")
    parser.add_argument("-WU", type=int, default=1, help="Weight of active S-boxes in EU")
    parser.add_argument("-WM", type=int, default=1, help="Weight of active S-boxes in the middle")
    parser.add_argument("-WL", type=int, default=1, help="Weight of active S-boxes in EL")
    
    parser.add_argument("-np", type=int, default=8, help="Number of parallel threads")
    parser.add_argument("-tl", "--timelimit", type=int, default=60, help="Time limit in seconds")
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

