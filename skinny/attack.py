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
import copy
from random import randint
from xml.dom.expatbuilder import parseString
logging.basicConfig(filename="minizinc-python.log", level=logging.DEBUG)
import time
import minizinc
import datetime
from argparse import ArgumentParser, RawTextHelpFormatter
from pathlib import Path
from differential import Differential
from linear import Linear
from drawdistinguisher import *
import yaml
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


class DL:
    ID_counter = 0

    def __init__(self, param) -> None:
        DL.ID_counter += 1
        self.id = DL.ID_counter
        self.name = "DL" + str(self.id)
        self.type = "DL"
        self.variant = param["variant"]
        self.cell_size = param["cellsize"]
        self.is_related_tweakey = param["is_related_tweakey"]
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
    #############################################################################################################################################    
    #  ____                      _        __              _____                           _           _   _   _                        ___                             _____          _ _     
    # / ___|  ___  __ _ _ __ ___| |__    / _| ___  _ __  |_   _| __ _   _ _ __   ___ __ _| |_ ___  __| | | | | |_ __  _ __   ___ _ __ / / |    _____      _____ _ __  |_   _| __ __ _(_) |___ 
    # \___ \ / _ \/ _` | '__/ __| '_ \  | |_ / _ \| '__|   | || '__| | | | '_ \ / __/ _` | __/ _ \/ _` | | | | | '_ \| '_ \ / _ \ '__/ /| |   / _ \ \ /\ / / _ \ '__|   | || '__/ _` | | / __|
    #  ___) |  __/ (_| | | | (__| | | | |  _| (_) | |      | || |  | |_| | | | | (_| (_| | ||  __/ (_| | | |_| | |_) | |_) |  __/ | / / | |__| (_) \ V  V /  __/ |      | || | | (_| | | \__ \
    # |____/ \___|\__,_|_|  \___|_| |_| |_|  \___/|_|      |_||_|   \__,_|_| |_|\___\__,_|\__\___|\__,_|  \___/| .__/| .__/ \___|_|/_/  |_____\___/ \_/\_/ \___|_|      |_||_|  \__,_|_|_|___/
    #                                                                                                          |_|   |_|                                                                    
    def search(self):
        """
        Search for a rectangle attack
        """

        if self.time_limit != -1:
            time_limit = datetime.timedelta(seconds=self.time_limit)
        else:
            time_limit = None
    
        start_time = time.time()
        ##########################
        ##########################
        # Step 1: find a truncated differential-linear trail
        print("Searching for a truncated differential-linear trail...")
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
        self.cp_inst["NPT"] = self.variant
        self.cp_inst["is_related_tweakey"] = self.is_related_tweakey
        self.cp_inst["cell_size"] = self.cell_size
        self.result = self.cp_inst.solve(timeout=time_limit, 
                                         processes=self.num_of_threads, 
                                         verbose=False, 
                                         debug_output=Path("./debug_output.txt",
                                         intermediate_solutions=True),
                                         random_seed=randint(0, 100))
                                        #  optimisation_level=2)                                              
        elapsed_time = time.time() - start_time
        print("Time used to find a truncated differential-linear trail: {:0.02f}".format(elapsed_time))
        print(self.result.status)
        if self.result.status == minizinc.Status.OPTIMAL_SOLUTION or self.result.status == minizinc.Status.SATISFIED or \
                            self.result.status == minizinc.Status.ALL_SOLUTIONS or self.result.status.UNKNOWN:
            self.attack_summary = self.find_concrete_distinguisher()
            self.Pm = self.result["Pm"]
            self.attack_summary += f"RU = {self.RU}, RM = {self.RM}, RL = {self.RL}, RMU = {self.RMU}, RML = {self.RML}\n"
            self.attack_summary += f"WU = {self.WU}, WM = {self.WM}, WL = {self.WL}\n"
            self.attack_summary += f"Variant = {self.variant}, Cell size = {self.cell_size}\n"
            self.no_active_sboxes_diff = self.result["P0"]
            self.no_common_active_sboxes = self.result["Pm"]
            self.no_active_sboxes_lin = self.result["P1"]            
            self.attack_summary += f"Number of active S-boxes in differential trail = {self.no_active_sboxes_diff}\n"
            self.attack_summary += f"Number of common active S-boxes                = {self.no_common_active_sboxes}\n"
            self.attack_summary += f"Number of active S-boxes in linear trail       = {self.no_active_sboxes_lin}\n"
            self.attack_summary += "log2(PU)                                       = {:0.2f}\n".format(self.P0)
            self.attack_summary += "log2(CM)                                       ~ -{:0.2f}\n".format(self.Pm/1.2)
            self.attack_summary += "log2(CL^2)                                     = {:0.2f}\n".format(self.P1)
            print(self.attack_summary)                    
            draw = DrawDL(self, output_file_name=self.output_file_name)
            draw.generate_attack_shape()
        elif self.result.status == minizinc.Status.UNSATISFIABLE:
            print("Model is unsatisfiable") 
        else:
            print("Solving process was interrupted")
    
    #############################################################################################################################################
    #############################################################################################################################################
    #############################################################################################################################################
    # _____  _             _    ____                                _          _____             _  _      
    # |  ___|(_) _ __    __| |  / ___| ___   _ __    ___  _ __  ___ | |_  ___  |_   _|_ __  __ _ (_)| | ___ 
    # | |_   | || '_ \  / _` | | |    / _ \ | '_ \  / __|| '__|/ _ \| __|/ _ \   | | | '__|/ _` || || |/ __|
    # |  _|  | || | | || (_| | | |___| (_) || | | || (__ | |  |  __/| |_|  __/   | | | |  | (_| || || |\__ \
    # |_|    |_||_| |_| \__,_|  \____|\___/ |_| |_| \___||_|   \___| \__|\___|   |_| |_|   \__,_||_||_||___/
    # Find concrete trails


    def find_concrete_distinguisher(self):
        """
        Find a concrete distinguisher for the given truncated differential-linear trails discovered by the word-based model
        """        

        params_default = {"rounds" : 0,
                "variant" : self.variant,
                "cellsize" : self.cell_size,
                "skipsb": 0,
                "upperbound1" : None,
                "upperbound2" : None,
                "start_round": 0,
                "end_round": None,
                "mode" : 0,
                "sweight" : 0,
                "endweight" : 384,              
                "timelimit" : 60,
                "fixedVariables" : {}}
        distinguisher_io = dict()
        distinguisher_io[f"tku"] = ""
        distinguisher_io[f"tkl"] = ""            
        print("#"*64)
        print("UPPER TRAIL\n")
        # instantiate the upper trail
        params = copy.deepcopy(params_default)
        params["rounds"] = self.RU + self.RM
        if self.cell_size == 4:
            # params["upperbound1"] = 4*self.result["P0m"]
            params["upperbound2"] = 2.2*self.result["P0"]
        else:
            # params["upperbound1"] = 5*self.result["P0m"]
            params["upperbound2"] = 4*self.result["P0"]
        params["end_round"] = self.RU
        for r in range(self.RU + self.RM + 1):
            for cell in range(16):
                if self.result["DXU"][r][cell] == 0:
                    params["fixedVariables"][f"x_{r}_{cell}"] = "0"
        for cell in range(16):
            if self.result["DXU"][0][cell] == 1:
                params["fixedVariables"][f"x_{0}_{cell}"] = "Y"
        for r in range(self.RU + self.RM):
            for cell in range(8):
                if self.result["DSTKU"][r][cell] == 0:
                    params["fixedVariables"][f"tk_{r}_{cell}"] = "0"
        UDiff = Differential(params, exact=True)
        UDiff.make_model()
        status = UDiff.solve(solution_limit=None, mip_focus=0)
        # compute the differential effect for upper trail
        if status == False:
            raise Exception("Failed to find a concrete upper trail!")
        self.upper_trail = UDiff.parse_solver_output()
        params = copy.deepcopy(params_default)
        params["rounds"] = self.RU
        params["mode"] = 2
        params["fixedVariables"] = {"x_0": self.upper_trail["x_0"], f"x_{self.RU}": self.upper_trail[f"x_{self.RU}"]}
        for z in range(self.variant):
            params["fixedVariables"][f"tk{z+1}_0"] = self.upper_trail[f"tk{z+1}_0"]
        if self.cell_size == 4:
            params["upperbound1"] = 3*self.result["P0"]
        else:
            params["upperbound1"] = 4*self.result["P0"]
        UDiffEffect = Differential(params)
        UDiffEffect.make_model()
        self.P0 = UDiffEffect.solve(log=0)
        # generate round tweakeys
        middle_part_up = dict()
        middle_part_up["tku"] = ""
        for z in range(self.variant):
            # fixed_round_tweakey[f"tk{z+1}_0"] = self.upper_trail[f"tk{z+1}_0"]
            middle_part_up["tku"] += self.upper_trail[f"tk{z+1}_{self.RU}"]
            distinguisher_io[f"tku"] += self.upper_trail[f"tk{z+1}_{0}"]
        # middle_part_up["tku"] = UDiff.generate_tweakey(self.RU + 1, fixed_round_tweakey, self.RU)
        middle_part_up["dxu"] = self.upper_trail[f"x_{self.RU}"]
        distinguisher_io["dxu"] = self.upper_trail[f"x_{0}"] 
        ##########################
        print("#"*64)
        print("LOWER TRAIL\n")
        # instantiate the lower trail
        params = copy.deepcopy(params_default)
        if self.RL == 0:
            params["rounds"] = self.RM + self.RL + 1
        else:
            params["rounds"] = self.RM + self.RL
        if self.cell_size == 4:
            params["upperbound2"] = 2.4*self.result["P1"]
        else:
            params["upperbound2"] = 4*self.result["P1"]
        params["start_round"] = self.RM
        params["end_round"] = self.RM + self.RL
        params["mode"] = 0
        for r in range(self.RM + self.RL + 1):
            for cell in range(16):
                if self.result["DXL"][r][cell] == 0:
                    params["fixedVariables"][f"x_{r}_{cell}"] = "0"
        for cell in range(16):
            if self.result["DXL"][self.RM + self.RL][cell] == 1:
                params["fixedVariables"][f"x_{self.RM + self.RL}_{cell}"] = "Y"
        for r in range(self.RM + self.RL):
            for cell in range(8):
                if self.result["DSTKL"][r][cell] == 0:
                    params["fixedVariables"][f"tk_{r}_{cell}"] = "0"
        LLinear = Linear(params, exact=True)
        LLinear.make_model()
        status = LLinear.solve(solution_limit=None, mip_focus=0)
        # compute the differential effect for the lower trail
        if status == False:
            raise Exception("Failed to find a concrete lower trail!")
        self.lower_trail = LLinear.parse_solver_output()
        params = copy.deepcopy(params_default)
        params["rounds"] = self.RL
        params["mode"] = 2
        params["fixedVariables"] = {"x_0": self.lower_trail[f"x_{self.RM}"], f"x_{self.RL}": self.lower_trail[f"x_{self.RM + self.RL}"]}
        if self.cell_size == 4:           
            params["upperbound1"] = 3*self.result["P1"]
        else:
            params["upperbound1"] = 4*self.result["P1"]
        LLinEffect = Linear(params)
        LLinEffect.make_model()
        self.P1 = LLinEffect.solve(log=0)      
        middle_part_low = dict()       
        middle_part_low["dxl"] = self.lower_trail[f"x_{self.RM}"]
        distinguisher_io["dxl"] = self.lower_trail[f"x_{self.RM + self.RL}"]
        ##########################
        attack_summary = "#"*64 + "\n"
        attack_summary += "Middle part\n"
        attack_summary += f"char dk_str[] = \"{middle_part_up['tku']}\";\n"
        attack_summary += f"char dp_str[]  = \"{middle_part_up['dxu']}\";\n"
        attack_summary += f"char dc_str[]  = \"{middle_part_low['dxl']}\";\n"
        attack_summary += "#"*64 + "\n"
        attack_summary += "Distinguisher\n"
        attack_summary += f"char dk_str[] = \"{distinguisher_io['tku']}\";\n"
        attack_summary += f"char dp_str[]  = \"{distinguisher_io['dxu']}\";\n"
        attack_summary += f"char dc_str[]  = \"{distinguisher_io['dxl']}\";\n"
        attack_summary += "#"*64 + "\n"
        return attack_summary

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
    params = {"variant"  : 2,
              "cellsize" : 4,
              "skipsb" : 0,
              "is_related_tweakey": 0,
              "RU" : 8,
              "RM" : 4,
              "RL" : 7,
              "RMU": 0,
              "RML": 0,
              "WU" : 2,
              "WM" : 1,
              "WL" : 2,
              "np" : 8,
              "t"  : 1800000,
              "solver"  : "gurobi",
              "output"  : "output.tex"}

    # Override parameters if they are set on command line
    if args.variant is not None:
        params["variant"] = args.variant
    if args.cellsize is not None:
        params["cellsize"] = args.cellsize
    if args.skipsb is not None:
        params["skipsb"] = args.skipsb
    if args.is_related_tweakey is not None:
        params["is_related_tweakey"] = args.is_related_tweakey
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
                                        "distinguisher for SKINNY family of block ciphers.",
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-v", "--variant", default=3, type=int, help="SKINNY variant (1, 2, 3, 4)")
    parser.add_argument("-c", "--cellsize", default=4, type=int, help="Cell size (4 or 8)")
    parser.add_argument("-sk", "--skipsb", default=0, type=int, help="Skip the first S-box layer (0 or 1)")
    parser.add_argument("-rt", "--is_related_tweakey", default=0, type=int, help="Use related tweakey (0 or 1)")

    parser.add_argument("-RU", type=int, default=0, help="Number of rounds in the first part of the DL distinguisher")
    parser.add_argument("-RM", type=int, default=7, help="Number of rounds in the middle part of the DL distinguisher")
    parser.add_argument("-RL", type=int, default=0, help="Number of rounds in the second part of the DL distinguisher")
    parser.add_argument("-RMU", type=int, default=0, help="Number of initial rounds in EM encoding probabilistically")
    parser.add_argument("-RML", type=int, default=1, help="Number of final rounds in EM encoding probabilistically")
    
    parser.add_argument("-WU", type=int, default=2, help="Weight of Sboxes through E0")
    parser.add_argument("-WM", type=int, default=1, help="Weight of Sboxes through Em")
    parser.add_argument("-WL", type=int, default=2, help="Weight of Sboxes through E1")
    parser.add_argument("-np", type=int, default=8, help="Number of parallel threads")
    parser.add_argument("-t", "--timelimit", type=int, default=360000, help="Time limit in seconds")
    parser.add_argument("-s", "--solver", default="ortools", type=str,
                        choices=['gecode', 'chuffed', 'coin-bc', 'gurobi', 'picat', 'scip', 'choco', 'ortools', 'cplex', 'cbc'],
                        help="choose a cp solver\n")    
    parser.add_argument("-o", "--output", default="output.tex", type=str, help="Output file name")

    # Parse command line arguments and construct parameter list
    args = parser.parse_args()
    params = loadparameters(args)
    bmd = DL(params)
    bmd.search()

if __name__ == "__main__":
    main()
