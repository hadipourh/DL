#!/usr/bin/env python3

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

Disclaimer: We acknowledge that the LBlock block cipher doesn't adhere to statistical assumptions 
in differential analysis, such as the random sub-key assumption
or Markov cipher assumption. The tool's primary function is to find some bounds
for the security of LBlock against differential linear and differential-linear cryptanalysis.
"""

from argparse import ArgumentParser, RawTextHelpFormatter
from truncdifflin import TruncatedDL
from diff import Diff
from lin import Lin 
# from plotdistinguisher import *

golden_value_diff = ["1", "0", "1", "0"]
# golden_value_diff = ["0", "0", "1", "1"]
golden_value_lin = ["0", "0", "0", "1"]

def main():

    parser = ArgumentParser(description="This tool finds the nearly optimum differential-linear distinguisher\n"
                                         "Example:\n"
                                         "python3 difflin.py -RU 4 -RM 8 -RL 4 -WU 6 -WM 3 -WL 6",
                            formatter_class=RawTextHelpFormatter)
                        
    parser.add_argument('-i', '--inputfile', type=str, help="Use an input file in yaml format")
    parser.add_argument('-RU', '--RU', type=int,
                        help="number of rounds covered by E0")
    parser.add_argument('-RM', '--RM', type=int,
                        help="number of rounds covered by Em")
    parser.add_argument('-RL', '--RL', type=int,
                        help="number of rounds covered by E1")
    parser.add_argument('-RMU', '--RMU', type=int,
                        help="number of rounds covered by E0m")
    parser.add_argument('-RML', '--RML', type=int,
                        help="number of rounds covered by E1m")
    parser.add_argument('-WU', '--WU', type=int,
                        help="cost of active S-boxes in E0")
    parser.add_argument('-WM', '--WM', type=int,
                        help="cost of active S-boxes in Em")
    parser.add_argument('-WL', '--WL', type=int,
                        help="cost of active S-boxes in E1")
    parser.add_argument('-tl', '--timelimit', type=int,
                        help="time limit in seconds")
    parser.add_argument('-ns', '--numofsols', type=int,
                        help="number of solutions (currently disabled)")

    # Parse command line arguments and construct parameter list.
    args = parser.parse_args()
    params = loadparameters(args)
    RU, RM, RL = params["RU"], params["RM"], params["RL"]
    RMU, RML = params["RMU"], params["RML"]
    WU, WM, WL = params["WU"], params["WM"], params["WL"]

    assert(RM > 0)
    # tex_content = tex_init()
    ##############################################################################################
    ##############################################################################################
    # Step1- Find a truncated differential-linear trail
    DL = TruncatedDL(RU=RU, RL=RL, RM=RM, RMU=RMU, RML=RML, WU=WU, WL=WL, WM=WM)
    DL.iterative = False
    DL.find_truncated_dl_trail()
    trunc_upper_trail, middle_part, trunc_lower_trail = DL.parse_solver_output()
    ##############################################################################################
    ##############################################################################################
    # Step2- Instantiate the upper/lower truncated trails with real differential trails
    upper_trail = None
    diff_effect_upper = 0
    if RU != 0:
        time_limit = 18000
        params = {"nrounds" : DL.RU,
                  "mode" : 0,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : time_limit,
                  "numberoftrails" : 1,
                  "fixedVariables" : {}}
        for nibble in range(16):
            if trunc_upper_trail[f"x_0"][nibble] == "0":
                for bit in range(4):
                    params["fixedVariables"][f"x_{0}_{nibble}_{bit}"] = "0"
            if trunc_upper_trail[f"x_{DL.RU}"][nibble] == "0":
                for bit in range(4):
                    params["fixedVariables"][f"x_{DL.RU}_{nibble}_{bit}"] = "0"
            if trunc_upper_trail[f"x_{DL.RU}"][nibble] == "1":
                for bit in range(4):
                    params["fixedVariables"][f"x_{DL.RU}_{nibble}_{bit}"] = golden_value_diff[bit]
        diff = Diff(params)
        diff.make_model()
        upper_trail = diff.solve()
        params["fixedVariables"] = {"x_0": upper_trail["x_0"], f"x_{DL.RU}": upper_trail[f"x_{DL.RU}"]}
        params["mode"] = 2
        diff = Diff(params)
        diff.make_model()
        diff_effect_upper = diff.solve()
    ##############################################################################################
    lower_trail = None
    lin_effect_lower = 0
    if RL != 0:
        time_limit = 18000
        params = {"nrounds" : DL.RL,
                  "mode" : 0,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : time_limit,
                  "numberoftrails" : 1,
                  "fixedVariables" : {}}
        for nibble in range(16):
            if trunc_lower_trail[f"x_{DL.RM}"][nibble] == "0":
                for bit in range(4):
                    params["fixedVariables"][f"x_{0}_{nibble}_{bit}"] = "0"
            if trunc_lower_trail[f"x_{DL.RM}"][nibble] == "1":
                for bit in range(4):
                    params["fixedVariables"][f"x_{0}_{nibble}_{bit}"] = golden_value_lin[bit]
            if trunc_lower_trail[f"x_{DL.R1}"][nibble] == "0":
                for bit in range(4):
                    params["fixedVariables"][f"x_{DL.RL}_{nibble}_{bit}"] = "0"
        lin = Lin(params)
        lin.make_model()
        lower_trail = lin.solve()
        params["fixedVariables"] = {"x_0": lower_trail["x_0"], f"x_{DL.RL}": lower_trail[f"x_{DL.RL}"]}
        params["mode"] = 0
        lin = Lin(params)
        lin.make_model()
        linear_trail = lin.solve()
        lin_effect_lower = -1*float(linear_trail['total_weight'])
    ##############################################################################################
    ##############################################################################################
    # print out a summary of result on terminal
    print("#"*27)
    print("Summary of the results:")
    print("Upper trail:")
    if upper_trail != None:
        diff.print_trail(diff_trail=upper_trail)
    print("#"*27)
    mactive_sboxes = middle_part["as"]
    print(f"Sandwich {RM} rounds in the middle with {mactive_sboxes} active S-boxes")
    print("#"*27)
    print("Lower trail:")
    if lower_trail != None:
        lin.print_trail(diff_trail=lower_trail)
    print("-"*27)
    total_weight = 0
    if diff_effect_upper != 0:
        print("differential effect of the upper trail: 2^(%0.02f)" % diff_effect_upper)
        total_weight += diff_effect_upper
    if lin_effect_lower != 0:
        print("correlation of the lower trail        : 2^(%0.02f)" % lin_effect_lower)
        total_weight += lin_effect_lower
    upper_bound =  total_weight + (-1)*mactive_sboxes
    lower_bound = total_weight + (-2)*mactive_sboxes
    print("Total correlation = p*r*q^2 = 2^({:.2f}) x r x 2^({:.2f})".format(diff_effect_upper, lin_effect_lower))
    print("2^({:.2f}) <= Total correlation <= 2^({:.2f})".format(lower_bound, upper_bound))
    print("To compute the accurate value of total correlation, evaluate 'r' either experimentally or by using the DLCT framework.")

    ##############################################################################################
    ##############################################################################################
    # plot distinguisher
    # To do: add the plotter

def loadparameters(args):
    """
    Get parameters from the argument list and inputfile.
    """

    # Load default values
    params = {"inputfile": "./input.yaml",
                "RU" : 4,
                "RM" : 9,
                "RL" : 4,

                "RMU": 0,
                "RML": 0,

                "WU" : 4,
                "WM" : 2.2,
                "WL" : 4,
                "timelimit" : 1200,
                "numofsols" : 1}

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
    
    if args.inputfile:
        params["inputfile"] = args.inputfile
    
    if args.RU != None:
        params["RU"] = args.RU

    if args.RM != None:
        params["RM"] = args.RM

    if args.RL != None:
        params["RL"] = args.RL
    
    if args.RMU != None:
        params["RMU"] = args.RMU

    if args.RML != None:
        params["RML"] = args.RML        

    if args.WU != None:
        params["WU"] = args.WU

    if args.WM != None:
        params["WM"] = args.WM

    if args.WL != None:
        params["WL"] = args.WL
    
    if args.timelimit != None:
        params["timelimit"] = args.timelimit

    if args.numofsols != None:
        params["numofsols"] = args.numofsols

    return params

if __name__ == "__main__":
    main()
