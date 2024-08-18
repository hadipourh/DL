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
from truncdifflin import TruncatedDiffLin
from diff import Diff
from lin import Lin
from plotdistinguisher import *

def main():

    # 7 rounds
    # RU, RM, RL = 1, 5, 1
    # WU, WM, WL = 3, 2, 3

    # 8 rounds
    # RU, RM, RL = 1, 5, 2
    # WU, WM, WL = 5, 4, 1

    # RU, RM, RL = 2, 5, 1
    # WU, WM, WL = 4, 4, 4

    # 9 rounds
    # RU, RM, RL = 1, 5, 3
    # WU, WM, WL = 4, 4, 3

    # RU, RM, RL = 2, 5, 2
    # WU, WM, WL = 4, 4, 4

    # 10 rounds
    # RU, RM, RL = 3, 4, 3
    # WU, WM, WL = 4, 4, 4

    # RU, RM, RL = 3, 4, 3
    # WU, WM, WL = 5, 7, 5

    parser = ArgumentParser(description="This tool finds the nearly optimum differential-linear distinguisher for CLEFIA\n"
                                         "Example:\n"
                                         "python3 attack.py -RU 2 -RM 4 -RL 1 -WU 3 -WM 2 -WL 3",
                            formatter_class=RawTextHelpFormatter)
                        
    parser.add_argument('-i', '--inputfile', type=str, help="Use an input file in yaml format")
    parser.add_argument('-RU', '--RU', type=int,
                        help="number of rounds covered by EU")
    parser.add_argument('-RM', '--RM', type=int,
                        help="number of rounds covered by Em")
    parser.add_argument('-RL', '--RL', type=int,
                        help="number of rounds covered by EL")
    parser.add_argument('-RMU', '--RMU', type=int,
                        help="number of rounds covered by EMU in the upper trail")
    parser.add_argument('-RML', '--RML', type=int,
                        help="number of rounds covered by EML in the lower trail")
    parser.add_argument('-WU', '--WU', type=int,
                        help="cost of active S-boxes in EU")
    parser.add_argument('-WM', '--WM', type=int,
                        help="cost of active S-boxes in Em")
    parser.add_argument('-WL', '--WL', type=int,
                        help="cost of active S-boxes in EL")
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
    
    tex_content = tex_init()
    ##############################################################################################
    ##############################################################################################
    # Step1- Find a truncated differential-linear trail
    DL = TruncatedDiffLin(RU=RU, RL=RL, RM=RM, RMU=RMU, RML=RML, WU=WU, WL=WL, WM=WM)    
    DL.find_truncated_difflin_trail()
    upper_trail, middle_part, lower_trail = DL.parse_solver_output()
    ##############################################################################################
    ##############################################################################################
    # Step2- Instantiate the upper/lower truncated trails with real differential trails
    diff_upper_trail = None
    diff_effect_upper = 0
    if RU > 0:
        time_limit = 10000
        params = {"nrounds" : DL.RU,
                  "mode" : 0,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : time_limit,
                  "numberoftrails" : 1,
                  "fixedVariables" : {}}
        for r in range(DL.RU + 1):
            for bn in range(4):
                for word in range(4):
                    byten = bn*4 + word
                    if upper_trail[f"x_{r}"][2*byten:2*(byten + 1)] == "00":
                        for bit in range(8):
                            params["fixedVariables"][f"x_{r}_{bn}_{word}_{bit}"] = "0"
        diff = Diff(params)
        diff.make_model()
        diff_upper_trail = diff.solve()
        params["fixedVariables"] = {"x_0": diff_upper_trail["x_0"], f"x_{DL.RU}": diff_upper_trail[f"x_{DL.RU}"]}
        params["mode"] = 2
        diff = Diff(params)
        diff.make_model()
        diff_effect_upper = diff.solve()
    ##############################################################################################
    lin_lower_trail = None
    lin_effect_lower = 0
    if RL > 0:
        time_limit = 10000
        params = {"nrounds" : DL.RL,
                  "mode" : 0,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : time_limit,
                  "numberoftrails" : 1,
                  "fixedVariables" : {}}
        for r in range(DL.RL + 1):
            for bn in range(4):
                for word in range(4):
                    byten = 4*bn + word
                    if lower_trail[f"x_{r + DL.RM}"][2*byten:2*(byten + 1)] == "00":
                        for bit in range(8):
                            params["fixedVariables"][f"x_{r}_{bn}_{word}_{bit}"] = "0"
        lin = Lin(params)
        lin.make_model()
        lin_lower_trail = lin.solve()
        params["fixedVariables"] = {"x_0": lin_lower_trail["x_0"], f"x_{DL.RL}": lin_lower_trail[f"x_{DL.RL}"]}
        params["mode"] = 2
        lin = Lin(params)
        lin.make_model()
        lin_effect_lower = lin.solve()
    ##############################################################################################
    ##############################################################################################
    # print out a summary of result on terminal
    print("#"*27)
    print("Summary of the results:")
    print("Upper trail:")
    if diff_upper_trail != None:
        diff.print_trail(trail=diff_upper_trail)
    print("#"*27)
    mactive_sboxes = middle_part["as"]
    print(f"Sandwich {RM} rounds in the middle with {mactive_sboxes} active S-boxes")
    print("#"*27)
    print("Lower trail:")
    if lin_lower_trail != None:
        lin.print_trail(trail=lin_lower_trail)
    print("-"*27)

    total_weight = 0
    if diff_effect_upper != 0:
        print("differential effect of the upper trail: 2^(%0.02f)" % diff_effect_upper)
        total_weight += diff_effect_upper
    if lin_effect_lower != 0:
        print("linear effect of the lower trail: 2^(%0.02f)" % lin_effect_lower)
        total_weight += lin_effect_lower
    upper_bound =  total_weight + (-4)*mactive_sboxes
    lower_bound = total_weight + (-6)*mactive_sboxes
    print("Total correlation = p*r*q^2 = 2^({:.2f}) x r x 2^({:.2f})".format(diff_effect_upper, lin_effect_lower))
    print("2^({:.2f}) <= Total correlation <= 2^({:.2f})".format(lower_bound, upper_bound))
    print("To compute the accurate value of total correlation, r should be evaluated experimentally or using the DLCT framework")

    ##############################################################################################
    ##############################################################################################
    # plot distinguisher
    # if diff_upper_trail != None:
    #     active_input_bits = diff.flatten_state([[4*i + j for j in range(4)] for i in range(16) if diff_upper_trail["x_0"][i] != "0"])
    #     tex_content += tikz_mark_input_bits(active_input_bits, color="red")
    #     tex_content += tex_diff_trail(trail=diff_upper_trail, markpattern="markupperpath", direction="->")
    # else:
    #     active_input_bits = []
    #     for i in range(16):
    #         if upper_trail[f"x_{0}"][i] != "0":
    #             active_input_bits.extend([j for j in range(4*i, 4*(i + 1))])
    #     tex_content += tikz_mark_input_bits(active_input_bits, color="red")

    # tex_content += tex_middle(upper_trail=upper_trail, midd_trail=middle_part, lower_trail=lower_trail, RU=RU, RM=RM, RL=RL)

    # # tex_content += tex_diff_trail(trail=lin_lower_trail, markpattern="marklowerpath", direction="<-")
    # if lin_lower_trail != None:
    #     tex_content += tex_diff_lower_trail(trail=lin_lower_trail, \
    #                                         upper_crossing_difference=[str(i) for i in range(16) if upper_trail[f"x_{RU + RM}"][i] != "0"],\
    #                                         markpattern="marklowerpath",\
    #                                         direction="<-")
    #     active_output_bits = diff.flatten_state([[4*i + j for j in range(4)] for i in range(16) if lin_lower_trail[f"x_{RL}"][i] != "0"])
    #     tex_content += tikz_mark_output_bits(active_output_bits, color="blue")
    # else:
    #     active_output_bits = []
    #     for i in range(16):
    #         if lower_trail[f"x_{RM + RL}"][i] != "0":
    #             active_output_bits.extend([j for j in range(4*i, 4*(i + 1))])
    #     tex_content += tikz_mark_output_bits(active_output_bits, color="blue")

    # tex_content += tex_fin(RU + RM + RL)
    # with open("bmd.tex", "w") as texfile:
    #     texfile.write(tex_content)

def loadparameters(args):
    """
    Get parameters from the argument list and inputfile.
    """

    # Load default values
    params = {"inputfile": "./input.yaml",
                "RU" : 2,
                "RM" : 5,
                "RL" : 1,
                "RMU": 0,
                "RML": 1,
                "WU" : 1,
                "WM" : 4,
                "WL" : 1,
                "timelimit" : 3200,
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
