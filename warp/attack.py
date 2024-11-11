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

Disclaimer: We acknowledge that the WARP block cipher doesn't adhere to statistical assumptions 
in differential analysis, such as the random sub-key assumption
or Markov cipher assumption. The tool's primary function is to find some bounds
for the security of WARP against differential, linear, and differential-linear cryptanalysis.
"""

from argparse import ArgumentParser, RawTextHelpFormatter
from truncdifflin import TruncatedDifflin
from diff import Diff
from lin import Lin
from plotdistinguisher import *
import time

# fixed_golden_value_diff = ["1", "0", "1", "0"]
# fixed_golden_value_linear = ["1", "0", "1", "0"]

fixed_golden_value_diff = ["0", "0", "1", "0"]
fixed_golden_value_linear = ["1", "0", "1", "1"]

def main():

    parser = ArgumentParser(description="This tool finds the nearly optimum differential-linear distinguishers for WARP block cipher.\n"
                                         "Example:\n"
                                         "python3 difflin.py -RU 6 -RM 10 -RL 6 -WU 2 -WM 1 -WL 1",
                            formatter_class=RawTextHelpFormatter)
                        
    parser.add_argument('-i', '--inputfile', type=str, help="Use an input file in yaml format")
    parser.add_argument('-RU', '--RU', type=int,
                        help="number of rounds covered by EU")
    parser.add_argument('-RM', '--RM', type=int,
                        help="number of rounds covered by EM")
    parser.add_argument('-RL', '--RL', type=int,
                        help="number of rounds covered by EL")
    parser.add_argument('-RMU', '--RMU', type=int,
                        help="number of rounds covered by EMU")
    parser.add_argument('-RML', '--RML', type=int, 
                        help="number of rounds covered by EML")
    parser.add_argument('-WU', '--WU', type=int,
                        help="cost of active S-boxes in EU")
    parser.add_argument('-WM', '--WM', type=int,
                        help="cost of active S-boxes in EM")
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
    RMU = params["RMU"]
    RML = params["RML"]
    WU, WM, WL = params["WU"], params["WM"], params["WL"]

    assert(RM > 0)
    tex_content = tex_init()
    start_time = time.time()
    ##############################################################################################
    ##############################################################################################
    # Step1- Find a truncated differential-linear trail
    dl = TruncatedDifflin(RU=RU, RM=RM, RL=RL, RMU=RMU, RML=RML, WU=WU, WM=WM, WL=WL)    
    dl.find_truncated_difflin_trail()
    upper_trail, middle_part, lower_trail = dl.parse_solver_output()
    ##############################################################################################
    ##############################################################################################
    # Step2- Instantiate the upper/lower truncated trails with real differential/linear trails
    diff_upper_trail = None
    diff_effect_upper = 0
    if RU != 0:
        time_limit = params["timelimit"]
        params = {"nrounds" : dl.RU,
                  "mode" : 0,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : time_limit,
                  "numberoftrails" : 1,
                  "fixedVariables" : {}}
        for nibble in range(32):
            if upper_trail[f"x_0"][nibble] == "0":
                for bit in range(4):
                    params["fixedVariables"][f"x_{0}_{nibble}_{bit}"] = "0"
            if upper_trail[f"x_{dl.RU}"][nibble] == "0":
                for bit in range(4):
                    params["fixedVariables"][f"x_{dl.RU}_{nibble}_{bit}"] = "0"
            if upper_trail[f"x_{dl.RU}"][nibble] == "1":
                pass
                for bit in range(4):
                    params["fixedVariables"][f"x_{dl.RU}_{nibble}_{bit}"] = fixed_golden_value_diff[bit]
        diff = Diff(params)
        diff.make_model()
        diff_upper_trail = diff.solve()
        params["fixedVariables"] = {"x_0": diff_upper_trail["x_0"], f"x_{dl.RU}": diff_upper_trail[f"x_{dl.RU}"]}
        params["mode"] = 2
        diff = Diff(params)
        diff.make_model()
        diff_effect_upper = diff.solve()
    ##############################################################################################
    lin_lower_trail = None
    lin_effect_lower = 0
    if RL != 0:
        time_limit = params["timelimit"]
        params = {"nrounds" : dl.RL,
                  "mode" : 0,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : time_limit,
                  "numberoftrails" : 1,
                  "fixedVariables" : {}}        
        for nibble in range(32):
            if lower_trail[f"x_{dl.RM}"][nibble] == "0":
                for bit in range(4):
                    params["fixedVariables"][f"x_{0}_{nibble}_{bit}"] = "0"
            if lower_trail[f"x_{dl.RM + dl.RL}"][nibble] == "0":
                for bit in range(4):
                    params["fixedVariables"][f"x_{dl.RL}_{nibble}_{bit}"] = "0"
            if lower_trail[f"x_{dl.RM}"][nibble] == "1":
                pass
                for bit in range(4):
                    params["fixedVariables"][f"x_{0}_{nibble}_{bit}"] = fixed_golden_value_linear[bit]                    
        lin = Lin(params)
        lin.make_model()
        lin_lower_trail = lin.solve()
        params["fixedVariables"] = {"x_0": lin_lower_trail["x_0"], f"x_{dl.RL}": lin_lower_trail[f"x_{dl.RL}"]}
        params["mode"] = 2
        lin = Lin(params)
        lin.make_model()
        lin_effect_lower = lin.solve()
    ##############################################################################################
    ##############################################################################################
    elapsed_time = time.time() - start_time
    # print out a summary of result in terminal    
    stroutput = "#"*55 + "\n" + "Summary of the results:\n"    
    if diff_upper_trail != None:
        stroutput += "A differential trail for EU:\n"
        stroutput += diff.print_trail(trail=diff_upper_trail)
    stroutput += "-"*55 + "\n"
    mactive_sboxes = middle_part["as"]
    stroutput += f"Sandwich {RM} rounds in the middle with {mactive_sboxes} active S-boxes\n"
    stroutput += "-"*55 + "\n"
    if lin_lower_trail != None:
        stroutput += "A linear trail for EL:\n"
        stroutput += lin.print_trail(trail=lin_lower_trail)
    total_weight = 0
    stroutput += "#"*55 + "\n"
    if diff_effect_upper != 0:
        stroutput += "differential effect of the upper trail: 2^(%0.02f)\n" % diff_effect_upper
        total_weight += diff_effect_upper
    if lin_effect_lower != 0:        
        stroutput += "squared correlation of the lower trail: 2^(%0.02f)\n" % lin_effect_lower
        total_weight += lin_effect_lower
    upper_bound =  total_weight + (-0.5)*mactive_sboxes
    lower_bound = total_weight + (-1.5)*mactive_sboxes
    stroutput += "#"*55 + "\n"
    stroutput += "\nTotal correlation = p*r*q^2 = 2^({:.2f}) x r x 2^({:.2f})".format(diff_effect_upper, lin_effect_lower)
    stroutput += "\n2^({:.2f}) <= Total correlation <= 2^({:.2f})".format(lower_bound, upper_bound)
    stroutput += "\nTo compute the accurate value of total probability, r should be evaluated experimentally or using the DLCT framework\n"
    stroutput += f"\nNumber of attacked rounds: {RU + RM + RL}"
    stroutput += f"\nConfiguration: RU={RU}, RM={RM}, RL={RL}, RMU={RMU}, RML={RML}, WU={WU}, WM={WM}, WL={WL}"
    print(stroutput)
    ##############################################################################################
    ##############################################################################################
    # plot distinguisher
    if diff_upper_trail != None:
        active_input_bits = diff.flatten_state([[4*i + j for j in range(4)] for i in range(32) if diff_upper_trail["x_0"][i] != "0"])
        tex_content += tikz_mark_input_bits(active_input_bits, color="tugred")
        tex_content += tex_diff_trail(trail=diff_upper_trail, markpattern="markupperpath", direction="->")
    else:
        active_input_bits = []
        for i in range(32):
            if upper_trail[f"x_{0}"][i] != "0":
                active_input_bits.extend([j for j in range(4*i, 4*(i + 1))])
        tex_content += tikz_mark_input_bits(active_input_bits, color="tugred")

    tex_content += tex_middle(upper_trail=upper_trail, midd_trail=middle_part, lower_trail=lower_trail, RU=RU, RM=RM, RL=RL)

    # tex_content += tex_diff_trail(trail=lin_lower_trail, markpattern="marklowerpath", direction="<-")
    if lin_lower_trail != None:
        tex_content += tex_lin_lower_trail(trail=lin_lower_trail, \
                                           upper_crossing_difference=[str(i) for i in range(32) if upper_trail[f"x_{RU + RM}"][i] != "0"],\
                                           markpattern="marklowerpath",\
                                           direction="<-")
        active_output_bits = lin.flatten_state([[4*i + j for j in range(4)] for i in range(32) if lin_lower_trail[f"x_{RL}"][i] != "0"])
        tex_content += tikz_mark_output_bits(active_output_bits, color="tugblue")
    else:
        active_output_bits = []
        for i in range(32):
            if lower_trail[f"x_{RM + RL}"][i] != "0":
                active_output_bits.extend([j for j in range(4*i, 4*(i + 1))])
        tex_content += tikz_mark_output_bits(active_output_bits, color="tugblue")
    tex_content += r"""\begin{comment}""" + "\n"
    tex_content += stroutput + "\n"
    tex_content += r"""\end{comment}""" + "\n"
    tex_content += tex_fin(RU + RM + RL)
    with open("output.tex", "w") as texfile:
        texfile.write(tex_content)
    # print the elapsed time
    print("Elapsed time: %0.02f seconds" % elapsed_time)

def loadparameters(args):
    """
    Get parameters from the argument list and inputfile.
    """

    # Load default values
    params = {"inputfile": "./input.yaml",
            "RU" : 7,
            "RM" : 9,
            "RL" : 8,

            "RMU": 0,
            "RML": 0,

            "WU" : 4,
            "WM" : 2,
            "WL" : 4, # 1.2 or 2
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
