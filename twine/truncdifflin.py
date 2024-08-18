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

Disclaimer: We acknowledge that the TWINE block cipher doesn't adhere to statistical assumptions 
in differential analysis, such as the random sub-key assumption
or Markov cipher assumption. The tool's primary function is to find some bounds
for the security of TWINE against differential and differential-linear cryptanalysis.
"""

from truncdiff import WordTwineDiff
from trunclin import WordTwineLin
import time
from gurobipy import *

class TruncatedDiffLin(WordTwineDiff, WordTwineLin):
    """
    This class is used to find a truncated difflin trail for TWINE block cipher
    """

    count = 0
    def __init__(self, RU, RM, RL, RMU, RML, WU=1, WM=1, WL=1):
        """
        Initialize the main parameters of the difflin trails

        :param RU int: length of EU
        :param RM int: length of EM
        :param RL int: length of EL
        :param WU int: cost of active S-boxes in the upper trail
        :param WL int: cost of active S-boxes in the lower trail
        :param WM int: cost of common active S-boxes between the upper and lower trails
        """

        super().__init__()
        self.lp_file_name = f"warp_{RU}_{RM}_{RL}.lp"
        self.RU = RU
        self.RM = RM
        self.RL = RL
        self.RMU = RMU
        self.RML = RML
        self.WU = WU
        self.WM = WM
        self.WL = WL        

    def constraint_by_xor_pr1(self, a, b, c):
        """
        operation:
        (a, b) |----> c = a + b
        model:
        c - a >= 0
        c - b >= 0
        a + b - c >= 0
        """

        constraints = ""
        constraints += f"{c} - {a} >= 0\n"
        constraints += f"{c} - {b} >= 0\n"
        constraints += f"{a} + {b} - {c} >= 0\n"
        return constraints
    
    def generate_constraint_for_differential_trail(self):
        """
        Generate the constraints describing the upper differential trail
        """

        constraints = ""
        for rn in range(self.RU + self.RM):
            x_in = self.generate_round_x_variables(rn, ul="u")
            x_out = self.generate_round_x_variables(rn + 1, ul="u")
            x_middle = self.inv_permute_nibbles(x_out)
            for nibble in range(8):
                constraints += self.constraints_by_equality(x_in[2*nibble], x_middle[2*nibble])
                if rn < self.RU + self.RMU:
                    constraints += self.constraint_by_trunc_xor(x_in[2*nibble], x_in[2*nibble + 1], x_middle[2*nibble + 1])
                else:
                    constraints += self.constraint_by_xor_pr1(x_in[2*nibble], x_in[2*nibble + 1], x_middle[2*nibble + 1])
        return constraints

    def generate_constraints_for_linear_trail(self):
        """
        Generate the constraints describing the propagation of lower linear trail
        """

        constraints = ""
        for rn in range(self.RM + self.RL):
            x_in = self.generate_round_x_variables(rn, ul="l")
            x_out = self.generate_round_x_variables(rn + 1, ul="l")
            x_middle = self.inv_permute_nibbles(x_out)
            for nibble in range(8):
                constraints += self.constraints_by_equality(x_in[2*nibble + 1], x_middle[2*nibble + 1])
                if rn < self.RM - self.RML:
                    constraints += self.constraint_by_xor_pr1(x_in[2*nibble + 1], x_middle[2*nibble], x_in[2*nibble])
                else:
                    constraints += self.constraint_by_trunc_xor(x_in[2*nibble + 1], x_middle[2*nibble], x_in[2*nibble])
        return constraints

    def generate_linking_vars(self, rn):
        """
        Generate linking variables to model the common active
        S-boxes between upper and lower trails
        """

        s = [f"s_{rn}_{n}" for n in range(8)]
        self.milp_variables.extend(s)
        return s

    def generate_objective_function(self):
        """
        Generate objective function of MILP model
        """

        upper_active_sboxes = []

        for r in range(0, self.RU):
            xu = self.generate_round_x_variables(rn=r, ul="u")
            for i in range(8):
                upper_active_sboxes.append(f"{self.WU} {xu[2*i]}")
        lower_active_sboxes = []
        for r in range(self.RM, self.RM + self.RL):
            xl = self.generate_round_x_variables(rn=r, ul="l")
            for i in range(8):
                lower_active_sboxes.append(f"{self.WL} {xl[2*i + 1]}")
        common_active_sboxes = []
        for r in range(self.RM):
            s = self.generate_linking_vars(r)
            for i in range(8):
                common_active_sboxes.append(f"{self.WM} {s[i]}")
        if upper_active_sboxes == [] and lower_active_sboxes == []:
            objective  = " + ".join(common_active_sboxes)
        elif upper_active_sboxes == []:
            objective  = " + ".join(lower_active_sboxes) + " + " + \
                         " + ".join(common_active_sboxes)
        elif lower_active_sboxes == []:
            objective  = " + ".join(upper_active_sboxes) + " + " + \
                         " + ".join(common_active_sboxes)
        else:
            objective  = " + ".join(upper_active_sboxes) + " + " + \
                         " + ".join(lower_active_sboxes) + " + " + \
                         " + ".join(common_active_sboxes)
        return objective


    def make_model(self):
        """
        Generate the main constrain of our MILP model
        describing the propagation of differential trails in upper and
        lower parts
        """

        constraints = "minimize\n"
        constraints += self.generate_objective_function()
        constraints += "\nsubject to\n"
        constraints += self.generate_constraint_for_differential_trail()
        constraints += self.exclude_trivial_solution(ul="u")
        constraints += self.generate_constraints_for_linear_trail()
        constraints += self.exclude_trivial_solution(ul="l")

        for rn in range(self.RM):
            s = self.generate_linking_vars(rn)
            xu = self.generate_round_x_variables(rn + self.RU, ul="u")
            xl = self.generate_round_x_variables(rn, ul="l")
            for i in range(8):
                constraints += f"{xu[2*i]} - {s[i]} >= 0\n"
                constraints += f"{xl[2*i + 1]} - {s[i]} >= 0\n"
                constraints += f"- {xu[2*i]} - {xl[2*i + 1]} + {s[i]} >= -1\n"
        constraints += self.declare_binary_vars()
        constraints += "end"
        with open(self.lp_file_name, "w") as lpfile:
            lpfile.write(constraints)

    def find_truncated_difflin_trail(self):
        """
        Solve the constructed model minimizing the number of active S-boxes
        """

        self.make_model()
        self.milp_model = read(self.lp_file_name)
        os.remove(self.lp_file_name)
        self.milp_model.setParam(GRB.Param.OutputFlag, True)
        start_time = time.time()

        # self.milp_model.Params.PoolSearchMode = 2
        # # Limit number of solutions
        # self.milp_model.Params.PoolSolutions = 2
        # # Choose solution number 1
        # self.milp_model.Params.SolutionNumber = 1
        ###################
        self.milp_model.optimize()
        ###################
        elapsed_time = time.time() - start_time
        time_line = "Total time to find the trail: %0.02f seconds\n".format(elapsed_time)
        objective_function = self.milp_model.getObjective()
        objective_value = objective_function.getValue()
        print(f"Number of active S-boxes: {objective_value}")

    def parse_solver_output(self):
        '''
        Extract the truncated differential characteristic from the solver output
        '''

        self.upper_trail = dict()
        self.lower_trail = dict()
        self.middle_part = dict()
        get_value_str = lambda t: str(int(self.milp_model.getVarByName(t).Xn))
        get_value_int = lambda t: int(self.milp_model.getVarByName(t).Xn)

        print("\nUpper Truncated Trail:\n")
        for r in range(self.RU + self.RM + 1):
            x_name = self.generate_round_x_variables(rn=r, ul="u")
            x_value = ''.join(list(map(get_value_str, x_name)))
            self.upper_trail[f"x_{r}"] = x_value
            print(x_value)
        print("\n%s\n%s" % ("+"*32, "#"*32))
        print("Lower Truncated Trail:\n")
        for r in range(self.RM + self.RL + 1):
            x_name = self.generate_round_x_variables(rn=r, ul="l")
            x_value = ''.join(list(map(get_value_str, x_name)))
            self.lower_trail[f"x_{r}"] = x_value
            print(x_value)
        print("\n%s\n%s" % ("#"*32, "#"*32))
        print("Middle Part:\n")
        for r in range(self.RM):
            s_name = self.generate_linking_vars(r)
            s_value = '*'.join(list(map(get_value_str, s_name))) + "*"
            self.middle_part[f"s_{r}"] = s_value
            print(s_value)
        s = []
        for r in range(self.RM):
            s.extend(self.generate_linking_vars(r))
        ncs = sum(list(map(get_value_int, s)))
        print(f"\nNumber of common active S-boxes: {ncs}")
        self.middle_part["as"] = ncs
        return self.upper_trail, self.middle_part, self.lower_trail

if __name__ == "__main__":
    RU, RM, RL = 0, 7, 0
    RMU, RML = 0, 0
    WU, WM, WL = 1, 1, 1
    bm = TruncatedDiffLin(RU=RU, RM=RM, RL=RL, RMU=RMU, RML=RML, WU=WU, WM=WM, WL=WL)
    bm.find_truncated_difflin_trail()
    bm.parse_solver_output()