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
for the security of LBlock against differential, linear and differential-linear cryptanalysis.
"""

from truncdiff import WordLBlockDiff
from trunclin import WordLBlockLin
import time
from gurobipy import *

class TruncatedDL(WordLBlockDiff, WordLBlockLin):
    """
    This class is used to find a truncated boomerang trail for LBlock block cipher
    """

    count = 0
    def __init__(self, RU, RL, RM, RMU, RML, WU=1, WL=1, WM=1):
        """
        Initialize the main parameters of the boomerang trails

        :param RU int: number of rounds covered by only the upper trail
        :param RL int: number of rounds covered by only the lower trail
        :param RM int: number of rounds covered by both the lower and upper trails (middle part)
        :param WU int: cost of active S-boxes in the upper trail
        :param WL int: cost of active S-boxes in the lower trail
        :param WM int: cost of common active S-boxes between the upper and lower trails
        """

        super().__init__()
        self.lp_file_name = f"lblock_{RU}_{RM}_{RL}.lp"
        self.RU = RU
        self.RMU = RMU
        self.R0 = RU + RM
        self.RL = RL
        self.RML = RML
        self.R1 = RL + RM
        self.RM = RM
        self.WU = WU
        self.WL = WL
        self.WM = WM
        self.iterative = False

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

    def generate_upper_constraints(self):
        """
        Generate the constraints describing the propagation of
        upper trail
        """
        constraints = ""
        for rn in range(self.R0):
            x_in = self.generate_round_x_variables(rn, ul="u")
            x_out = self.generate_round_x_variables(rn + 1, ul="u")
            x_middle = self.swap(x_out)
            sbo = self.apply_permutation(x_in[0:8])
            for n in range(8):
                constraints += self.constraints_by_equality(x_in[n], x_middle[n])
                if rn < self.RU + self.RMU:
                    constraints += self.constraint_by_trunc_xor(sbo[n], x_in[8 + (n + 2)%8], x_middle[8 + n])
                else:
                    constraints += self.constraint_by_xor_pr1(sbo[n], x_in[8 + (n + 2)%8], x_middle[8 + n])
                    # constraints += self.constraint_by_trunc_xor(sbo[n], x_in[8 + (n + 2)%8], x_middle[8 + n])
        return constraints

    def generate_lower_constraints(self):
        """
        Generate the constraints describing the propagation of
        lower trail
        """

        constraints = ""
        for rn in range(self.R1):
            x_in = self.generate_round_x_variables(rn, ul="l")
            x_out = self.generate_round_x_variables(rn + 1, ul="l")
            x_middle = self.swap(x_out)
            sbi = self.apply_inv_permutation([x_in[8 + (n + 2)%8] for n in range(8)])
            for n in range(8):
                constraints += self.constraints_by_equality(x_in[8 + (n + 2)%8], x_middle[8 + n])
                if rn < self.RM - self.RML:
                    constraints += self.constraint_by_xor_pr1(x_middle[n], sbi[n], x_in[n])
                    # constraints += self.constraints_by_trunc_fork(x_in[n], sbi[n], x_middle[n])
                else:
                    constraints += self.constraint_by_trunc_xor(x_middle[n], sbi[n], x_in[n])
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
                upper_active_sboxes.append(f"{self.WU} {xu[i]}")

        lower_active_sboxes = []
        for r in range(self.RM, self.R1):
            xl = self.generate_round_x_variables(rn=r, ul="l")
            xl = self.apply_inv_permutation([xl[8 + (n + 2)%8] for n in range(8)])
            for i in range(8):
                lower_active_sboxes.append(f"{self.WL} {xl[i]}")

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
        constraints += self.generate_upper_constraints()
        constraints += self.exclude_trivial_solution(ul="u")
        constraints += self.generate_lower_constraints()
        constraints += self.exclude_trivial_solution(ul="l")

        for rn in range(self.RM):
            s = self.generate_linking_vars(rn)
            xu = self.generate_round_x_variables(rn + self.RU, ul="u")
            xl = self.generate_round_x_variables(rn, ul="l")
            xl = self.apply_inv_permutation([xl[8 + (n + 2)%8] for n in range(8)])
            for i in range(8):
                constraints += f"{xu[i]} - {s[i]} >= 0\n"
                constraints += f"{xl[i]} - {s[i]} >= 0\n"
                constraints += f"- {xu[i]} - {xl[i]} + {s[i]} >= -1\n"
        if self.iterative == True:
            x_in = self.generate_round_x_variables(0, ul="u")
            # x_out = self.generate_round_x_variables(self.R1, ul="l")
            x_out = self.generate_round_x_variables(self.RM, ul="u")
            for i in range(16):
                constraints += f"{x_in[i]} - {x_out[i]} = 0\n"
        constraints += self.declare_binary_vars()
        constraints += "end"
        with open(self.lp_file_name, "w") as lpfile:
            lpfile.write(constraints)

    def find_truncated_dl_trail(self):
        """
        Solve the constructed model minimizing the number of active S-boxes
        """

        self.make_model()
        self.milp_model = read(self.lp_file_name)
        os.remove(self.lp_file_name)
        self.milp_model.setParam(GRB.Param.OutputFlag, True)

        self.milp_model.Params.PoolSearchMode = 2
        # Limit number of solutions
        self.milp_model.Params.PoolSolutions = 10
        # Choose solution number 1
        self.milp_model.Params.SolutionNumber = 0

        start_time = time.time()
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
        Extract the discovered truncated trail from the solver's output
        '''

        self.upper_trail = dict()
        self.lower_trail = dict()
        self.middle_part = dict()
        get_value_str = lambda t: str(round(self.milp_model.getVarByName(t).Xn))
        get_value_int = lambda t: int(round(self.milp_model.getVarByName(t).Xn))

        print("\nUpper Truncated Trail:\n")
        for r in range(self.R0 + 1):
            x_name = self.generate_round_x_variables(rn=r, ul="u")
            x_value = ''.join(list(map(get_value_str, x_name)))
            self.upper_trail[f"x_{r}"] = x_value
            print(x_value)
        print("\n%s\n%s" % ("+"*16, "#"*16))
        print("Lower Truncated Trail:\n")
        for r in range(self.R1 + 1):
            x_name = self.generate_round_x_variables(rn=r, ul="l")
            x_value = ''.join(list(map(get_value_str, x_name)))
            self.lower_trail[f"x_{r}"] = x_value
            print(x_value)
        print("\n%s\n%s" % ("#"*16, "#"*16))
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
    RU, RM, RL = 1, 9, 0
    RMU, RML = 0, 0
    WU, WM, WL = 1, 1, 1
    bm = TruncatedDL(RU=RU, RL=RL, RM=RM, RMU=RMU, RML=RML, WU=WU, WL=WL, WM=WM)
    bm.iterative = False
    bm.find_truncated_dl_trail()
    bm.parse_solver_output()