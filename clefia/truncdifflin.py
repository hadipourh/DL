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

from truncdiff import TruncDiffClefia
from trunclin import TruncLinClefia
import time
import uuid
import random
from gurobipy import *

class TruncatedDiffLin(TruncDiffClefia, TruncLinClefia):
    """
    This class is used to find a truncated differential-linear trail for CLEFIA block cipher
    """

    count = 0
    def __init__(self, RU, RL, RM, RMU=0, RML=0, WU=1, WL=1, WM=1):
        """
        Initialize the main parameters of the boomerang trails

        :param RU int: number of rounds covered by only the differential trail
        :param RL int: number of rounds covered by only the linear trail
        :param RM int: number of rounds covered by both the linear and differential trails (middle part)
        :param RMU int: number of rounds passing probabilsitically through EM in the upper part
        :param RML int: number of rounds passing probabilsitically through EM in the lower part
        :param WU int: cost of active S-boxes in the differential trail
        :param WL int: cost of active S-boxes in the linear trail
        :param WM int: cost of common active S-boxes between the differential and linear trails
        """

        super().__init__()                
        self.RU = RU
        self.RM = RM
        self.RL = RL
        self.RMU = RMU
        self.RML = RML
        self.total_nrounds = RU + RM + RL
        self.WU = WU
        self.WM = WM
        self.WL = WL
        self.lp_file_name = f"clefia_{self.total_nrounds}_{RU}_{RM}_{RL}_{uuid.uuid4().hex}.lp"

    def constraint_by_xor(self, a, b, c):
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

    def constraints_by_mds_prob_one(self, dx, dy):
        """
        Generate constraints describing the propagation of truncated
        differential trail through the MDS matrix with probability one

        :param dx list[4]: input of MDS
        :param dy list[4]: output of MDS
        """

        constraints = ""
        for y in dy:
            for x in dx:
                constraints += f"{y} - {x} >= 0\n"
        self.mds_counter += 1
        return constraints

    def generate_constraints_for_differential_trail(self):
        """
        Generate the constraints describing the propagation of
        differential trail throughout EMoEU
        """
        constraints = ""
        for rn in range(self.RU + self.RM):
            x_in = self.generate_round_x_variables(rn, "u")
            x_out = self.generate_round_x_variables(rn + 1, "u")
            x_middle = self.apply_inv_permutation(x_out)
            z = self.generate_round_z_variables(rn, "u")
            if rn < self.RU + self.RMU:
                constraints += self.constraints_by_mds(dx=x_in[0], dy=z[0], pr_one=False)
                constraints += self.constraints_by_mds(dx=x_in[2], dy=z[1], pr_one=False)
            else:                
                constraints += self.constraints_by_mds(dx=x_in[0], dy=z[0], pr_one=True)
                constraints += self.constraints_by_mds(dx=x_in[2], dy=z[1], pr_one=True)
            for n in range(4):
                if rn < self.RU + self.RMU:
                    constraints += self.constraint_by_trunc_xor(z[0][n], x_in[1][n], x_middle[1][n])
                    constraints += self.constraint_by_trunc_xor(z[1][n], x_in[3][n], x_middle[3][n])
                else:
                    constraints += self.constraint_by_xor(z[0][n], x_in[1][n], x_middle[1][n])
                    constraints += self.constraint_by_xor(z[1][n], x_in[3][n], x_middle[3][n])
                constraints += self.constraints_by_equality(x_in[0][n], x_middle[0][n])
                constraints += self.constraints_by_equality(x_in[2][n], x_middle[2][n])
            if rn >= 4:
                constraints += self.diffusion_switching_mechanism_diff(rn, "u")                
        return constraints

    def generate_constraints_for_linear_trail(self):
        """
        Generate the constraints describing the propagation of
        linear trail throughout ELoEM
        """

        constraints = ""
        for rn in range(self.RM + self.RL):
            x_in = self.generate_round_x_variables(rn, "l")
            x_out = self.generate_round_x_variables(rn + 1, "l")
            x_middle = self.apply_inv_permutation(x_out)
            z = self.generate_round_z_variables(rn, "l")
            if rn < self.RM - self.RML:
                constraints += self.constraints_by_mds(dx=x_in[1], dy=z[0], pr_one=True)
                constraints += self.constraints_by_mds(dx=x_in[3], dy=z[1], pr_one=True)
            else:
                constraints += self.constraints_by_mds(dx=x_in[1], dy=z[0], pr_one=False)
                constraints += self.constraints_by_mds(dx=x_in[3], dy=z[1], pr_one=False)
            for n in range(4):
                if rn < self.RM - self.RML:
                    constraints += self.constraint_by_xor(x_middle[0][n], z[0][n], x_in[0][n])
                    constraints += self.constraint_by_xor(x_middle[2][n], z[1][n], x_in[2][n])
                else:
                    constraints += self.constraint_by_trunc_xor(x_middle[0][n], z[0][n], x_in[0][n])
                    constraints += self.constraint_by_trunc_xor(x_middle[2][n], z[1][n], x_in[2][n])
                constraints += self.constraints_by_equality(x_in[1][n], x_middle[1][n])
                constraints += self.constraints_by_equality(x_in[3][n], x_middle[3][n])
            if rn >= 2:
                constraints += self.diffusion_switching_mechanism_lin(rn, "l")                
        return constraints

    def generate_linking_vars(self, rn):
        """
        Generate linking variables to model the common active
        S-boxes between differential and linear trails
        """

        s = [[f"s_{rn}_{bn}_{n}" for n in range(4)] for bn in range(2)]
        self.milp_variables.extend(self.flatten_byte_state(s))
        return s

    def generate_objective_function(self):
        """
        Generate objective function of MILP model
        """

        upper_active_sboxes = []

        for r in range(0, self.RU):
            xu = self.generate_round_x_variables(rn=r, ul="u")
            for i in range(4):
                if i % 2 == 0:
                    upper_active_sboxes.append(f"{self.WU*4.67} {xu[0][i]}")
                    upper_active_sboxes.append(f"{self.WU*6} {xu[2][i]}")
                else:
                    upper_active_sboxes.append(f"{self.WU*6} {xu[0][i]}")
                    upper_active_sboxes.append(f"{self.WU*4.67} {xu[2][i]}")
        lower_active_sboxes = []
        for r in range(self.RM, self.RM + self.RL):
            zl = self.generate_round_z_variables(rn=r, ul="l")
            for i in range(4):
                if i % 2 == 0:
                    lower_active_sboxes.append(f"{self.WL*4.83} {zl[0][i]}")
                    lower_active_sboxes.append(f"{self.WL*6} {zl[1][i]}")
                else:
                    lower_active_sboxes.append(f"{self.WL*6} {zl[0][i]}")
                    lower_active_sboxes.append(f"{self.WL*4.83} {zl[1][i]}")
        common_active_sboxes = []
        for r in range(self.RM):
            s = self.generate_linking_vars(r)
            for i in range(4):
                if i % 2 == 0:
                    common_active_sboxes.append(f"{self.WM*3.41} {s[0][i]}")
                    common_active_sboxes.append(f"{self.WM*4} {s[1][i]}")
                else:
                    common_active_sboxes.append(f"{self.WM*4} {s[0][i]}")
                    common_active_sboxes.append(f"{self.WM*3.41} {s[1][i]}")
        if upper_active_sboxes == [] and lower_active_sboxes == []:
            objective  = " + ".join(common_active_sboxes)
        elif upper_active_sboxes == [] and lower_active_sboxes != [] and common_active_sboxes != []:
            objective  = " + ".join(lower_active_sboxes) + " + " + \
                         " + ".join(common_active_sboxes)
        elif lower_active_sboxes == [] and upper_active_sboxes != [] and common_active_sboxes != []:
            objective  = " + ".join(upper_active_sboxes) + " + " + \
                         " + ".join(common_active_sboxes)
        elif common_active_sboxes == [] and upper_active_sboxes != [] and lower_active_sboxes != []:
            objective  = " + ".join(upper_active_sboxes) + " + " + \
                         " + ".join(lower_active_sboxes)
        else:
            objective  = " + ".join(upper_active_sboxes) + " + " + \
                         " + ".join(lower_active_sboxes) + " + " + \
                         " + ".join(common_active_sboxes)
        return objective


    def make_model(self):
        """
        Generate the main constrain of our MILP model
        describing the propagation of differential trails in differential and
        linear parts
        """

        constraints = "minimize\n"
        constraints += self.generate_objective_function()
        constraints += "\nsubject to\n"
        constraints += self.generate_constraints_for_differential_trail()
        constraints += self.exclude_trivial_solution(ul="u")
        constraints += self.generate_constraints_for_linear_trail()
        constraints += self.exclude_trivial_solution(ul="l")

        for rn in range(self.RM):
            s = self.generate_linking_vars(rn)
            xu = self.generate_round_x_variables(rn + self.RU, ul="u")            
            zl = self.generate_round_z_variables(rn, ul="l")
            for i in range(4):
                constraints += f"{xu[0][i]} - {s[0][i]} >= 0\n"
                constraints += f"{xu[2][i]} - {s[1][i]} >= 0\n"
                constraints += f"{zl[0][i]} - {s[0][i]} >= 0\n"
                constraints += f"{zl[1][i]} - {s[1][i]} >= 0\n"
                constraints += f"- {xu[0][i]} - {zl[0][i]} + {s[0][i]} >= -1\n"
                constraints += f"- {xu[2][i]} - {zl[1][i]} + {s[1][i]} >= -1\n"
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

        # self.milp_model.Params.PoolSearchMode = 2
        # # Limit number of solutions
        # self.milp_model.Params.PoolSolutions = 2
        # # Choose solution number 1
        # self.milp_model.Params.SolutionNumber = 0        
        self.milp_model.Params.Seed = random.randint(0, 100000)
        start_time = time.time()
        ###############################################################
        #  ____          _               __  __             _        _ 
        # / ___|   ___  | |__   __ ___  |  \/  |  ___    __| |  ___ | |
        # \___ \  / _ \ | |\ \ / // _ \ | |\/| | / _ \  / _` | / _ \| |
        #  ___) || (_) || | \ V /|  __/ | |  | || (_) || (_| ||  __/| |
        # |____/  \___/ |_|  \_/  \___| |_|  |_| \___/  \__,_| \___||_|
        #
        self.milp_model.optimize()
        ###############################################################
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
        get_value_str = lambda t: str(int(self.milp_model.getVarByName(t).Xn)).zfill(2)
        get_value_int = lambda t: int(self.milp_model.getVarByName(t).Xn)

        print("\nUpper Truncated Trail:\n")
        for r in range(self.RU + self.RM + 1):
            x_name = self.generate_round_x_variables(rn=r, ul="u")
            x_name = self.flatten_byte_state(x_name)
            x_value = ''.join(list(map(get_value_str, x_name)))
            self.upper_trail[f"x_{r}"] = x_value
            print(x_value)
        print("\n%s\n%s" % ("+"*32, "#"*32))
        print("Lower Truncated Trail:\n")
        for r in range(self.RM + self.RL + 1):
            x_name = self.generate_round_x_variables(rn=r, ul="l")
            x_name = self.flatten_byte_state(x_name)
            x_value = ''.join(list(map(get_value_str, x_name)))
            self.lower_trail[f"x_{r}"] = x_value
            print(x_value)
        print("\n%s\n%s" % ("#"*32, "#"*32))
        print("Middle Part:\n")
        for r in range(self.RM):
            s_name = self.generate_linking_vars(r)
            s_name = self.flatten_byte_state(s_name)
            s_value = ''.join(list(map(get_value_str, s_name[0:4]))) + "*"*8 + \
                      ''.join(list(map(get_value_str, s_name[4::]))) + "*"*8
            self.middle_part[f"s_{r}"] = s_value
            print(s_value)
        s = []
        for r in range(self.RM):
            s.extend(self.flatten_byte_state(self.generate_linking_vars(r)))
        ncs = sum(list(map(get_value_int, s)))
        print(f"\nNumber of common active S-boxes: {ncs}")
        self.middle_part["as"] = ncs
        return self.upper_trail, self.middle_part, self.lower_trail

if __name__ == "__main__":
    RU, RM, RL = 0, 4, 0
    RMU, RML = 0, 0
    WU, WM, WL = 1, 1, 1    
    bm = TruncatedDiffLin(RU=RU, RL=RL, RM=RM, RMU=RMU, RML=RML, WU=WU, WL=WL, WM=WM)
    bm.find_truncated_difflin_trail()
    bm.parse_solver_output()