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

import time
import os
import uuid
from gurobipy import *

class TruncLinClefia:
    """
    This class can be used to find a truncated linear trail
    with minimum number of active S-boxes for CLEFIA block cipher.

    x_roundNumber_branchNumber_byteNumber

    Variable mapping in truncated differential cryptanalysis:

    ... x_r_0                                   --- x_r_1 ...
    ... |                                       |     |
    ... |--------> | S |---|4*4 MDS|--- z_r_0---+---->+    ...
    ... |                                       |    ...

    Variable mapping in truncated linear cryptanalysis:

    ... x_r_0                                   --- x_r_1 ...
    ... |                                       |     |
    ... |----- z_r_0 ---> | S |---|4*4 MDS|-----+---->+    ...
    ... |                                       |    ...

    """
    count = 0

    def __init__(self, nrounds=1) -> None:
        TruncLinClefia.count += 1
        self.xor_counter = 0
        self.mds_counter = 0
        self.DSM_counter = 0
        self.dummy_var = "d"
        self.nrounds = nrounds
        self.milp_variables = []        
        self.lp_file_name = f"clefia_{nrounds}r_{uuid.uuid4()}.lp"        
        self.permute_branches = [3, 0, 1, 2]

    @staticmethod
    def ordered_set(seq):
        """
        This method eliminates duplicated elements in a given list,
        and returns a list in which each elements appears only once
        """

        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    def apply_permutation(self, state):
        """
        Apply the permute branches on the state array
        """

        output = [0]*4
        for i in range(4):
            output[self.permute_branches[i]] = state[i]
        return output

    def apply_inv_permutation(self, state):
        """
        Apply the inverse of branch permutation on the state array
        """

        output = [0]*4
        for i in range(4):
            output[i] = state[self.permute_branches[i]]
        return output

    @staticmethod
    def flatten_byte_state(s):
        state_bytes = [s[bn][byten] for bn in range(len(s)) for byten in range(len(s[0]))]
        return state_bytes

    def generate_round_x_variables(self, rn, ul="u"):
        """
        Generate the input variables of rn'th round

        :param rn int: round number
        :param ul str: 'u' or 'l' denoting whether it is a variable in upper or lower trail
        """

        x = [[f"x{ul}_{rn}_{bn}_{byten}" for byten in range(4)] for bn in range(4)]
        self.milp_variables.extend(self.flatten_byte_state(x))
        return x

    def generate_round_z_variables(self, rn, ul="u"):
        """
        Generate the variables representing the activeness at the output of MDS matrix

        :param rn int: round number
        :param ul str: 'u' or 'l' denoting whether it is a variable in upper or lower trail
        """

        z = [[f"z{ul}_{rn}_{bn}_{byten}" for byten in range(4)] for bn in range(2)]
        self.milp_variables.extend(self.flatten_byte_state(z))
        return z

    def constraint_by_trunc_xor(self, a, b, c, model=1):
        """
        operation:
        (a, b) |----> c = a + b
        model 1:
        a + b + c >= 2 d
        d >= a
        d >= b
        d >= c
        model 2:
        a + b - c >= 0
        a - b + c >= 0
        - a + b + c >= 0
        """

        constraints = ""
        if model == 1:
            d = f"{self.dummy_var}_{self.xor_counter}"
            self.milp_variables.append(d)
            constraints += f"{a} + {b} + {c} -  2 {d} >= 0\n"
            constraints += f"{d} - {a} >= 0\n"
            constraints += f"{d} - {b} >= 0\n"
            constraints += f"{d} - {c} >= 0\n"
            self.xor_counter += 1
        elif model == 2:
            constraints += f"{a} + {b} - {c} >= 0\n"
            constraints += f"{a} - {b} + {c} >= 0\n"
            constraints += f"- {a} + {b} + {c} >= 0\n"
        return constraints

    def constraints_by_equality(self, a, b):
        """
        a = b
        """

        constraint = f"{a} - {b} = 0\n"
        return constraint

    def constraints_by_mds(self, dx, dy, pr_one=False):
        """
        Generate constraints describing the propagation of truncated
        linear trail through the MDS matrix with branch number 5

        :param dx list[4]: input of MDS
        :param dy list[4]: output of MDS
        """

        iodiffs = dx + dy
        iodiffs_sum = " + ".join(iodiffs)
        constraints = f"{iodiffs_sum} - 5 dm_{self.mds_counter} >= 0\n"

        for x in iodiffs:
            constraints += f"dm_{self.mds_counter} - {x} >= 0\n"
        self.mds_counter += 1

        # Model the propagation of linear trails through the MDS with probability one
        if pr_one:
            for y in dy:
                for x in dx:
                    constraints += f"{y} - {x} >= 0\n"
        return constraints

    def diffusion_switching_mechanism_lin(self, rn, ul="u"):
        """
        Generate some conditions to model the diffusion switching mechanism (DSM)
        in CLEFIA. Switching between two different MDS matrices in the diffusion layer of
        CLEFIA makes it stronger against the linear and linear attacks, since the concatenation of
        these two MDS matrices satisfies an special property which guarantess a certain number of active S-boxes
        over 6 consecutive rounds of CLEFIA.

        Reference:
        - On Feistel Structures Using a Diffusion Switching Mechanism (https://www.iacr.org/archive/fse2006/40470042/40470042.pdf)

        :param rn int: round number
        :rtype: string
        :return: constraints modeling the DSM in rounds rn - 2, ..., rn:
        """

        constraints = ""
        
        z_0 = self.generate_round_z_variables(rn - 2, ul)
        z_1 = self.generate_round_z_variables(rn - 1, ul)
        z_2 = self.generate_round_z_variables(rn, ul)

        temp1 = " + ".join(z_0[0] + z_1[0] + z_2[1])
        constraints += f"{temp1} - 5 dsm_{self.DSM_counter} >= 0\n"
        temp2 = " + ".join(z_0[0] + z_2[1])
        constraints += f"{temp2} - dsm_{self.DSM_counter} >= 0\n"
        constraints += f"{temp2} - 8 dsm_{self.DSM_counter} <= 0\n"
        temp3 = " + ".join(z_1[0])
        # constraints += f"{temp3} - dsm_{self.DSM_counter} >= 0\n"
        constraints += f"{temp3} - 4 dsm_{self.DSM_counter} <= 0\n"
        self.milp_variables.append(f"dsm_{self.DSM_counter}")

        self.DSM_counter += 1

        temp1 = " + ".join(z_0[1] + z_1[1] + z_2[0])
        constraints += f"{temp1} - 5 dsm_{self.DSM_counter} >= 0\n"
        temp2 = " + ".join(z_0[1] + z_2[0])
        constraints += f"{temp2} - dsm_{self.DSM_counter} >= 0\n"
        constraints += f"{temp2} - 8 dsm_{self.DSM_counter} <= 0\n"
        temp3 = " + ".join(z_1[1])        
        # constraints += f"{temp3} - dsm_{self.DSM_counter} >= 0\n"
        constraints += f"{temp3} - 4 dsm_{self.DSM_counter} <= 0\n"
        self.milp_variables.append(f"dsm_{self.DSM_counter}")

        self.DSM_counter += 1

        return constraints

    def generate_constraints(self, ul="u"):
        """
        Generate the constraints of MILP model
        """

        constraints = ""
        for rn in range(self.nrounds):
            x_in = self.generate_round_x_variables(rn, ul)
            x_out = self.generate_round_x_variables(rn + 1, ul)
            x_middle = self.apply_inv_permutation(x_out)
            z = self.generate_round_z_variables(rn)            
            constraints += self.constraints_by_mds(dx=z[0], dy=x_in[1])
            constraints += self.constraints_by_mds(dx=z[1], dy=x_in[3])
            for n in range(4):
                constraints += self.constraint_by_trunc_xor(x_in[0][n], z[0][n], x_middle[0][n])
                constraints += self.constraint_by_trunc_xor(x_in[2][n], z[1][n], x_middle[2][n])
                constraints += self.constraints_by_equality(x_in[1][n], x_middle[1][n])
                constraints += self.constraints_by_equality(x_in[3][n], x_middle[3][n])
            if rn >= 2:
                constraints += self.diffusion_switching_mechanism_lin(rn)
        return constraints

    def declare_binary_vars(self):
        """
        Declare binary variables of MILP model
        """

        self.milp_variables = self.ordered_set(self.milp_variables)
        for n in range(self.xor_counter):
            self.milp_variables.append(f"d_{n}")
        for n in range(self.mds_counter):
            self.milp_variables.append(f"dm_{n}")
        constraints = "Binary\n"
        constraints += "\n".join(self.milp_variables) + "\n"
        return constraints

    def generate_objective_function(self):
        """
        Generate the objective function of MILP model
        """

        sbox_inputs = []
        for r in range(self.nrounds):
            round_input = self.generate_round_z_variables(r)
            sbox_inputs.extend(round_input[0] + round_input[1])
        objective = " + ".join(sbox_inputs) + "\n"
        return objective

    def exclude_trivial_solution(self, ul="u"):
        """
        Exclude all-zero solution from the solutions space
        """
        x_0 = self.generate_round_x_variables(0, ul)
        x_0 = self.flatten_byte_state(x_0)
        constraint = " + ".join(x_0) + " >= 1\n"
        return constraint

    def make_model(self):
        """
        Generate the MILP model describing propagation of a truncated linear
        trail through CLEFIA block cipher
        """

        lp_contents = f"\\ Truncated linear trail for {self.nrounds} rounds of CLEFIA\n"
        lp_contents = "minimize\n"
        lp_contents += self.generate_objective_function()
        lp_contents += "subject to\n"
        lp_contents += self.generate_constraints()
        lp_contents += self.exclude_trivial_solution()
        lp_contents += self.declare_binary_vars()
        lp_contents += "End"
        with open(self.lp_file_name, "w") as lp_file:
            lp_file.write(lp_contents)

    def find_truncated_differential_trail(self):
        """
        Solve the constructed model minimizing the number of active S-boxes
        """

        self.make_model()
        self.milp_model = read(self.lp_file_name)
        os.remove(self.lp_file_name)
        self.milp_model.setParam(GRB.Param.OutputFlag, True)
        start_time = time.time()
        ###################
        self.milp_model.optimize()
        ###################
        elapsed_time = time.time() - start_time
        time_line = "Total time to find the trail: %0.02f seconds\n".format(elapsed_time)
        objective_function = self.milp_model.getObjective()
        objective_value = objective_function.getValue()
        print(f"Number of active S-boxes: {objective_value}")
    
    def print_trail(self):
        """
        Print the discovered truncated trail
        """

        for rn in range(self.nrounds):
            x_variables = self.flatten_byte_state(self.generate_round_x_variables(rn))            
            x_variables_0 = x_variables[:4]
            x_variables_1 = x_variables[4:8]
            x_variables_2 = x_variables[8:12]
            x_variables_3 = x_variables[12:16]            
            x_value = "".join(list(map(lambda t: str(int(self.milp_model.getVarByName(t).Xn)), x_variables_0)))
            x_value += " " + "".join(list(map(lambda t: str(int(self.milp_model.getVarByName(t).Xn)), x_variables_1)))
            x_value += " " + "".join(list(map(lambda t: str(int(self.milp_model.getVarByName(t).Xn)), x_variables_2)))
            x_value += " " + "".join(list(map(lambda t: str(int(self.milp_model.getVarByName(t).Xn)), x_variables_3)))
            z_variables = self.flatten_byte_state(self.generate_round_z_variables(rn))
            z_variables_0 = z_variables[:4]
            z_variables_1 = z_variables[4:8]
            z_value = "".join(list(map(lambda t: str(int(self.milp_model.getVarByName(t).Xn)), z_variables_0)))
            z_value += " "*6
            z_value += "".join(list(map(lambda t: str(int(self.milp_model.getVarByName(t).Xn)), z_variables_1)))
            print(x_value)
            print(z_value)

if __name__ == "__main__":
    nrounds = 21
    clefia_upper = TruncLinClefia(nrounds=nrounds)
    clefia_upper.find_truncated_differential_trail()
    clefia_upper.print_trail()