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
in linear analysis, such as the random sub-key assumption
or Markov cipher assumption. The tool's primary function is to find some bounds
for the security of LBlock against linear and linear cryptanalysis.
"""

"""
This module can find truncated linear trail
with minimum number of active S-boxes for LBlock block cipher.
"""

from statistics import mode
import time
import os
from gurobipy import *

class WordLBlockLin:
    """
    This class can be used to find a truncated linear trail
    with minimum number of active S-boxes for LBlock block cipher.

    x_roundNumber_nibbleNumber_bitNumber
    x_roundNumber_nibbleNumber_0: msb
    x_roundNumber_nibbleNumber_3: lsb
    Variable mapping:

    ... x_r_0                            ---  x_r_1  ...
    ... |                                 |     |
    ... |-----x_r_1---> | S | ---x_r_1----+---->+    ...
    ... |                                       |    ...
    """
    count = 0

    def __init__(self, nrounds=1) -> None:
        WordLBlockLin.count += 1
        self.xor_counter = 0
        self.dummy_var = "d"
        self.nrounds = nrounds
        self.milp_variables = []
        self.lp_file_name = f"lblock_{nrounds}r.lp"
        self.permute_nibbles = [2, 0, 3, 1, 6, 4, 7, 5]

    @staticmethod
    def ordered_set(seq):
        """
        This method eliminates duplicated elements in a given list,
        and returns a list in which each elements appears only once
        """

        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    def swap(self, state):
        temp = [0]*16
        for i in range(8):
            temp[i] = state[8 + i]
            temp[8 + i] = state[i]
        return temp

    def apply_permutation(self, half_state):
        """
        Apply the permute nibble operation on a half state array
        """
        output = [0]*8
        for i in range(8):
            output[self.permute_nibbles[i]] = half_state[i]
        return output

    def apply_inv_permutation(self, half_state):
        """
        Apply the inverse permute nibble operation on a half state array
        """

        output = [0]*8
        for i in range(8):
            output[i] = half_state[self.permute_nibbles[i]]
        return output

    def generate_round_x_variables(self, rn, ul="u"):
        """
        Generate the input variables of rn'th round

        :param rn int: round number
        :param ul str: 'u' or 'l' denoting whether it is a variable in upper or lower trail
        """

        x = [f"x{ul}_{rn}_{nibble}" for nibble in range(16)]
        self.milp_variables.extend(x)
        return x

    def constraints_by_trunc_fork(self, a, b, c, model=2):
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

    def generate_constraints(self, ul="u"):
        """
        Generate the constraints of MILP model
        """

        constraints = ""
        for rn in range(self.nrounds):
            x_in = self.generate_round_x_variables(rn, ul)
            x_out = self.generate_round_x_variables(rn + 1, ul)
            x_middle = self.swap(x_out)
            sbi = self.apply_inv_permutation([x_in[8 + (n + 2)%8] for n in range(8)])            
            for n in range(8):
                constraints += self.constraints_by_trunc_fork(x_in[n], sbi[n], x_middle[n])
                constraints += self.constraints_by_equality(x_in[8 + (n + 2)%8], x_middle[8 + n])
        return constraints

    def declare_binary_vars(self):
        """
        Declare binary variables of MILP model
        """

        self.milp_variables = self.ordered_set(self.milp_variables)
        constraints = "Binary\n"
        constraints += "\n".join(self.milp_variables) + "\n"
        return constraints

    def generate_objective_function(self):
        """
        Generate the objective function of MILP model
        """

        sbox_output_masks = []
        for r in range(self.nrounds):
            round_input = self.generate_round_x_variables(r)
            sbox_output_masks.extend([round_input[8 + i] for i in range(8)])
        objective = " + ".join(sbox_output_masks) + "\n"
        return objective

    def exclude_trivial_solution(self, ul="u"):
        """
        Exclude all-zero solution from the solutions space
        """
        x_0 = self.generate_round_x_variables(0, ul)
        constraint = " + ".join(x_0) + " >= 1\n"
        return constraint

    def make_model(self):
        """
        Generate the MILP model describing propagation of a truncated linear
        trail through LBlock block cipher
        """

        lp_contents = f"\\ Truncated linear trail for {self.nrounds} rounds of LBlock\n"
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
        milp_model = read(self.lp_file_name)
        os.remove(self.lp_file_name)
        milp_model.setParam(GRB.Param.OutputFlag, True)
        start_time = time.time()
        ###################
        milp_model.optimize()
        ###################
        elapsed_time = time.time() - start_time
        time_line = "Total time to find the trail: %0.02f seconds\n".format(elapsed_time)
        objective_function = milp_model.getObjective()
        objective_value = objective_function.getValue()
        print(f"Number of active S-boxes: {objective_value}")

if __name__ == "__main__":
    nrounds = 4
    lblock_upper = WordLBlockLin(nrounds=nrounds)
    lblock_upper.find_truncated_differential_trail()