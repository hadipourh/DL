#!/usr/bin/env python3

"""
MIT License

Copyright (c) 2024 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to lo so, subject to the following conditions:

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

from argparse import ArgumentParser, RawTextHelpFormatter
from numpy import diff
import yaml
import time
from gurobipy import *
import math
import os

class Lin:
    """
    This class is used to find linear trail as well as
    computing the linear effect of SIMECK block cipher.
    """

    diff_count = 0

    def __init__(self, params) -> None:
        Lin.diff_count += 1

        self.nrounds = params["nrounds"]
        self.block_size = params["blocksize"]
        self.half_block_size = self.block_size // 2
        self.time_limit = params["timelimit"]
        self.start_weight = params["startweight"]
        self.end_weight = params["endweight"]
        self.fixed_variables = params['fixedVariables']
        self.mode = params['mode']
        self.number_of_trails = params["numberoftrails"]
        self.eps = 1e-3
        
        # SIMECK:
        # self.left_rotation_a0 = 2
        # self.left_rotation_a1 = 8
        # self.left_rotation_a2 = 1

        # SIMECK:
        self.left_rotation_a0 = 1
        self.left_rotation_a1 = 5
        self.left_rotation_a2 = 0
        self.lp_file_name = f"simeck_nr_{self.nrounds}.lp"
        self.result_file_name = f"simeck_nr_{self.nrounds}.txt"

        self.milp_variables = []



    @staticmethod
    def ordered_set(seq):
        """
        This method eliminates duplicated elements in a given list,
        and returns a list in which each elements appears only once
        """

        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    @staticmethod
    def flatten_state(state):
        """
        Flatten a 2D list into a 1D list
        """

        return [item for sublist in state for item in sublist]

    def convert_str_to_binarystatevector(self, str_hex):
        state = list(map(int, list(bin(int(str_hex, base=16))[2:].zfill(self.half_block_size))))        
        return state

    def generate_round_half_x_variables(self, prefix, rn):
        """
        Generate the left or right input variables of the rn'th round
        
        :param int rn: round number
        :param int lr: left or right half of the state
        """

        x = [f"{prefix}_{rn}_{bit}" for bit in range(self.half_block_size)]
        self.milp_variables.extend(x)
        return x

    def generate_round_dummy_variables(self, prefix, rn):
        """
        Generate the variables corresponding to the
        output of S-boxes in rn'th round
        :param int rn: round number
        """

        y = [f"{prefix}_{rn}_{bit}" for bit in range(self.half_block_size)]
        self.milp_variables.extend(y)
        return y

    def generate_round_pr_variables(self, rn):
        """
        Generate the variables encoding the probability of S-boxes
        """

        pr = [[f"pr_{rn}_{bit_position}"] for bit_position in range(self.half_block_size)]
        self.milp_variables.extend(self.flatten_state(pr))
        return pr

    def constraints_by_equality(self, a, b):
        """
        Generate constraints for equality
        a = b
        """
        constraint = f"{a} - {b} = 0\n"
        return constraint

    def xor(self, a, b, c):
        '''
        Generate the constraints of a binary XOR
        a xor b = c can be modeled with 4 inequalities (without definition of dummy variable) by removing all impossible vectors (a, b, c)
        '''

        lp_contents = ""
        lp_contents += f"{a} + {b} - {c} >= 0\n"                
        lp_contents += f"{a} - {b} + {c} >= 0\n"
        lp_contents += f"-1 {a} + {b} + {c} >= 0\n"        
        lp_contents += f"-1 {a} - {b} - {c} >= -2\n"        
        return lp_contents

    def xor3(self, b, a2, a1, a0):
        '''
        Generate the constraints of a three-input XOR  (b = a0 xor a1 xor a2)    
        b - a2 - a1 - a0 >= -2
        - b + a2 - a1 - a0 >= -2
        - b - a2 + a1 - a0 >= -2
        b + a2 + a1 - a0 >= 0
        - b - a2 - a1 + a0 >= -2
        b + a2 - a1 + a0 >= 0
        b - a2 + a1 + a0 >= 0
        - b + a2 + a1 + a0 >= 0
        The above inequalities are derived with QuineMcCluskey algorithm
        '''

        lp_contents = ""
        lp_contents += f"{b} - {a2} - {a1} - {a0} >= -2\n"
        lp_contents += f"-1 {b} + {a2} - {a1} - {a0} >= -2\n"
        lp_contents += f"-1 {b} - {a2} + {a1} - {a0} >= -2\n"
        lp_contents += f"{b} + {a2} + {a1} - {a0} >= 0\n"
        lp_contents += f"-1 {b} - {a2} - {a1} + {a0} >= -2\n"
        lp_contents += f"{b} + {a2} - {a1} + {a0} >= 0\n"
        lp_contents += f"{b} - {a2} + {a1} + {a0} >= 0\n"
        lp_contents += f"-1 {b} + {a2} + {a1} + {a0} >= 0\n"
        return lp_contents
    
    def generate_constraints_for_fork_4(self, a0, a1, a2, a3, b):
        """
        Generate the constraints of a Fork with 4 outputs
        The conditions are the same as the linear constraints for XOR with 4 inputs
        sage: sb = [0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0]
        sage: sa = SboxAnalyzer(sb)
        sage: ddt
        [16  0]
        [ 0 16]
        [ 0 16]
        [16  0]
        [ 0 16]
        [16  0]
        [16  0]
        [ 0 16]
        [ 0 16]
        [16  0]
        [16  0]
        [ 0 16]
        [16  0]
        [ 0 16]
        [ 0 16]
        [16  0]
        sage: cnf, milp = sa.minimized_diff_constraints()
        Simplifying the MILP/SAT constraints ...
        Time used to simplify the constraints: 0.00 seconds
        Number of constraints: 16
        Input:	a0||a1||a2||a3; a0: msb
        Output:	b0; b0: msb
        Weight: 0

        ['- a0 + a1 + a2 + a3 + b0 >= 0',
        'a0 - a1 + a2 + a3 + b0 >= 0',
        'a0 + a1 - a2 + a3 + b0 >= 0',
        '- a0 - a1 - a2 + a3 + b0 >= -2',
        'a0 + a1 + a2 - a3 + b0 >= 0',
        '- a0 - a1 + a2 - a3 + b0 >= -2',
        '- a0 + a1 - a2 - a3 + b0 >= -2',
        'a0 - a1 - a2 - a3 + b0 >= -2',
        'a0 + a1 + a2 + a3 - b0 >= 0',
        '- a0 - a1 + a2 + a3 - b0 >= -2',
        '- a0 + a1 - a2 + a3 - b0 >= -2',
        'a0 - a1 - a2 + a3 - b0 >= -2',
        '- a0 + a1 + a2 - a3 - b0 >= -2',
        'a0 - a1 + a2 - a3 - b0 >= -2',
        'a0 + a1 - a2 - a3 - b0 >= -2',
        '- a0 - a1 - a2 - a3 - b0 >= -4']
        """

        lp_contents = ""
        lp_contents += f"- {a0} + {a1} + {a2} + {a3} + {b} >= 0\n"
        lp_contents += f"{a0} - {a1} + {a2} + {a3} + {b} >= 0\n"
        lp_contents += f"{a0} + {a1} - {a2} + {a3} + {b} >= 0\n"
        lp_contents += f"- {a0} - {a1} - {a2} + {a3} + {b} >= -2\n"
        lp_contents += f"{a0} + {a1} + {a2} - {a3} + {b} >= 0\n"
        lp_contents += f"- {a0} - {a1} + {a2} - {a3} + {b} >= -2\n"
        lp_contents += f"- {a0} + {a1} - {a2} - {a3} + {b} >= -2\n"
        lp_contents += f"{a0} - {a1} - {a2} - {a3} + {b} >= -2\n"
        lp_contents += f"{a0} + {a1} + {a2} + {a3} - {b} >= 0\n"
        lp_contents += f"- {a0} - {a1} + {a2} + {a3} - {b} >= -2\n"
        lp_contents += f"- {a0} + {a1} - {a2} + {a3} - {b} >= -2\n"
        lp_contents += f"{a0} - {a1} - {a2} + {a3} - {b} >= -2\n"
        lp_contents += f"- {a0} + {a1} + {a2} - {a3} - {b} >= -2\n"
        lp_contents += f"{a0} - {a1} + {a2} - {a3} - {b} >= -2\n"
        lp_contents += f"{a0} + {a1} - {a2} - {a3} - {b} >= -2\n"
        lp_contents += f"- {a0} - {a1} - {a2} - {a3} - {b} >= -4\n"
        return lp_contents


    def constraints_by_and(self, li, lo, pr):
        """
        Generate constraints modeling the DDT of S-box

        :param str[4] li: input difference
        :param str[4] lo: output difference
        :param str[3] pr: probability of (li --> lo) such that
                          hamming_weight(pr) = -log2(pr(li --> lo))
        :return constraints encoding the LAT squared of S-box (here it is a simple AND operation):
        :rtype str:
        sage: sb = [0, 0, 0, 1]
        sage: sa = SboxAnalyzer(sb)
        sage: lat = sa.linear_approximation_table(scale='correlation')
        sage: lat
        [   1  1/2]
        [   0  1/2]
        [   0  1/2]
        [   0 -1/2]
        sage: cnf, milp = sa.minimized_linear_constraints()
        Simplifying the MILP/SAT constraints ...
        Time used to simplify the constraints: 0.00 seconds
        Number of constraints: 4
        Input:	a0||a1; a0: msb
        Output:	b0; b0: msb
        Weight: 2.0000 p0
        sage: milp
        ['- a0 + p0 >= 0', '- a1 + p0 >= 0', '- b0 + p0 >= 0', 'b0 - p0 >= 0']

        """
        self.and_inequalities = ['- a0 + p0 >= 0', '- a1 + p0 >= 0', '- b0 + p0 >= 0', 'b0 - p0 >= 0']

        constraints = ""
        for ineq in self.and_inequalities:
            temp = ineq
            for i in range(2):
                temp = temp.replace(f"a{i}", li[i])
            for i in range(1):
                temp = temp.replace(f"b{i}", lo[i])
            for i in range(1):
                temp = temp.replace(f"p{i}", pr[i])
            constraints += temp + "\n"
        return constraints
    
    def generate_objective_function(self):
        """
        Generate the objective function of MILP model
        The objective is minimizing the summation of variables
        which encode the weight (or probability exponent) the
        linear trail
        """

        objective_function = "minimize\n"
        weight = []
        for r in range(self.nrounds):
            weight += [f"2 pr_{r}_{bit_position}" for bit_position in range(self.half_block_size)]
        weight = " + ".join(weight)
        objective_function += weight + "\n"
        return objective_function

    def generate_constraints(self):
        """
        Generate the constraints describing the propagation
        of linear trails through a reduced-round SIMECK
        """

        constraints = "subject to\n"
        for rn in range(self.nrounds):
            xl_in = self.generate_round_half_x_variables(prefix='xl', rn=rn)
            xr_in = self.generate_round_half_x_variables(prefix='xr', rn=rn)
            pr = self.generate_round_pr_variables(rn)
            m0 = self.generate_round_half_x_variables(prefix='m0', rn=rn)
            m1 = self.generate_round_half_x_variables(prefix='m1', rn=rn)
            m2 = self.generate_round_half_x_variables(prefix='m2', rn=rn)
            xl_out = self.generate_round_half_x_variables(prefix='xl', rn=rn + 1)
            xr_out = self.generate_round_half_x_variables(prefix='xr', rn=rn + 1)
            for bit_position in range(self.half_block_size):
                a0 = m0[(bit_position + self.left_rotation_a0)%self.half_block_size]
                a1 = m1[(bit_position + self.left_rotation_a1)%self.half_block_size]
                a2 = m2[(bit_position + self.left_rotation_a2)%self.half_block_size]
                b0 = xr_in[bit_position]
                constraints += self.constraints_by_equality(xr_in[bit_position], a0)
                constraints += self.generate_constraints_for_fork_4(m0[bit_position], m1[bit_position], m2[bit_position], xr_out[bit_position], xl_in[bit_position])
                constraints += self.constraints_by_equality(xr_in[bit_position], xl_out[bit_position])
                constraints += self.constraints_by_and([a1, a2], [b0], pr[bit_position])
        return constraints

    def declare_binary_vars(self):
        """
        Declare binary variables of MILP model
        """

        self.milp_variables = self.ordered_set(self.milp_variables)
        constraints = "Binary\n"
        constraints += "\n".join(self.milp_variables) + "\n"
        return constraints

    def exclude_trivial_trail(self):
        """
        Exclude all-zero solution from the solution space
        """

        input_diff = self.generate_round_half_x_variables(prefix='xl', rn=0) +\
              self.generate_round_half_x_variables(prefix='xr', rn=0)
        input_diff = " + ".join(input_diff)
        constraint = f"{input_diff} >= 1\n"
        return constraint

    def declare_fixed_variables(self):
        lp_contents = ""
        for cond in self.fixed_variables.items():
            var = cond[0]
            val = cond[1]
            var = var.split('_')
            if var[0] not in ['xl', 'xr']:
                print("Error: The variable name should start with xl or xr!")
                return
            if len(var) == 2:
                assert(var[0] in ["xl", "xr"])
                state_vars = self.generate_round_half_x_variables(prefix=var[0], rn=var[1])               
                state_values = list(bin(int(val, base=16))[2:].zfill(self.half_block_size))
                for i in range(self.half_block_size):
                    lp_contents += f"{state_vars[i]} = {state_values[i]}\n"                  
            elif len(var) == 3:
                assert(var[0] in ["xl", "xr"])
                state_vars = f"x_{var[1]}_{var[2]}"
                lp_contents += f"{state_vars} = {val}\n"
            else:
                pass
        return lp_contents

    def make_model(self):
        """
        Build the MILP model to find the best linear trail
        """

        lp_contents = "\\ Linear attack on {} rounds of SIMECK\n".format(self.nrounds)
        lp_contents += self.generate_objective_function()
        lp_contents += self.generate_constraints()
        lp_contents += self.exclude_trivial_trail()
        lp_contents += self.declare_fixed_variables()
        lp_contents += self.declare_binary_vars()
        lp_contents += "end"
        with open(self.lp_file_name, "w") as lp_file:
            lp_file.write(lp_contents)

    def exclude_the_previous_sol(self):
        '''
        Let x{S} be the binary variables. Suppose you have a binary
        solution x* in available from the most recent optimization.
        Let N be the subset of S such that x*[n] = 1 for all n in N
        Then, add the following constraint:
        sum{n in N} x[n] - sum{s in S-N} x[s] <= |N|-1
        '''

        all_vars = self.milp_model.getVars()
        nonzero_vars = [v for v in all_vars if v.x == 1]
        zero_vars = [v for v in all_vars if v.x == 0]
        support = len(nonzero_vars)
        first_term = sum(nonzero_vars)
        second_term = sum(zero_vars)
        lhs = first_term - second_term
        self.milp_model.addConstr(lhs <= support - 1)

    def solve(self):
        output = None
        self.milp_model = read(self.lp_file_name)
        os.remove(self.lp_file_name)
        if self.mode == 0:
            output = self.find_characteristic()
        elif self.mode == 1:
            self.find_multiple_characteristics(self.number_of_trails)
        elif self.mode == 2:
            output = self.compute_linear_effect()
            # self.compute_differential_effect_classic_method()
        else:
            print('Enter a number in [0, 1, 2], for the mode parameter please!')
        return output

    def parse_solver_output(self):
        """
        Extract the linear characteristic from the solver output
        """

        get_bit_value = lambda t: str(int(self.milp_model.getVarByName(t).Xn))
        characteristic = dict()
        for r in range(self.nrounds + 1):
            x = self.generate_round_half_x_variables(prefix='xl', rn=r)
            x_value = hex(int("0b" + "".join(list(map(lambda t: str(int(self.milp_model.getVarByName(t).Xn)), x))), 2))[2:].zfill(self.half_block_size//4)
            characteristic[f"xl_{r}"] = x_value
            x = self.generate_round_half_x_variables(prefix='xr', rn=r)
            x_value = hex(int("0b" + "".join(list(map(lambda t: str(int(self.milp_model.getVarByName(t).Xn)), x))), 2))[2:].zfill(self.half_block_size//4)
            characteristic[f"xr_{r}"] = x_value
        for r in range(self.nrounds):
            round_probability = 2*sum([int(self.milp_model.getVarByName(f"pr_{r}_{bit_position}").Xn) for bit_position in range(self.half_block_size)])
            characteristic[f"pr_{r}"] = f"-{round_probability}"
        characteristic["total_weight"] = "%0.02f" % self.total_weight
        characteristic["nrounds"] = self.nrounds
        return characteristic

    @staticmethod
    def print_trail(trail):
        """
        Print out the discovered linear characteristic
        """

        header = ['xl', 'xr', 'cr']
        # Print everthing
        trail_values = map(str, trail.values())
        col_width = max(len(s) for s in trail_values) + 2
        header_str = "Rounds\t"
        data_str = ""
        current_row = 0
        for entry in header[0:-2]:
            header_str += entry.ljust(col_width)
        header_str += header[-2].ljust(col_width)
        header_str += header[-1].ljust(7)
        for r in range(trail["nrounds"] + 1):
            data_str += str(current_row) + '\t'
            data_str += trail.get(f"xl_{r}", 'none').ljust(col_width)
            data_str += trail.get(f"xr_{r}", 'none').ljust(col_width)
            data_str += trail.get(f"pr_{r}", 'none').ljust(col_width)
            data_str += '\n'
            current_row += 1
        stroutput = header_str
        stroutput += "\n" + "-"*len(header_str) + "\n"
        stroutput += data_str
        total_weight = trail["total_weight"]
        stroutput += f"Weight: -{total_weight}" + "\n"
        print(stroutput)
        return stroutput

    def find_characteristic(self):
        """
        Find the best linear trail for reduced-round SIMECK
        """
        diff_trail = None
        self.milp_model.Params.OutputFlag = False
        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        obj = self.milp_model.getObjective()
        # Consider the start_weight
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        time_start = time.time()
        #m.setParam(GRB.Param.Threads, 16)
        self.milp_model.optimize()
        # Gurobi syntax: m.Status == 2 represents the model is feasible.
        if (self.milp_model.Status == GRB.OPTIMAL or self.milp_model.Status == GRB.TIME_LIMIT or \
            self.milp_model.Status == GRB.INTERRUPTED):
            self.total_weight = self.milp_model.objVal
            print(f"\nThe squared correlation of the best linear characteristic: 2^-({self.total_weight})")
            print("\nLinear trail:\n")
            diff_trail = self.parse_solver_output()
            self.print_trail(trail=diff_trail)
        # Gurobi syntax: m.Status == 3 represents the model is infeasible. (GRB.Status.INFEASIBLE)
        elif self.milp_model.Status == GRB.INFEASIBLE:
            print("The model is infeasible!")
        else:
            print("Unknown error!")
        elapsed_time = time.time() - time_start
        print("Time used: %0.02f" % elapsed_time)
        return diff_trail

    def find_multiple_characteristics(self, number_of_trails=2):
        """
        Find multiple linear trails for reduced-round of SIMECK
        """
        self.milp_model.Params.PreSolve = 1
        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        obj = self.milp_model.getObjective()
        # Consider the start_weight
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        self.milp_model.Params.OutputFlag = False
        self.milp_model.Params.PoolSearchMode = 2
        # Limit number of solutions
        self.milp_model.Params.PoolSolutions = number_of_trails
        time_start = time.time()
        self.milp_model.optimize()
        if (self.milp_model.Status == GRB.OPTIMAL or self.milp_model.Status == GRB.TIME_LIMIT or \
            self.milp_model.Status == GRB.INTERRUPTED):
            # First Method:
            for sol_number in range(number_of_trails):
                if (self.milp_model.Status == GRB.OPTIMAL):
                    self.total_weight = self.milp_model.PoolObjVal
                    diff_trail = self.parse_solver_output()
                    self.print_trail(trail=diff_trail)
                elif (self.milp_model.Status == GRB.TIME_LIMIT or self.milp_model.Status == GRB.INTERRUPTED):
                    self.total_weight = self.milp_model.PoolObjVal
                    diff_trail = self.parse_solver_output()
                    self.print_trail(trail=diff_trail)
                    break
                else:
                    break
                self.exclude_the_previous_sol()
                print("#"*50)
                self.milp_model.optimize()
            # Second Method:
            # number_of_trails = self.milp_model.SolCount
            # for sol_number in range(number_of_trails):
            #     self.milp_model.Params.SolutionNumber = sol_number
            #     # PoolObjVal : This attribute is used to query the objective value of the <span>$</span>k<span>$</span>-th solution stored in the pool of feasible solutions found so far for the problem
            #     self.total_weight = self.milp_model.PoolObjVal
            #     diff_trail = self.parse_solver_output()
            #     self.print_trail(trail=diff_trail)
        # Gurobi syntax: m.Status == 3 represents the model is infeasible. (GRB.INFEASIBLE)
        elif self.milp_model.Status == GRB.INFEASIBLE:
            print("The model is infeasible!")
        else:
            print("Unknown error!")
        elapsed_time = time.time() - time_start
        print("Total time to find %s linear trails: %0.02f" % (number_of_trails, elapsed_time))

    def compute_linear_effect(self):
        """
        Compute the linear effect for a given input/output differences
        Some general information about Gurobi:
        PoolSolutions: It controls the size of the solution pool.
        Changing this parameter won't affect the number of solutions that are found -
        it simply determines how many of those are retained
        You can use the PoolSearchMode parameter to control the approach used to find solutions.
        In its default setting (0), the MIP search simply aims to find one optimal solution.
        Setting the parameter to 2 causes the MIP to lo a systematic search for the n best solutions.
        With a setting of 2, it will find the n best solutions,
        where n is determined by the value of the PoolSolutions parameter
        SolCount: Number of solutions found during the most recent optimization.

        Model status:
        LOADED	1	Model is loaded, but no solution information is available.
        OPTIMAL	2	Model was solved to optimality (subject to tolerances), and an optimal solution is available.
        INFEASIBLE	3	Model was proven to be infeasible.
        """

        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        #self.milp_model.Params.PreSolve = 0 # Activating this flag causes the performance to be decreased, but the accuracy will be increased
        self.milp_model.Params.PoolSearchMode = 2
        self.milp_model.Params.PoolSolutions = 1
        self.milp_model.Params.OutputFlag = False

        self.milp_model.printStats()

        # Consider the start_weight
        obj = self.milp_model.getObjective()
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        time_start = time.time()
        self.milp_model.optimize()
        current_probability = 0
        if (self.milp_model.Status == GRB.OPTIMAL):
            self.total_weight = self.milp_model.objVal
            diff_prob = 0
            print('\n')
            while (self.milp_model.Status == GRB.OPTIMAL and self.total_weight <= self.end_weight):
                self.total_weight = self.milp_model.PoolObjVal
                self.milp_model.Params.PoolSolutions = 2000000000 #GRB.MAXINT
                temp_constraint = self.milp_model.addConstr(obj == self.total_weight, name='temp_constraint')
                # self.milp_model.Params.PoolGap = 0
                # self.milp_model.Params.PreSolve = 0
                # self.milp_model.printStats()
                self.milp_model.update()
                self.milp_model.optimize()
                diff_prob += math.pow(2, -self.total_weight) * self.milp_model.SolCount
                print(f"Current weight: {self.total_weight}")
                print(f"Number of trails: {self.milp_model.SolCount}")
                current_probability = math.log(diff_prob, 2)
                print(f"\tCurrent Probability: 2^({current_probability})")
                elapsed_time = time.time() - time_start
                print("Time used = %0.04f seconds\n" % elapsed_time)
                self.milp_model.remove(temp_constraint)
                self.milp_model.Params.PoolSolutions = 1
                self.milp_model.addConstr(obj >= (self.total_weight + self.eps), name='temp_cond')
                #self.milp_model.Params.PreSolve = 0
                self.milp_model.optimize()
        elif (self.milp_model.Status == GRB.INFEASIBLE):
            print("The model is infeasible!")
        else:
            print("Unknown Error!")
        return current_probability

    def compute_differential_effect_classic_method(self):
        """
        Compute linear effect by enumerating all possible linear trails
        """

        if self.time_limit != None:
            self.milp_model.Params.TIME_LIMIT = self.time_limit
        self.milp_model.Params.OutputFlag = False
        # self.milp_model.printStats()
        # Consider the start_weight
        obj = self.milp_model.getObjective()
        if self.start_weight != None:
            self.milp_model.addConstr(obj >= self.start_weight, 'start_weight_constraint')
        time_start = time.time()
        self.milp_model.optimize()
        # self.milp_model.Params.Quad = 1
        sol_dict = dict()
        if (self.milp_model.Status == GRB.OPTIMAL):
            self.total_weight = self.milp_model.objVal
            diff_prob = 0
            print('\n')
            while (self.milp_model.Status == GRB.OPTIMAL and self.total_weight <= self.end_weight):
                self.total_weight = self.milp_model.objVal
                diff_prob += math.pow(2, -self.total_weight)
                total_weight_st = 'ntrails_%0.2f' % self.total_weight
                sol_dict[total_weight_st] = sol_dict.get(total_weight_st, 0) + 1
                print('Current weight: %s' % str(self.total_weight))
                print('Number of trails: %d' % sol_dict[total_weight_st])
                print('\tCurrent Probability: 2^(' + str(math.log(diff_prob, 2)) + ')')
                time_end = time.time()
                print('Time used = %0.4f seconds\n' % (time_end - time_start))
                self.exclude_the_previous_sol()
                self.milp_model.optimize()
        elif (self.milp_model.Status == GRB.INFEASIBLE):
            print('The model is infeasible!')
        else:
            print('Unknown Error!')

def loadparameters(args):
        """
        Get parameters from the argument list and inputfile.
        """
        # Load default values
        params = {"nrounds" : 1,
                  "blocksize" : 32,
                  "mode" : 0,
                  "startweight" : 0,
                  "endweight" : 128,
                  "timelimit" : 3600,
                  "numberoftrails" : 1,
                  "fixedVariables" : {}}
    
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
        if args.nrounds:
            params["nrounds"] = args.nrounds[0]

        if args.blocksize:
            params["blocksize"] = args.blocksize[0]

        if args.startweight:
            params["startweight"] = args.startweight[0]

        if args.endweight:
            params["endweight"] = args.endweight[0]

        if args.mode:
            params["mode"] = args.mode[0]

        if args.timelimit:
            params["timelimit"] = args.timelimit[0]

        if args.numberoftrails:
            params["numberoftrails"] = args.numberoftrails[0]
        return params

def main():
    """
    Parse the arguments and start the request functionality with the provided
    parameters.
    """

    parser = ArgumentParser(description="This tool finds the best linear"
                                        "trail in a cryptographic primitive"
                                        "using Gurobi",
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('--nrounds', type=int, nargs=1,
                        help="The number of rounds for the cipher")
    
    parser.add_argument('--blocksize', type=int, nargs=1,
                        choices=[32, 48, 64, 96, 128],
                        help="The blocksize of the cipher")
    
    parser.add_argument('--startweight', type=int, nargs=1,
                        help="Starting weight for the trail search.")
    parser.add_argument('--endweight', type=int, nargs=1,
                        help="Stop search after reaching endweight.")
    parser.add_argument('--mode', type=int, nargs=1,
                        choices=[0, 1, 2], help=
                        "0 = search characteristic for fixed round\n"
                        "1 = determine the probability of the linear\n")
    parser.add_argument('--timelimit', type=int, nargs=1,
                        help="Set a timelimit for the search in seconds.")
    parser.add_argument('--inputfile', help="Use an yaml input file to read the parameters.", nargs=1)
    parser.add_argument('--numberoftrails', type=int, nargs=1, 
                        help="Number of trails.")

    # Parse command line arguments and construct parameter list.
    args = parser.parse_args()
    params = loadparameters(args)
    present = Lin(params)
    present.make_model()
    present.solve()

if __name__ == "__main__":
    main()