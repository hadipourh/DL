/*
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
*/

% Strategy:
% Divide E into 4 parts: EU, EMU0, EMU1
%                            EML0, EML1, EL
% Use DDT (LAT) to encode the propagation of differential (resp. linear) trails through EU (resp. EL).
% Use *-DDT (*-LAT) to encode the propagation of differential (resp. linear) trails through EMU0 (resp. EML1).
% Enode the propagation of differential (resp. linear) trails through EMU1 (resp. EML0) with probability one.
% Minimize the probability of transtions through EU, EL as well as the number of common active S-boxes trhough
% the middle part, i.e., EMU0, EMU1, EML0, EML1.

include "eml.mzn";

int: blocksize;
int: RU;
int: RMU;
int: RM;
int: RML;
int: RL;
int: WU;
int: WM;
int: WL;
int: halfblocksize = blocksize div 2;

constraint assert(blocksize in {32, 48, 64}, "Invalid value for blocksize: " ++
        "blocksize must be 32, 48, or 64");
constraint assert(RU >= 0, "Invalid value for RU: " ++
        "RU must be greater than or equal to 0");
constraint assert(RMU >= 0, "Invalid value for RMU: " ++
        "RMU must be greater than or equal to 0");
constraint assert(RM >= 0, "Invalid value for RM: " ++
        "RM must be greater than or equal to 0");
constraint assert(RML >= 0, "Invalid value for RML: " ++
        "RML must be greater than or equal to 0");
constraint assert(RL >= 0, "Invalid value for RL: " ++
        "RL must be greater than or equal to 0");
constraint assert(RMU <= RM, "Invalid value for RMU: " ++
        "RMU must be less than RM");
constraint assert(RML <= RM, "Invalid value for RML: " ++ 
        "RML must be less than RM");

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%  ____          __  _                _                   _                    _____       _      _            
% |  _ \   ___  / _|(_) _ __    ___  | |     ___    ___  | | __ _   _  _ __   |_   _|__ _ | |__  | |  ___  ___ 
% | | | | / _ \| |_ | || '_ \  / _ \ | |    / _ \  / _ \ | |/ /| | | || '_ \    | | / _` || '_ \ | | / _ \/ __|
% | |_| ||  __/|  _|| || | | ||  __/ | |___| (_) || (_) ||   < | |_| || |_) |   | || (_| || |_) || ||  __/\__ \
% |____/  \___||_|  |_||_| |_| \___| |_____|\___/  \___/ |_|\_\ \__,_|| .__/    |_| \__,_||_.__/ |_| \___||___/
%                                                                     |_|                                      
% define lookup tables and constants

int: left_rotation_a0 = 1;
int: left_rotation_a1 = 5;
int: left_rotation_a2 = 0;

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____|| | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____| \___/ 
% Constraints for EU

array[0..RU, 0..<halfblocksize] of var 0..1: xu_left;
array[0..RU, 0..<halfblocksize] of var 0..1: xu_right;
array[0..<RU, 0..<halfblocksize] of var 0..1: and_output;
array[0..<RU, 0..<halfblocksize] of var 0..1: xor_output;
array[0..<RU, 0..<halfblocksize] of var 0..1: pu;

% exclude all-zero input
constraint sum(bit in 0..<halfblocksize)(xu_left[0, bit] + xu_right[0, bit]) != 0;

% constraints to model the round function
constraint forall(rn in 0..(RU - 1), bit in 0..<halfblocksize)(
    xu_left[rn, bit] = xu_right[rn + 1, bit]
    /\
    and_diff(
        xu_left[rn, (bit + left_rotation_a1) mod halfblocksize],
        xu_left[rn, (bit + left_rotation_a2) mod halfblocksize],
        and_output[rn, bit],
        pu[rn, bit]
    )
    /\
    xor2(
        xu_left[rn, (bit + left_rotation_a0) mod halfblocksize],
        and_output[rn, bit],
        xor_output[rn, bit]
    )
    /\
    xor2(
        xu_right[rn, bit],
        xor_output[rn, bit],
        xu_left[rn + 1, bit]
    )
    % We can also use the following sbox instead of using (and_diff, xor2, xor2)
    % diff_sbox(
    %     xu_left[rn, (bit + left_rotation_a0) mod halfblocksize],
    %     xu_left[rn, (bit + left_rotation_a1) mod halfblocksize],
    %     xu_left[rn, (bit + left_rotation_a2) mod halfblocksize],
    %     xu_right[rn, bit],
    %     xu_left[rn + 1, bit],
    %     pu[rn, bit]
    % )
);

var 0..256: PU;
constraint PU = 1*sum(rn in 0..<RU, bit in 0..<halfblocksize)(pu[rn, bit]);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  __  __  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____||  \/  || | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |\/| || | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |  | || |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_|  |_| \___/ 
% Constraints for EMU              

array[0..RM, 0..<halfblocksize] of var -1..1: xmu_left;
array[0..RM, 0..<halfblocksize] of var -1..1: xmu_right;

% link the input of EMU to the output of EU
constraint forall(bit in 0..<halfblocksize)(xmu_left[0, bit] = xu_left[RU, bit] /\ xmu_right[0, bit] = xu_right[RU, bit]);

% constraints to model the round function
constraint forall(rn in 0..<RM, bit in 0..<halfblocksize)(
    xmu_left[rn, bit] = xmu_right[rn + 1, bit]
    /\
    if rn < RMU then 
    (
        diff_sbox_star(
            xmu_left[rn, (bit + left_rotation_a0) mod halfblocksize],
            xmu_left[rn, (bit + left_rotation_a1) mod halfblocksize],
            xmu_left[rn, (bit + left_rotation_a2) mod halfblocksize],
            xmu_right[rn, bit],
            xmu_left[rn + 1, bit]
        )    
    ) else (
        diff_sbox_pr1(
            xmu_left[rn, (bit + left_rotation_a0) mod halfblocksize],
            xmu_left[rn, (bit + left_rotation_a1) mod halfblocksize],
            xmu_left[rn, (bit + left_rotation_a2) mod halfblocksize],
            xmu_right[rn, bit],
            xmu_left[rn + 1, bit]
        )
    ) endif
);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  __  __  _     
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____||  \/  || |    
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |\/| || |    
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |  | || |___ 
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_|  |_||_____|
% Constraints for EML                                                                                                     

array[0..RM, 0..<halfblocksize] of var -1..1: xml_left;
array[0..RM, 0..<halfblocksize] of var -1..1: xml_right;
array[0..<RML, 0..<halfblocksize] of var 0..1: mm0;
array[0..<RML, 0..<halfblocksize] of var 0..1: mm1;
array[0..<RML, 0..<halfblocksize] of var 0..1: mm2;

% link the output of EML the input of EL
constraint forall(bit in 0..<halfblocksize)(xml_left[RM, bit] = xl_left[0, bit] /\ xml_right[RM, bit] = xl_right[0, bit]);

% constraints to model the round function
int: offset = RM - RML;
constraint if blocksize = 32 then eml32(offset) elseif blocksize = 48 then eml48(offset) else eml64(offset) endif;
constraint forall(rn in offset..<RM, bit in 0..<halfblocksize)(
    xml_right[rn, bit] = mm0[rn - offset, (bit + left_rotation_a0) mod halfblocksize]
    /\
    fork4(mm0[rn - offset, bit], mm1[rn - offset, bit], mm2[rn - offset, bit], xml_right[rn + 1, bit], xml_left[rn, bit])
    /\
    xml_right[rn, bit] = xml_left[rn + 1, bit]
    /\
    and_linear_star(mm1[rn - offset, (bit + left_rotation_a1) mod halfblocksize], 
                    mm2[rn - offset, (bit + left_rotation_a2) mod halfblocksize], 
                    xml_right[rn, bit])
);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  _     
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____|| |    
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |    
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |___ 
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_____|
% constraints for EL

array[0..RL, 0..<halfblocksize] of var 0..1: xl_left;
array[0..RL, 0..<halfblocksize] of var 0..1: xl_right;
array[0..<RL, 0..<halfblocksize] of var 0..1: cl;
array[0..<RL, 0..<halfblocksize] of var 0..1: m0;
array[0..<RL, 0..<halfblocksize] of var 0..1: m1;
array[0..<RL, 0..<halfblocksize] of var 0..1: m2;

% exlude all-zero input 
constraint sum(bit in 0..<halfblocksize)(xl_left[0, bit] + xl_right[0, bit]) != 0;

% constraints to model the round function
constraint forall(rn in 0..<RL, bit in 0..<halfblocksize)(
    xl_right[rn, bit] = m0[rn, (bit + left_rotation_a0) mod halfblocksize]
    /\
    fork4(m0[rn, bit], m1[rn, bit], m2[rn, bit], xl_right[rn + 1, bit], xl_left[rn, bit])
    /\
    xl_right[rn, bit] = xl_left[rn + 1, bit]
    /\
    and_linear(m1[rn, (bit + left_rotation_a1) mod halfblocksize], 
                m2[rn, (bit + left_rotation_a2) mod halfblocksize], 
                xl_right[rn, bit], 
                cl[rn, bit])
);

var 0..256: CL;
constraint CL = 2*sum(rn in 0..<RL, bit in 0..<halfblocksize)(cl[rn, bit]);
% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%  ____                           _       ____   _                _                      
% / ___|   ___   __ _  _ __  ___ | |__   / ___| | |_  _ __  __ _ | |_  ___   __ _  _   _ 
% \___ \  / _ \ / _` || '__|/ __|| '_ \  \___ \ | __|| '__|/ _` || __|/ _ \ / _` || | | |
%  ___) ||  __/| (_| || |  | (__ | | | |  ___) || |_ | |  | (_| || |_|  __/| (_| || |_| |
% |____/  \___| \__,_||_|   \___||_| |_| |____/  \__||_|   \__,_| \__|\___| \__, | \__, |
%                                                                           |___/  |___/ 
% Search strategy
ann: search_ann;
% search_ann = seq_search([
%     int_search([yl[RL - 1, 0, i] | i in 0..63], indomain_min, complete),
%     int_search([yl[RL - 1, 1, i] | i in 0..63], indomain_min, complete),
%     int_search([yl[RL - 1, 2, i] | i in 0..63], indomain_min, complete),
%     int_search([yl[RL - 1, 3, i] | i in 0..63], indomain_min, complete),
%     int_search([yl[RL - 1, 3, i] | i in 0..63], indomain_min, complete)]);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ___   _      _              _    _               _____                     _    _               
%  / _ \ | |__  (_)  ___   ___ | |_ (_)__   __ ___  |  ___|_   _  _ __    ___ | |_ (_)  ___   _ __  
% | | | || '_ \ | | / _ \ / __|| __|| |\ \ / // _ \ | |_  | | | || '_ \  / __|| __|| | / _ \ | '_ \ 
% | |_| || |_) || ||  __/| (__ | |_ | | \ V /|  __/ |  _| | |_| || | | || (__ | |_ | || (_) || | | |
%  \___/ |_.__/_/ | \___| \___| \__||_|  \_/  \___| |_|    \__,_||_| |_| \___| \__||_| \___/ |_| |_|
%             |__/                                                                                  
% Objective function

% #############################################################################################################################################
var 0..4*256: CM;
constraint CM = 
                1*sum(rn in 0..<RM, bit in 0..<halfblocksize)(
                bool2int(xmu_left[rn, bit] == -1 /\ xml_left[rn, bit] == -1)
                +
                bool2int(xmu_left[rn, bit] == -1 /\ xml_left[rn, bit] == 1)
                +                 
                bool2int(xmu_left[rn, bit] == 1 /\ xml_left[rn, bit] == -1)
                +
                bool2int(xmu_right[rn, bit] == -1 /\ xml_right[rn, bit] == -1)
                +
                bool2int(xmu_right[rn, bit] == -1 /\ xml_right[rn, bit] == 1)
                +
                bool2int(xmu_right[rn, bit] == 1 /\ xml_right[rn, bit] == -1));


solve minimize (WU*PU + WM*CM + WL*CL);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%     _                 _  _  _                       _____                     _    _                    
%    / \   _   _ __  __(_)| |(_)  __ _  _ __  _   _  |  ___|_   _  _ __    ___ | |_ (_)  ___   _ __   ___ 
%   / _ \ | | | |\ \/ /| || || | / _` || '__|| | | | | |_  | | | || '_ \  / __|| __|| | / _ \ | '_ \ / __|
%  / ___ \| |_| | >  < | || || || (_| || |   | |_| | |  _| | |_| || | | || (__ | |_ | || (_) || | | |\__ \
% /_/   \_\\__,_|/_/\_\|_||_||_| \__,_||_|    \__, | |_|    \__,_||_| |_| \___| \__||_| \___/ |_| |_||___/
%                                             |___/       
% auxiliary functions

predicate xor2(var 0..1: a0, var 0..1: a1, var 0..1: b0) = 
    % Generate the constraints of a XOR with 2 inputssage: sb = SboxAnalyzer(sb)
    % sage: sb
    % (0, 1, 1, 0)
    % sage: cnf, milp = sb.minimized_diff_constraints()
    % Simplifying the MILP/SAT constraints ...
    % Time used to simplify the constraints: 0.00 seconds
    % Number of constraints: 4
    % Input:	a0||a1; a0: msb
    % Output:	b0; b0: msb
    % Weight: 0
    % sage: milp
    % ['- a0 + a1 + b0 >= 0',
    % 'a0 - a1 + b0 >= 0',
    % 'a0 + a1 - b0 >= 0',
    % '- a0 - a1 - b0 >= -2']
    - a0 + a1 + b0 >= 0
    /\
    a0 - a1 + b0 >= 0
    /\
    a0 + a1 - b0 >= 0
    /\
    - a0 - a1 - b0 >= -2;

predicate and_diff(var 0..1: a0, var 0..1: a1, var 0..1: b0, var 0..1: p0) = 
    % Generate the constraints of a AND with 2 inputs
    % sage: sb = [0, 0, 0, 1]
    % sage: sb = SboxAnalyzer(sb)
    % sage: cnf, milp = sb.minimized_diff_constraints()
    % Simplifying the MILP/SAT constraints ...
    % Time used to simplify the constraints: 0.00 seconds
    % Number of constraints: 4
    % Input:	a0||a1; a0: msb
    % Output:	b0; b0: msb
    % Weight: 1.0000 p0
    % sage: milp
    % ['a0 + a1 - b0 >= 0', 'a0 + a1 - p0 >= 0', '- a1 + p0 >= 0', '- a0 + p0 >= 0']
    a0 + a1 - b0 >= 0
    /\
    a0 + a1 - p0 >= 0
    /\
    - a1 + p0 >= 0
    /\
    - a0 + p0 >= 0;

predicate fork4(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3, var 0..1: b0) = 
    % Generate the constraints of a Fork with 4 outputs
    % The conditions are the same as the linear constraints for XOR with 4 inputs
    % sage: sb = [0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0]
    % sage: sa = SboxAnalyzer(sb)
    % sage: ddt
    % [16  0]
    % [ 0 16]
    % [ 0 16]
    % [16  0]
    % [ 0 16]
    % [16  0]
    % [16  0]
    % [ 0 16]
    % [ 0 16]
    % [16  0]
    % [16  0]
    % [ 0 16]
    % [16  0]
    % [ 0 16]
    % [ 0 16]
    % [16  0]
    % sage: cnf, milp = sa.minimized_diff_constraints()
    % Simplifying the MILP/SAT constraints ...
    % Time used to simplify the constraints: 0.00 seconds
    % Number of constraints: 16
    % Input:	a0||a1||a2||a3; a0: msb
    % Output:	b0; b0: msb
    % Weight: 0

    % ['- a0 + a1 + a2 + a3 + b0 >= 0',
    % 'a0 - a1 + a2 + a3 + b0 >= 0',
    % 'a0 + a1 - a2 + a3 + b0 >= 0',
    % '- a0 - a1 - a2 + a3 + b0 >= -2',
    % 'a0 + a1 + a2 - a3 + b0 >= 0',
    % '- a0 - a1 + a2 - a3 + b0 >= -2',
    % '- a0 + a1 - a2 - a3 + b0 >= -2',
    % 'a0 - a1 - a2 - a3 + b0 >= -2',
    % 'a0 + a1 + a2 + a3 - b0 >= 0',
    % '- a0 - a1 + a2 + a3 - b0 >= -2',
    % '- a0 + a1 - a2 + a3 - b0 >= -2',
    % 'a0 - a1 - a2 + a3 - b0 >= -2',
    % '- a0 + a1 + a2 - a3 - b0 >= -2',
    % 'a0 - a1 + a2 - a3 - b0 >= -2',
    % 'a0 + a1 - a2 - a3 - b0 >= -2',
    % '- a0 - a1 - a2 - a3 - b0 >= -4']
    - a0 + a1 + a2 + a3 + b0 >= 0
    /\
    a0 - a1 + a2 + a3 + b0 >= 0
    /\
    a0 + a1 - a2 + a3 + b0 >= 0
    /\
    - a0 - a1 - a2 + a3 + b0 >= -2
    /\
    a0 + a1 + a2 - a3 + b0 >= 0
    /\
    - a0 - a1 + a2 - a3 + b0 >= -2
    /\
    - a0 + a1 - a2 - a3 + b0 >= -2
    /\
    a0 - a1 - a2 - a3 + b0 >= -2
    /\
    a0 + a1 + a2 + a3 - b0 >= 0
    /\
    - a0 - a1 + a2 + a3 - b0 >= -2
    /\
    - a0 + a1 - a2 + a3 - b0 >= -2
    /\
    a0 - a1 - a2 + a3 - b0 >= -2
    /\
    - a0 + a1 + a2 - a3 - b0 >= -2
    /\
    a0 - a1 + a2 - a3 - b0 >= -2
    /\
    a0 + a1 - a2 - a3 - b0 >= -2
    /\
    - a0 - a1 - a2 - a3 - b0 >= -4;

predicate diff_sbox(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3, var 0..1: b0, var 0..1: p0) = 
    % sage: sb = SBox([0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1])
    % sage: sa = SboxAnalyzer(sb)
    % sage: cnf, milp = sa.minimized_diff_constraints()
    % Simplifying the MILP/SAT constraints ...
    % Time used to simplify the constraints: 0.00 seconds
    % Number of constraints: 7
    % Input:	a0||a1||a2||a3; a0: msb
    % Output:	b0; b0: msb
    % Weight: 1.0000 p0
    - a0 + a1 + a2 - a3 - b0 >= -2
    /\
    a0 + a1 + a2 + a3 - b0 >= 0
    /\
    a0 + a1 + a2 - a3 + b0 >= 0
    /\
    - a0 + a1 + a2 + a3 + b0 >= 0
    /\
    - a2 + p0 >= 0
    /\
    - a1 + p0 >= 0
    /\
    a1 + a2 - p0 >= 0;

predicate diff_sbox_star(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3, var 0..1: b0) = 
    % sage: cnf, milp = sa.minimized_diff_constraints(subtable='star')
    % Simplifying the MILP/SAT constraints ...
    % product term(s) span more than one line (warning only)
    % Time used to simplify the constraints: 0.00 seconds
    % Number of constraints: 4
    % Input:	a0||a1||a2||a3; a0: msb
    % Output:	b0; b0: msb
    % sage: milp
    % ['- a0 + a1 + a2 + a3 + b0 >= 0',
    % 'a0 + a1 + a2 - a3 + b0 >= 0',
    % 'a0 + a1 + a2 + a3 - b0 >= 0',
    % '- a0 + a1 + a2 - a3 - b0 >= -2']
    - a0 + a1 + a2 + a3 + b0 >= 0
    /\
    a0 + a1 + a2 - a3 + b0 >= 0
    /\
    a0 + a1 + a2 + a3 - b0 >= 0
    /\
    - a0 + a1 + a2 - a3 - b0 >= -2;

predicate diff_sbox_pr1(var -1..1: x3, var -1..1: x2, var -1..1: x1, var -1..1: x0, var -1..1: y0) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y0 = 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y0 = 1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y0 = 1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y0 = 0)
    else (y0 = -1)
    endif;


predicate and_linear(var 0..1: a0, var 0..1: a1, var 0..1: b0, var 0..1: p0) = 
    % sage: sb = [0, 0, 0, 1]
    % sage: sa = SboxAnalyzer(sb)
    % sage: lat = sa.linear_approximation_table(scale='correlation')
    % sage: lat
    % [   1  1/2]
    % [   0  1/2]
    % [   0  1/2]
    % [   0 -1/2]
    % sage: cnf, milp = sa.minimized_linear_constraints()
    % Simplifying the MILP/SAT constraints ...
    % Time used to simplify the constraints: 0.00 seconds
    % Number of constraints: 4
    % Input:	a0||a1; a0: msb
    % Output:	b0; b0: msb
    % Weight: 2.0000 p0
    % sage: milp
    % ['- a0 + p0 >= 0', '- a1 + p0 >= 0', '- b0 + p0 >= 0', 'b0 - p0 >= 0']
    - a0 + p0 >= 0
    /\
    - a1 + p0 >= 0
    /\
    - b0 + p0 >= 0
    /\
    b0 - p0 >= 0;

predicate and_linear_star(var 0..1: a0, var 0..1: a1, var 0..1: b0) = 
    % sage: sb = [0, 0, 0, 1]
    % sage: sa = SboxAnalyzer(sb)
    % sage: lat = sa.linear_approximation_table(scale='correlation')
    % sage: lat
    % [   1  1/2]
    % [   0  1/2]
    % [   0  1/2]
    % [   0 -1/2]
    % sage: cnf, milp = sa.minimized_linear_constraints(subtable='star')
    % Simplifying the MILP/SAT constraints ...
    % Time used to simplify the constraints: 0.00 seconds
    % Number of constraints: 2
    % Input:	a0||a1; a0: msb
    % Output:	b0; b0: msb
    % sage: milp
    % ['- a1 + b0 >= 0', '- a0 + b0 >= 0']
    - a1 + b0 >= 0
    /\
    - a0 + b0 >= 0;

predicate eml32(int: nrounds) = 
	let {
		array[0..<nrounds, 0..<16] of var -1..1: ydl;
		array[0..<nrounds, 0..<16] of var -1..1: zdl;
		array[0..<nrounds, 0..<16] of var -1..1: wdl;
		array[0..<nrounds, 0..<16] of var -1..1: zdl1;
		array[0..<nrounds, 0..<16] of var -1..1: pd;
	} in 
    % Branching
    forall(round in 0..<nrounds, i in 0..<16)
    (
        branching(xml_right[round + 1, i], zdl1[round, i], ydl[round, i])
    )
    /\
    % Diffusion layer
    forall(round in 0..<nrounds, i in 0..<16)
    (
        zdl[round, i] = zdl1[round, left_rotate_1_32[i]]
    )
    /\
    % Xor
    forall(round in 0..<nrounds, i in 0..<16)
    (
        xml_left[round + 1, i] = zdl[round, i]
    )
    /\
    % S-box layer
    forall(round in 0..<nrounds, i in 0..4)
    (
        sbox(ydl[round, i], ydl[round, 5 + i], ydl[round, 10 + i], xml_left[round + 1, i], xml_left[round + 1, 5 + i],
             pd[round, i], xml_left[round, 5 + i], pd[round, 10 + i], wdl[round, i], wdl[round, 5 + i])
    )
    /\
    forall(round in 0..<nrounds)
    (
        sbox(pd[round, 10], ydl[round, 15], pd[round, 4], xml_left[round + 1, 10], xml_left[round + 1, 15],
             xml_left[round, 10], xml_left[round, 15], xml_left[round, 4], wdl[round, 10], wdl[round, 15])
    )
    /\
    forall(round in 0..<nrounds, i in 0..3)
    (
        sbox1(pd[round, i], pd[round, 11 + i], xml_left[round + 1, 11 + i], xml_left[round, i], xml_left[round, 11 + i], wdl[round, 11 + i])
    )
    /\
    % Xor
    forall(round in 0..<nrounds, i in 0..<16)
    (
        zdl[round, i] = wdl[round, i]
    )
    /\
    forall(round in 0..<nrounds, i in 0..<16)
    (
        xml_right[round, i] = wdl[round, i]
    );

predicate eml48(int: nrounds) = 
   let {
        array[0..<nrounds, 0..<24] of var -1..1: ydl;
        array[0..<nrounds, 0..<24] of var -1..1: zdl;
        array[0..<nrounds, 0..<24] of var -1..1: wdl;
        array[0..<nrounds, 0..<24] of var -1..1: zdl1;
        array[0..<nrounds, 0..<24] of var -1..1: pd;
   } in 
   % Branching
   forall(round in 0..<nrounds, i in 0..<24)
   (
   branching(xml_right[round + 1, i], zdl1[round, i], ydl[round, i])
   ) 
   /\
   % Diffusion layer
   forall(round in 0..<nrounds, i in 0..<24)
   (
      zdl[round, i] = zdl1[round, left_rotate_1_48[i]]
   )
   /\
   % Xor
   forall(round in 0..<nrounds, i in 0..<24)
   (
       xml_left[round + 1, i] = zdl[round, i]
   )
   /\
    % S-box layer
    forall(round in 0..<nrounds, i in 0..4)
    (
        sbox(ydl[round, i], ydl[round, 5 + i], ydl[round, 10 + i], xml_left[round + 1, i], xml_left[round + 1, 5 + i],
            pd[round, i], xml_left[round, 5 + i], pd[round, 10 + i], wdl[round, i], wdl[round, 5 + i])
    )
    /\
    forall(round in 0..<nrounds, i in 0..3)
    (
        sbox(pd[round, 10 + i], ydl[round, 15 + i], ydl[round, 20 + i], xml_left[round + 1, 10 + i], xml_left[round + 1, 15 + i],
            xml_left[round, 10 + i], xml_left[round, 15 + i], pd[round, 20 + i], wdl[round, 10 + i], wdl[round, 15 + i])
    )
    /\
    forall(round in 0..<nrounds)
    (
        sbox(pd[round, 14], ydl[round, 19], pd[round, 0], xml_left[round + 1, 14], xml_left[round + 1, 19],
            xml_left[round, 14], xml_left[round, 19], xml_left[round, 0], wdl[round, 14], wdl[round, 19])
    )
    /\
    forall(round in 0..<nrounds, i in 1..4)
    (
        sbox1(pd[round, i], pd[round, 19 + i], xml_left[round + 1, 19 + i], xml_left[round, i], xml_left[round, 19 + i], wdl[round, 19 + i])
    )
    /\
    % Xor
    forall(round in 0..<nrounds, i in 0..<24)
    (
        zdl[round, i] = wdl[round, i]
    )
    /\
    forall(round in 0..<nrounds, i in 0..<24)
    (
        xml_right[round, i] = wdl[round, i]
    );

predicate eml64(int: nrounds) = 
    let {
        array[0..<nrounds, 0..<32] of var -1..1: ydl;
        array[0..<nrounds, 0..<32] of var -1..1: zdl;
        array[0..<nrounds, 0..<32] of var -1..1: wdl;
        array[0..<nrounds, 0..<32] of var -1..1: zdl1;
        array[0..<nrounds, 0..<32] of var -1..1: pd;
    } in 
    % Branching
    forall(round in 0..<nrounds, i in 0..<32)
    (
    branching(xml_right[round + 1, i], zdl1[round, i], ydl[round, i])
    )
    /\
    % Diffusion layer
    forall(round in 0..<nrounds, i in 0..<32)
    (
        zdl[round, i] = zdl1[round, left_rotate_1_64[i]]
    )
    /\
    % Xor
    forall(round in 0..<nrounds, i in 0..<32)
    (
        xml_left[round + 1, i] = zdl[round, i]  
    )
    /\
    % S-box layer
    forall(round in 0..<nrounds, i in 0..4)
    (
        sbox(ydl[round, i], ydl[round, 5 + i], ydl[round, 10 + i], xml_left[round + 1, i], xml_left[round + 1, 5 + i],
            pd[round, i], xml_left[round, 5 + i], pd[round, 10 + i], wdl[round, i], wdl[round, 5 + i])
    )
    /\
    forall(round in 0..<nrounds, i in 0..4)
    (
        sbox(pd[round, 10 + i], ydl[round, 15 + i], ydl[round, 20 + i], xml_left[round + 1, 10 + i], xml_left[round + 1, 15 + i],
            xml_left[round, 10 + i], xml_left[round, 15 + i], pd[round, 20 + i], wdl[round, 10 + i], wdl[round, 15 + i])
    )
    /\
    forall(round in 0..<nrounds, i in 0..1)
    (
        sbox(pd[round, 20 + i], ydl[round, 25 + i], ydl[round, 30 + i], xml_left[round + 1, 20 + i], xml_left[round + 1, 25 + i],
            xml_left[round, 20 + i], xml_left[round, 25 + i], pd[round, 30 + i], wdl[round, 20 + i], wdl[round, 25 + i])
    )
    /\
    forall(round in 0..<nrounds, i in 0..2)
    (
        sbox(pd[round, 22 + i], ydl[round, 27 + i], pd[round, i], xml_left[round + 1, 22 + i], xml_left[round + 1, 27 + i],
            xml_left[round, 22 + i], xml_left[round, 27 + i], xml_left[round, i], wdl[round, 22 + i], wdl[round, 27 + i])
    )
    /\
    forall(round in 0..<nrounds, i in 0..1)
    (
        sbox1(pd[round, 30 + i], pd[round, 3 + i], xml_left[round + 1, 30 + i], xml_left[round, 30 + i], xml_left[round, 3 + i], wdl[round, 30 + i])
    )
    /\
    % Xor
    forall(round in 0..<nrounds, i in 0..<32)
    (    
        zdl[round, i] = wdl[round, i]  
    )
    /\
    forall(round in 0..<nrounds, i in 0..<32)
    (
        xml_right[round, i] = wdl[round, i]  
    );