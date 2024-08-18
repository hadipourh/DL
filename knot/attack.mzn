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

int: RU;
int: RMU;
int: RM;
int: RML;
int: RL;
int: offset;
int: is_limited;
% nc: number of columns in each row of state array. For KNOT-(4*nc) it is nc
int: nc;

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
        "RMU must be less than or equal to RM");
constraint assert(RML <= RM, "Invalid value for RML: " ++ 
        "RML must be less than or equal to RM");
% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%  ____          __  _                _                   _                    _____       _      _            
% |  _ \   ___  / _|(_) _ __    ___  | |     ___    ___  | | __ _   _  _ __   |_   _|__ _ | |__  | |  ___  ___ 
% | | | | / _ \| |_ | || '_ \  / _ \ | |    / _ \  / _ \ | |/ /| | | || '_ \    | | / _` || '_ \ | | / _ \/ __|
% | |_| ||  __/|  _|| || | | ||  __/ | |___| (_) || (_) ||   < | |_| || |_) |   | || (_| || |_) || ||  __/\__ \
% |____/  \___||_|  |_||_| |_| \___| |_____|\___/  \___/ |_|\_\ \__,_|| .__/    |_| \__,_||_.__/ |_| \___||___/
%                                                                     |_|                                      
% define lookup tables

array[0..3] of int: rotation_dictionary = array1d(0..3,[0, 1, 8, 25]);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____|| | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____| \___/ 
% Constraints for EU

array[0..RU, 0..3, 0..<nc] of var 0..1: xu;
array[0..(RU - 1), 0..3, 0..<nc] of var 0..1: yu;
array[0..(RU - 1), 0..1, 0..<nc] of var 0..1: pu;

% S-box layer
constraint forall(round in 0..(RU - 1), column in 0..<nc)
(
    sbox_diff(xu[round, 0, column], xu[round, 1, column], xu[round, 2, column], xu[round, 3, column],
              yu[round, 0, column], yu[round, 1, column], yu[round, 2, column], yu[round, 3, column],
              pu[round, 0, column], pu[round, 1, column])
);

% Diffusion layer
constraint forall(round in 0..(RU - 1), row in 0..3, column in 0..<nc)
(
   xu[round + 1, row, column] = yu[round, row, ((column + rotation_dictionary[row]) mod nc)]
);

constraint sum(i in 0..3, j in 0..<nc) (xu[0, i, j]) != 0;

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  __  __  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____||  \/  || | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |\/| || | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |  | || |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_|  |_| \___/ 
% Constraints for EMU                                                                                                                                                             

array[0..RM, 0..3, 0..<nc] of var -1..1: xmu;
array[0..(RM - 1), 0..3, 0..<nc] of var -1..1: ymu;
array[0..(RM - 1), 0..<nc] of var 0..1: dmu;

% link the input of EM to the output of EU
constraint forall(row in 0..3, column in 0..<nc)(xmu[0, row, column] = xu[RU, row, column]);

% S-box layer
constraint forall(round in 0..(RM - 1), column in 0..<nc)
(
    if round < RMU then 
    (
        sbox_diff_star(xmu[round, 0, column], xmu[round, 1, column], xmu[round, 2, column], xmu[round, 3, column],
                       ymu[round, 0, column], ymu[round, 1, column], ymu[round, 2, column], ymu[round, 3, column],
                       dmu[round, column])
    ) else 
    (
        sbox_diff_m(xmu[round, 0, column], xmu[round, 1, column], xmu[round, 2, column], xmu[round, 3, column],
                    ymu[round, 0, column], ymu[round, 1, column], ymu[round, 2, column], ymu[round, 3, column],
                    dmu[round, column])
    ) endif
);

% Diffusion layer
constraint forall(round in 0..(RM - 1), row in 0..3, column in 0..<nc)
(
    xmu[round + 1, row, column] = ymu[round, row, ((column + rotation_dictionary[row]) mod nc)]
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

array[0..RM, 0..3, 0..<nc] of var -1..1: xml;
array[0..(RM - 1), 0..3, 0..<nc] of var -1..1: yml;
array[0..(RM - 1), 0..<nc] of var 0..1: dml;

% S-box layer
constraint forall(round in 0..(RM - 1), column in 0..<nc)
(
    if round < (RM - RML) then 
    (
        sbox_lin_m(yml[round, 0, column], yml[round, 1, column], yml[round, 2, column], yml[round, 3, column],
                   xml[round, 0, column], xml[round, 1, column], xml[round, 2, column], xml[round, 3, column],
                   dml[round, column])
    ) else
    (
        sbox_lin_star(yml[round, 0, column], yml[round, 1, column], yml[round, 2, column], yml[round, 3, column],
                      xml[round, 0, column], xml[round, 1, column], xml[round, 2, column], xml[round, 3, column],
                      dml[round, column])
    ) endif
);

% Diffusion layer
constraint forall(round in 0..(RM - 1), row in 0..3, column in 0..<nc)
(
    xml[round + 1, row, column] = yml[round, row, ((column + rotation_dictionary[row]) mod nc)]
    % yml[round, row, column] = xml[round + 1, row, ((column - rotation_dictionary[row] + nc) mod nc)]
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

array[0..RL, 0..3, 0..<nc] of var 0..1: xl;
array[0..(RL - 1), 0..3, 0..<nc] of var 0..1: yl;
array[0..(RL - 1), 0..1, 0..<nc] of var 0..1: pl;

% link the output of EML to the input of EL
constraint forall(row in 0..3, column in 0..<nc)(xl[0, row, column] = xml[RM, row, column]);

% S-box layer
constraint forall(round in 0..(RL - 1), column in 0..<nc)
(
    sbox_lin(yl[round, 0, column], yl[round, 1, column], yl[round, 2, column], yl[round, 3, column],
             xl[round, 0, column], xl[round, 1, column], xl[round, 2, column], xl[round, 3, column],
             pl[round, 0, column], pl[round, 1, column])
);

% Diffusion layer
constraint forall(round in 0..(RL - 1), row in 0..3, column in 0..<nc)
(
    xl[round + 1, row, column] = yl[round, row, ((column + rotation_dictionary[row]) mod nc)]
    % yl[round, row, column] = xl[round + 1, row, ((column - rotation_dictionary[row] + nc) mod nc)]
);

constraint sum(row in 0..3, column in 0..<nc)(xl[RL, row, column]) != 0;

% #############################################################################################################################################
% #############################################################################################################################################
var 0..(4*nc + 5000): PU;
var 0..(4*nc + 5000): CM;
var 0..(4*nc + 5000): CL;
var 0..(4*nc + 5000): NASM;
constraint PU = 3*sum(r in offset..(RU - 1), column in 0..<nc)(pu[r, 0, column]) + 
                2*sum(r in offset..(RU - 1), column in 0..<nc)(pu[r, 1, column]);

constraint CM = 1*sum(r in RM..(RM), row in 0..3, column in 0..<nc)(
                bool2int(xmu[r, row, column] = -1 /\ xml[r, row, column] = -1)
                +                
                bool2int(xmu[r, row, column] = -1 /\ xml[r, row, column] = 1)
                +                 
                bool2int(xmu[r, row, column] = 1 /\ xml[r, row, column] = -1));
    
% constraint CM = sum(r in 0..(RM - 1), row in 0..3, column in 0..<nc)(
%                 bool2int((ymu[r, row, column] != 0 /\ yml[r, row, column] != 0)));

constraint NASM = 2*sum(r in 0..(RM - 1), column in 0..<nc)(bool2int(dmu[r, column] + dml[r, column] == 2));
constraint CL = 4*sum(r in 0..(RL - 1), column in 0..<nc)(pl[r, 0, column]) +               
                2*sum(r in 0..(RL - 1), column in 0..<nc)(pl[r, 1, column]);

constraint if(is_limited == 1) then
(
    forall(row in 2..3, column in 0..<nc)(xu[0, row, column] = 0)
    /\
    forall(row in 1..3, column in 0..<nc)(xl[RL, row, column] = 0)
) else
    true
endif;
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
%     int_search([yl[RL - 1, 0, i] | i in 0..<nc], indomain_min, complete),
%     int_search([yl[RL - 1, 1, i] | i in 0..<nc], indomain_min, complete),
%     int_search([yl[RL - 1, 2, i] | i in 0..<nc], indomain_min, complete),
%     int_search([yl[RL - 1, 3, i] | i in 0..<nc], indomain_min, complete),
%     int_search([yl[RL - 1, 3, i] | i in 0..<nc], indomain_min, complete)]);

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

% constraint forall(row in 2..3, column in 0..<nc)(xu[0, row, column] = 0);
% constraint forall(column in 0..<nc)((xu[0, 0, column] = 1 /\ xu[0, 3, column] = 1) \/ (xu[0, 0, column] = 0 /\ xu[0, 3, column] = 0));
% constraint forall(row in 1..3, column in 0..<nc)(xl[RL, row, column] = 0);

% solve :: restart_constant(10000)
%       :: search_ann minimize (PU + NASM + CL);

solve minimize (nc*PU + CM + nc*CL);
% solve minimize (PU + CM + CL);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################

% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xu[r, i, j]) = 1 then "1" else show(fix(xu[r, i, j])) endif | j in 0..<nc]) ++ "\"\n"
%     | r in 0..0, i in 0..3];

% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "yu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(yu[r, i, j]) = 1 then "1" else show(fix(yu[r, i, j])) endif | j in 0..<nc]) ++ "\"\n"
%     | r in 0..0 where RU > 0, i in 0..3];

% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xu[r, i, j]) = 1 then "1" else show(fix(xu[r, i, j])) endif | j in 0..<nc]) ++ "\"\n"
%     | r in 1..RU, i in 0..3];

% output ["\n"];
% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xmu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xmu[r, i, j]) = -1 then "?" else show(fix(xmu[r, i, j])) endif | j in 0..<nc]) ++ "\"\n"
%     | r in 0..RM, i in 0..3];

% output ["\n"];
% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xml_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xml[r, i, j]) = -1 then "?" else show(fix(xml[r, i, j])) endif | j in 0..<nc]) ++ "\"\n"
%     | r in 0..RM, i in 0..3];

% output ["\n"];
% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xl_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xl[r, i, j]) = 1 then "1" else show(fix(xl[r, i, j])) endif | j in 0..<nc]) ++ "\"\n"
%     | r in 0..RL, i in 0..3];

% output ["\n\nCost of EU:                " ++ show(PU)];
% output ["\n\nNo. common active S-boxes: " ++ show(NASM)];
% output ["\n\nNo. common active bits:    " ++ show(CM)];
% output ["\n\nCost of EL:                " ++ show(CL)];

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

% We use S-box Analyzer to encode DDT: 
% https://github.com//sboxanalyzer
% Simplifying the MILP/SAT constraints ...
% Time used to simplify the constraints: 0.01 seconds
% Number of constraints: 37
% Input:	a0||a1||a2||a3; a0: msb
% Output:	b0||b1||b2||b3; b0: msb
% Weight: 3.0000 p0 + 2.0000 p1
predicate sbox_diff(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                    var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                    var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - a0 - a2 + b0 + p0 >= -1
    /\
    a0 - b0 - b1 + p0 >= -1
    /\
    - a3 - b0 - b1 + p0 >= -2
    /\
    a2 + b0 - b1 + p0 >= 0
    /\
    a0 - a2 + b1 + p0 >= 0
    /\
    - a2 - b0 - b3 + p0 >= -2
    /\
    a2 + b1 - b2 + p1 >= 0
    /\
    - a0 + b1 + b2 + p1 >= 0
    /\
    - a0 - a1 + a2 + a3 - b0 >= -2
    /\
    a0 - a1 + a2 + a3 + b0 >= 0
    /\
    a0 + a1 + a2 - a3 + b1 >= 0
    /\
    a0 - a1 - a2 - b1 - b2 >= -3
    /\
    a0 + a1 - a2 - b1 + b2 >= -1
    /\
    a0 + a2 + b1 + b2 - b3 >= 0
    /\
    - a0 - b0 + b1 - b2 + b3 >= -2
    /\
    - a1 - a3 + b0 + b1 + p0 >= -1
    /\
    b0 + b1 - b2 - b3 - p1 >= -2
    /\
    - a1 + a2 - a3 - b0 + p1 >= -2
    /\
    - a0 + a1 + a2 - b1 + p1 >= -1
    /\
    a1 + a2 + a3 - b1 + p1 >= 0
    /\
    a0 - a2 + b0 - b1 + p1 >= -1
    /\
    - a0 - a3 - b0 - b3 + p1 >= -3
    /\
    - b0 + b1 - b2 - b3 + p1 >= -2
    /\
    - b0 + b1 + b2 + b3 + p1 >= 0
    /\
    a1 + a3 + b0 - b1 + b2 - b3 >= -1
    /\
    - a0 - a2 + a3 - b0 - b1 + b3 >= -3
    /\
    a1 - a2 - a3 + b0 + b2 + b3 >= -1
    /\
    a1 + a2 + a3 + b1 + p0 - p1 >= 0
    /\
    a2 + b1 + b2 + b3 + p0 - p1 >= 0
    /\
    a0 + b0 + b1 + b3 - p0 + p1 >= 0
    /\
    a1 - a2 - a3 + b0 - b1 - b2 - b3 >= -4
    /\
    - a1 - a2 - a3 + b0 - b1 + b2 - b3 >= -4
    /\
    a1 - a2 + a3 + b0 - b1 - b2 + b3 >= -2
    /\
    - a1 - a2 + a3 + b0 - b1 + b2 + b3 >= -2
    /\
    - a1 - a2 - a3 + b0 - b1 - b2 + b3 >= -4
    /\
    - a1 - a2 + a3 + b0 - b1 - b2 - b3 >= -4
;

% We use S-box Analyzer to encode DDT: 
% https://github.com//sboxanalyzer
% We use S-box Analyzer to encode *-DDT:
% Simplifying the MILP/SAT constraints ...
% Time used to simplify the constraints: 0.00 seconds
% Number of constraints: 30
% Input:	a0||a1||a2||a3; a0: msb
% Output:	b0||b1||b2||b3; b0: msb
predicate sbox_diff_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                         var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                         var 0..1: d) =
    - a0 - a1 + a2 + a3 - b0 >= -2
    /\
    a0 + a1 + a2 + a3 - b0 >= 0
    /\
    - a0 + a1 + a2 + b0 - b1 >= -1
    /\
    a0 + a2 + a3 + b0 - b1 >= 0
    /\
    a0 + a1 + a2 - a3 + b1 >= 0
    /\
    - a0 + a1 + a2 + a3 + b1 >= 0
    /\
    - a1 + a2 - a3 + b0 + b1 >= -1
    /\
    a0 - a1 - a2 - b1 - b2 >= -3
    /\
    a0 + a1 - a2 - b1 + b2 >= -1
    /\
    - a0 - a3 - b0 - b1 - b3 >= -4
    /\
    - a2 - b0 + b1 - b2 - b3 >= -3
    /\
    a2 + b0 + b1 - b2 - b3 >= -1
    /\
    - a0 - a2 + b1 + b2 - b3 >= -2
    /\
    a0 + a2 + b1 + b2 - b3 >= 0
    /\
    - a0 + a2 - a3 - b0 + b3 >= -2
    /\
    - a0 - b0 + b1 - b2 + b3 >= -2
    /\
    a0 + b0 + b1 - b2 + b3 >= 0
    /\
    - a1 + a2 + b1 + b2 + b3 >= 0
    /\
    a0 - b0 + b1 + b2 + b3 >= 0
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    a0 - a1 + a2 - a3 - b0 - b1 >= -3
    /\
    a1 + a3 + b0 - b1 + b2 - b3 >= -1
    /\
    - a0 - a2 + a3 - b0 - b1 + b3 >= -3
    /\
    - a0 + a1 - a3 + b0 + b2 + b3 >= -1
    /\
    a1 - a2 - a3 + b0 - b1 - b2 - b3 >= -4
    /\
    - a1 - a2 - a3 + b0 - b1 + b2 - b3 >= -4
    /\
    a1 - a2 + a3 + b0 - b1 - b2 + b3 >= -2
    /\
    - a1 - a2 + a3 + b0 - b1 + b2 + b3 >= -2
    /\
    - a1 - a2 - a3 + b0 - b1 - b2 + b3 >= -4
    /\
    - a1 - a2 + a3 + b0 - b1 - b2 - b3 >= -4
    /\
    if exists([a0 == 1, a1 == 1, a2 == 1, a3 == 1]) then d = 1
    else d = 0 endif
;

% We use S-box Analyzer to encode LAT: 
% https://github.com//sboxanalyzer
% Simplifying the MILP/SAT constraints ...
% Time used to simplify the constraints: 0.01 seconds
% Number of constraints: 34
% Input:	a0||a1||a2||a3; a0: msb
% Output:	b0||b1||b2||b3; b0: msb
% Weight: 4.0000 p0 + 2.0000 p1
predicate sbox_lin(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                   var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                   var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - a3 - b0 - b1 + p0 >= -2
    /\
    - a2 + b1 - b3 + p0 >= -1
    /\
    - a3 - b1 + b3 + p0 >= -1
    /\
    - a1 + a2 + a3 + p1 >= 0
    /\
    a2 + a3 - b0 + p1 >= 0
    /\
    a2 + b1 - b3 + p1 >= 0
    /\
    - a2 + a3 + b3 + p1 >= 0
    /\
    - a3 + b1 + b3 + p1 >= 0
    /\
    a2 + a3 - p0 + p1 >= 0
    /\
    - a0 - a1 + a2 + a3 - b0 >= -2
    /\
    - a0 + a1 + a2 + a3 + b0 >= 0
    /\
    - a0 + a2 - b0 + b1 + b3 >= -1
    /\
    a0 + b0 + b1 - b2 + b3 >= 0
    /\
    - a0 - a2 + a3 - b1 + p0 >= -2
    /\
    a0 + a1 - b0 - b3 + p0 >= -1
    /\
    a0 - a3 - b0 - b3 + p0 >= -2
    /\
    a0 + a3 + b0 - b3 + p0 >= 0
    /\
    a2 + b0 - b1 - b3 + p0 >= -1
    /\
    - a1 - a2 + b2 - b3 + p0 >= -2
    /\
    a0 + a2 - b1 + b3 + p0 >= 0
    /\
    a1 + a2 + b1 + b2 - p1 >= 0
    /\
    a2 + a3 + b1 + b3 - p1 >= 0
    /\
    a0 + b1 + b2 + b3 - p1 >= 0
    /\
    b0 + b1 + b2 + b3 - p1 >= 0
    /\
    - a2 - a3 - b1 - b3 + p1 >= -3
    /\
    a0 - a1 - a2 + a3 - b2 + p0 >= -2
    /\
    - a1 - a2 + a3 - b0 - b2 + p0 >= -3
    /\
    a1 - a2 - a3 + b0 - b2 + p0 >= -2
    /\
    a1 + a3 - b0 + b1 - b2 + p0 >= -1
    /\
    - a1 - a3 + b1 - b2 - b3 + p0 >= -3
    /\
    - a0 - a2 - a3 - b2 + b3 + p0 >= -3
    /\
    a1 - a2 - b1 + b2 + b3 + p0 >= -1
    /\
    - a0 - a3 + b0 + b1 - b3 - p1 >= -3
;


% We use S-box Analyzer to encode *-LAT: 
% https://github.com//sboxanalyzer
% Simplifying the MILP/SAT constraints ...
% Time used to simplify the constraints: 0.00 seconds
% Number of constraints: 25
% Input:	a0||a1||a2||a3; a0: msb
% Output:	b0||b1||b2||b3; b0: msb

predicate sbox_lin_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                        var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                        var 0..1: d) = 
    - a0 - a1 + a2 + a3 - b0 >= -2
    /\
    - a0 + a1 + a2 + a3 + b0 >= 0
    /\
    a1 + a2 - a3 + b1 + b2 >= 0
    /\
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    a0 + a2 + a3 + b0 - b3 >= 0
    /\
    a2 + a3 + b0 - b1 - b3 >= -1
    /\
    a1 + a2 + a3 + b1 - b3 >= 0
    /\
    - a0 - a2 + a3 - b1 + b3 >= -2
    /\
    a0 + a2 + a3 - b1 + b3 >= 0
    /\
    - a0 + a2 - b0 + b1 + b3 >= -1
    /\
    a0 - a1 + a3 - b2 + b3 >= -1
    /\
    a2 + a3 + b1 - b2 + b3 >= 0
    /\
    a3 - b0 + b1 - b2 + b3 >= -1
    /\
    a0 + b0 + b1 - b2 + b3 >= 0
    /\
    a0 - b0 + b1 + b2 + b3 >= 0
    /\
    - a1 + b0 + b1 + b2 + b3 >= 0
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    a0 + a2 - a3 - b0 + b1 - b3 >= -2
    /\
    - a0 + a2 - a3 + b0 + b1 - b3 >= -2
    /\
    - a1 + a2 - a3 + b1 - b2 - b3 >= -3
    /\
    - a0 - a2 - a3 + b1 - b2 + b3 >= -3
    /\
    a1 - a2 + a3 - b1 + b2 + b3 >= -1
    /\
    a1 - a2 - a3 - b1 - b2 - b3 >= -4
    /\
    - a1 - a2 - a3 - b1 + b2 - b3 >= -4
    /\
    - a2 - a3 - b0 - b1 - b3 >= -4
    /\
    if exists([a0 == 1, a1 == 1, a2 == 1, a3 == 1]) then d = 1
    else d = 0 endif
;

% We use S-box Analyzer to encode deterministic differential behaviour:
% x0: lsb
predicate sbox_diff_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                      var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                      var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y3 == 1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y3 == -1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;

% We use S-box Analyzer to encode deterministic linear behaviour:
% x0: lsb
predicate sbox_lin_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                     var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                     var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y3 == 1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;