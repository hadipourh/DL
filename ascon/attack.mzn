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
% define lookup tables

array[0..4, 0..1] of int: rotation_dictionary = array2d(0..4, 0..1,
[
    19, 28,
    61, 39,
    1, 6,
    10, 17,
    7, 41
]);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____|| | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____| \___/ 
% Constraints for EU

array[0..RU, 0..4, 0..63] of var 0..1: xu;
array[0..(RU - 1), 0..4, 0..63] of var 0..1: yu;
array[0..(RU - 1), 0..2, 0..63] of var 0..1: pu;

% S-box layer
constraint forall(round in 0..(RU - 1), column in 0..63)
(
    sbox_diff(xu[round, 0, column], xu[round, 1, column], xu[round, 2, column], xu[round, 3, column], xu[round, 4, column], 
              yu[round, 0, column], yu[round, 1, column], yu[round, 2, column], yu[round, 3, column], yu[round, 4, column],
              pu[round, 0, column], pu[round, 1, column], pu[round, 2, column])
    % sbox_diff(xu[round, 4, column], xu[round, 3, column], xu[round, 2, column], xu[round, 1, column], xu[round, 0, column],
    %           yu[round, 4, column], yu[round, 3, column], yu[round, 2, column], yu[round, 1, column], yu[round, 0, column],
    %           pu[round, 0, column], pu[round, 1, column], pu[round, 2, column])
);

% Diffusion layer
constraint forall(round in 0..(RU - 1), row in 0..4, column in 0..63)
(
    xor3(yu[round, row, column], yu[round, row, (column - rotation_dictionary[row, 0] + 64) mod 64], yu[round, row, (column - rotation_dictionary[row, 1] + 64) mod 64], 
         xu[round + 1, row, column])
);

constraint sum(i in 0..4, j in 0..63) (xu[0, i, j]) != 0;

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  __  __  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____||  \/  || | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |\/| || | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |  | || |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_|  |_| \___/ 
% Constraints for EMU                                                                                                                                                             

array[0..RM, 0..4, 0..63] of var -1..1: xmu;
array[0..(RM - 1), 0..4, 0..63] of var -1..1: ymu;
array[0..(RM - 1), 0..63] of var 0..1: dmu;

% link the input of EM to the output of EU
constraint forall(row in 0..4, column in 0..63)(xmu[0, row, column] = xu[RU, row, column]);

% S-box layer
constraint forall(round in 0..(RM - 1), column in 0..63)
(
    if round < RMU then 
    (
        sbox_diff_star(xmu[round, 0, column], xmu[round, 1, column], xmu[round, 2, column], xmu[round, 3, column], xmu[round, 4, column], 
                       ymu[round, 0, column], ymu[round, 1, column], ymu[round, 2, column], ymu[round, 3, column], ymu[round, 4, column],
                       dmu[round, column])
        % sbox_diff_star(xmu[round, 4, column], xmu[round, 3, column], xmu[round, 2, column], xmu[round, 1, column], xmu[round, 0, column], 
        %                ymu[round, 4, column], ymu[round, 3, column], ymu[round, 2, column], ymu[round, 1, column], ymu[round, 0, column],
        %                dmu[round, column])
    ) else 
    (
        sbox_diff_m(xmu[round, 0, column], xmu[round, 1, column], xmu[round, 2, column], xmu[round, 3, column], xmu[round, 4, column], 
                    ymu[round, 0, column], ymu[round, 1, column], ymu[round, 2, column], ymu[round, 3, column], ymu[round, 4, column],
                    dmu[round, column])
        % sbox_diff_m(xmu[round, 4, column], xmu[round, 3, column], xmu[round, 2, column], xmu[round, 1, column], xmu[round, 0, column], 
        %             ymu[round, 4, column], ymu[round, 3, column], ymu[round, 2, column], ymu[round, 1, column], ymu[round, 0, column],
        %             dmu[round, column])
    ) endif
);

% Diffusion layer
constraint forall(round in 0..(RM - 1), row in 0..4, column in 0..63)
(
    if round < RMU then 
    (
        xor3(ymu[round, row, column], ymu[round, row, ((column - rotation_dictionary[row, 0] + 64) mod 64)], ymu[round, row, ((column - rotation_dictionary[row, 1] + 64) mod 64)], 
             xmu[round + 1, row, column])
    ) else 
    (
        xor3_pr1(ymu[round, row, column], ymu[round, row, ((column - rotation_dictionary[row, 0] + 64) mod 64)], ymu[round, row, ((column - rotation_dictionary[row, 1] + 64) mod 64)], 
                 xmu[round + 1, row, column])
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

array[0..RM, 0..4, 0..63] of var -1..1: xml;
array[0..(RM - 1), 0..4, 0..63] of var -1..1: yml;
array[0..(RM - 1), 0..63] of var 0..1: dml;

% S-box layer
constraint forall(round in 0..(RM - 1), column in 0..63)
(
    if round < (RM - RML) then 
    (
        sbox_lin_m(yml[round, 0, column], yml[round, 1, column], yml[round, 2, column], yml[round, 3, column], yml[round, 4, column], 
                   xml[round, 0, column], xml[round, 1, column], xml[round, 2, column], xml[round, 3, column], xml[round, 4, column],
                   dml[round, column])
        % sbox_lin_m(yml[round, 4, column], yml[round, 3, column], yml[round, 2, column], yml[round, 1, column], yml[round, 0, column], 
        %            xml[round, 4, column], xml[round, 3, column], xml[round, 2, column], xml[round, 1, column], xml[round, 0, column],
        %            dml[round, column])
    ) else
    (
        sbox_lin_star(yml[round, 0, column], yml[round, 1, column], yml[round, 2, column], yml[round, 3, column], yml[round, 4, column], 
                      xml[round, 0, column], xml[round, 1, column], xml[round, 2, column], xml[round, 3, column], xml[round, 4, column],
                      dml[round, column])
        % sbox_lin_star(yml[round, 4, column], yml[round, 3, column], yml[round, 2, column], yml[round, 1, column], yml[round, 0, column], 
        %               xml[round, 4, column], xml[round, 3, column], xml[round, 2, column], xml[round, 1, column], xml[round, 0, column],
        %               dml[round, column])
    ) endif
);

% Diffusion layer
constraint forall(round in 0..(RM - 1), row in 0..4, column in 0..63)
(
    if round < (RM - RML) then 
    (
        xor3_pr1(xml[round + 1, row, column], xml[round + 1, row, ((column + rotation_dictionary[row, 0] + 64) mod 64)], xml[round + 1, row, ((column + rotation_dictionary[row, 1] + 64) mod 64)], 
                 yml[round, row, column])
    ) else
    (
        xor3(xml[round + 1, row, column], xml[round + 1, row, ((column + rotation_dictionary[row, 0] + 64) mod 64)], xml[round + 1, row, ((column + rotation_dictionary[row, 1] + 64) mod 64)], 
                 yml[round, row, column])
    ) endif

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

array[0..RL, 0..4, 0..63] of var 0..1: xl;
array[0..(RL - 1), 0..4, 0..63] of var 0..1: yl;
array[0..(RL - 1), 0..1, 0..63] of var 0..1: pl;

% link the output of EML to the input of EL
constraint forall(row in 0..4, column in 0..63)(xl[0, row, column] = xml[RM, row, column]);

% S-box layer
constraint forall(round in 0..(RL - 1), column in 0..63)
(
    sbox_lin(yl[round, 0, column], yl[round, 1, column], yl[round, 2, column], yl[round, 3, column], yl[round, 4, column], 
             xl[round, 0, column], xl[round, 1, column], xl[round, 2, column], xl[round, 3, column], xl[round, 4, column],
             pl[round, 0, column], pl[round, 1, column])
    % sbox_lin(yl[round, 4, column], yl[round, 3, column], yl[round, 2, column], yl[round, 1, column], yl[round, 0, column], 
    %          xl[round, 4, column], xl[round, 3, column], xl[round, 2, column], xl[round, 1, column], xl[round, 0, column],
    %          pl[round, 0, column], pl[round, 1, column])
);

% Diffusion layer
constraint forall(round in 0..(RL - 1), row in 0..4, column in 0..63)
(
    xor3(xl[round + 1, row, column], xl[round + 1, row, ((column + rotation_dictionary[row, 0] + 64) mod 64)], xl[round + 1, row, ((column + rotation_dictionary[row, 1] + 64) mod 64)], 
         yl[round, row, column])
);

constraint sum(row in 0..4, column in 0..63)(xl[RL, row, column]) != 0;

% #############################################################################################################################################
% #############################################################################################################################################
var 0..320: PU;
var 0..(320 + 4*320): CM;
var 0..320: CL;
var 0..320: NASM;
constraint PU = 4*sum(r in offset..(RU - 1), column in 0..63)(pu[r, 0, column]) + 
                3*sum(r in offset..(RU - 1), column in 0..63)(pu[r, 1, column]) +
                2*sum(r in offset..(RU - 1), column in 0..63)(pu[r, 2, column]);

constraint CM = 
                1*sum(r in 0..(RM - 1), row in 0..4, column in 0..63)(
                bool2int(xmu[r, row, column] = -1 /\ xml[r, row, column] = -1)
                +
                bool2int(xmu[r, row, column] = -1 /\ xml[r, row, column] = 1)
                +                 
                bool2int(xmu[r, row, column] = 1 /\ xml[r, row, column] = -1));
                % +
                % 1*sum(r in 0..(RM - 1), row in 0..4, column in 0..63)(
                % bool2int(ymu[r, row, column] = -1 /\ yml[r, row, column] = -1))
                % +
                % 1*sum(r in 0..(RM - 1), row in 0..4, column in 0..63)(
                % bool2int(ymu[r, row, column] = -1 /\ yml[r, row, column] = 1))
                % + 
                % 1*sum(r in 0..(RM - 1), row in 0..4, column in 0..63)(
                % bool2int(ymu[r, row, column] = 1 /\ yml[r, row, column] = -1));
    
% constraint CM = sum(r in 0..(RM - 1), row in 0..4, column in 0..63)(
%                 bool2int((ymu[r, row, column] != 0 /\ yml[r, row, column] != 0)));

% constraint CM = sum(r in 0..(RM - 1), row in 0..4, column in 0..63)(
%                 bool2int((xmu[r, row, column] = -1 /\ xml[r, row, column] = -1) 
%                          \/
%                          (xmu[r, row, column] = -1 /\ xml[r, row, column] = 1)
%                          \/
%                          (xmu[r, row, column] = 1 /\ xml[r, row, column] = -1)
%                          \/ 
%                          (xmu[r, row, column] = 1 /\ xml[r, row, column] = 1)));

constraint NASM = sum(r in 0..(RM - 1), column in 0..63)(bool2int(dmu[r, column] + dml[r, column] == 2));
constraint CL = 4*sum(r in 0..(RL - 1), column in 0..63)(pl[r, 0, column]) +               
                2*sum(r in 0..(RL - 1), column in 0..63)(pl[r, 1, column]);

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

% constraint forall(row in 4..4, column in 0..63)(xu[0, row, column] = 0);
% constraint forall(column in 0..63)((xu[0, 0, column] = 1 /\ xu[0, 4, column] = 1) \/ (xu[0, 0, column] = 0 /\ xu[0, 4, column] = 0));
% constraint forall(row in 2..4, column in 0..63)(xl[RL, row, column] = 0);
% constraint PU <= 45;
% constraint CM <= 20;
% constraint CL <= 25;

% solve :: restart_constant(10000)
%       :: search_ann minimize (PU + NASM + CL);

% constraint forall(row in 0..0, column in 0..0)(xu[0, row, column] = 1);
% constraint forall(row in 0..0, column in 1..1)(xu[0, row, column] = 0);
% constraint forall(row in 0..0, column in 2..63)(xu[0, row, column] = 0);
% constraint forall(row in 1..4, column in 0..63)(xu[0, row, column] = 0);

% constraint forall(row in 0..0, column in 0..53)(yml[3, row, column] = 0);
% constraint forall(row in 0..0, column in 54..54)(yml[3, row, column] = 1);
% constraint forall(row in 0..0, column in 55..63)(yml[3, row, column] = 0);
% constraint forall(row in 1..4, column in 0..63)(yml[3, row, column] = 0);

% constraint forall(row in 0..0, column in 0..63)(xu[0, row, column] = 0);
% constraint sum(column in 0..63)(xl[RL, 0, column]) != 0;
% constraint sum(column in 0..63)(xl[RL, 1, column]) != 0;
% constraint forall(row in 1..4, column in 0..63)(xl[RL, row, column] = 0);
solve minimize (PU + CM + CL);
% solve minimize (PU + NASM + CL);

% #############################################################################################################################################
% Constraints to find a DL distinguisher usable in the mode of operation 
% Ascon-AEAD
% ####################################################
% Difference can be in X[3] and X[4]
% Output mask should be zero in X[1], X[2], X[3], X[4]

% constraint forall(row in 0..2, column in 0..63)(xu[0, row, column] = 0);
% constraint forall(row in 1..4, column in 0..63)(xl[RL, row, column] = 0);
% ####################################################

% Ascon-MAC
% ####################################################
% Difference can be in X[0], X[1], X[2], X[3]
% Output linear mask should be zero in X[2], X[3], X[4]

% constraint forall(column in 0..63)(xu[0, 4, column] = 0);
% constraint forall(row in 2..4, column in 0..63)(xl[RL, row, column] = 0);
% #############################################################################################################################################

% #############################################################################################################################################
% Artificial constraints to reproduce the 4-round DL distinguisher of Ascon from the paper: https://eprint.iacr.org/2015/030.pdf
% constraint forall(row in 0..2, column in 0..<64)(xu[0, row, column] = 0);
% constraint forall(row in 3..4, column in 1..<64)(xu[0, row, column] = 0);
% constraint forall(row in 3..4, column in 0..0)(xu[0, row, column] = 1);

% constraint forall(row in 1..4, column in 0..<64)(yml[1, row, column] = 0);
% constraint forall(row in 0..0, column in 0..<53)(yml[1, row, column] = 0);
% constraint forall(row in 0..0, column in 55..<64)(yml[1, row, column] = 0);
% constraint yml[1, 0, 54] = 1;
% #############################################################################################################################################

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################

% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xu[r, i, j]) = 1 then "1" else show(fix(xu[r, i, j])) endif | j in 0..63]) ++ "\"\n"
%     | r in 0..0, i in 0..4];

% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "yu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(yu[r, i, j]) = 1 then "1" else show(fix(yu[r, i, j])) endif | j in 0..63]) ++ "\"\n"
%     | r in 0..0 where RU > 0, i in 0..4];

% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xu[r, i, j]) = 1 then "1" else show(fix(xu[r, i, j])) endif | j in 0..63]) ++ "\"\n"
%     | r in 1..RU, i in 0..4];

% output ["\n"];
% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xmu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xmu[r, i, j]) = -1 then "?" else show(fix(xmu[r, i, j])) endif | j in 0..63]) ++ "\"\n"
%     | r in 0..RM, i in 0..4];

% output ["\n"];
% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xml_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xml[r, i, j]) = -1 then "?" else show(fix(xml[r, i, j])) endif | j in 0..63]) ++ "\"\n"
%     | r in 0..RM, i in 0..4];

% output ["\n"];
% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xl_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xl[r, i, j]) = 1 then "1" else show(fix(xl[r, i, j])) endif | j in 0..63]) ++ "\"\n"
%     | r in 0..RL, i in 0..4];

% output ["\n\nCost of EU: " ++ show(4*sum(r in 0..(RU - 1), column in 0..63)(pu[r, 0, column]) + 
%                                    3*sum(r in 0..(RU - 1), column in 0..63)(pu[r, 1, column]) +
%                                    2*sum(r in 0..(RU - 1), column in 0..63)(pu[r, 2, column]))];
% output ["\n\nCost of EM: " ++ show(sum(r in 0..(RM - 1), column in 0..63)(bool2int(dmu[r, column] + dml[r, column] == 2)))];
% output ["\n\nCost of EM: " ++ show(sum(r in 0..(RM - 1), row in 0..4, column in 0..63)(
%                 bool2int((xmu[r, row, column] = -1 /\ xml[r, row, column] = -1) 
%                          \/
%                          (xmu[r, row, column] = -1 /\ xml[r, row, column] = 1)
%                          \/
%                          (xmu[r, row, column] = 1 /\ xml[r, row, column] = -1))))];
% output ["\n\nCost of EL: " ++ show(sum(r in 0..(RL - 1), column in 0..63)(4*pl[r, 0, column]) + 
%                                    sum(r in 0..(RL - 1), column in 0..63)(2*pl[r, 1, column]))];

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

predicate xor3(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: b) = 
    b - a2 - a1 - a0 >= -2
    /\        
    - b + a2 - a1 - a0 >= -2
    /\
    - b - a2 + a1 - a0 >= -2
    /\
    b + a2 + a1 - a0 >= 0
    /\
    - b - a2 - a1 + a0 >= -2
    /\
    b + a2 - a1 + a0 >= 0
    /\
    b - a2 + a1 + a0 >= 0
    /\
    - b + a2 + a1 + a0 >= 0
;

predicate xor3_pr1(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: y) = 
    if exists([x0 == -1, x1 == -1, x2 == -1]) then y = -1
    else y = (sum([x0, x1, x2]) mod 2) endif
;

predicate sbox_diff(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3, var 0..1: a4, 
                   var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3, var 0..1: b4, 
                   var 0..1: p0, var 0..1: p1, var 0..1: p2) = 
    - p0 - p1 >= -1
    /\
    - p0 - p2 >= -1
    /\
    - p1 - p2 >= -1
    /\
    a0 + a2 - p2 >= 0
    /\
    - a1 + a2 - a3 + p0 >= -1
    /\
    - a0 - a2 + a3 + p1 >= -1
    /\
    - a0 - a2 - a4 + p1 >= -2
    /\
    - a2 - a4 + p0 + p1 >= -1
    /\
    - a3 - b0 + p0 + p1 >= -1
    /\
    a1 - a4 + p1 + p2 >= 0
    /\
    - a1 - a2 + a3 + a4 - b2 >= -2
    /\
    - a1 + a2 + a3 + a4 + b2 >= 0
    /\
    - a0 - a3 - a4 - b2 - b3 >= -4
    /\
    a0 + a3 + a4 - b2 + b3 >= 0
    /\
    - a0 - a3 - a4 + b2 + b3 >= -2
    /\
    a0 + a1 + a3 + a4 - b4 >= 0
    /\
    a0 - a1 + a3 - a4 + p0 >= -1
    /\
    a1 - a3 + a4 + b0 + p0 >= 0
    /\
    a0 - a1 - a4 + b1 + p0 >= -1
    /\
    a1 - a2 + a4 + b2 + p0 >= 0
    /\
    a0 - a1 - a2 - a3 + p1 >= -2
    /\
    - a0 - a1 + a3 - a4 + p1 >= -2
    /\
    a0 + a1 + a2 - b4 + p1 >= 0
    /\
    - a4 + b0 - b4 + p0 + p1 >= -1
    /\
    a0 - a2 + a3 + a4 + p2 >= 0
    /\
    - a0 + a1 + a2 + p0 + p2 >= 0
    /\
    a1 + a3 - b1 + p1 + p2 >= 0
    /\
    - a0 - a1 - a2 - a3 + a4 + b1 >= -3
    /\
    a0 - a1 - a2 + a3 + a4 - b3 >= -2
    /\
    a0 - a1 + a2 - a3 - a4 + b3 >= -2
    /\
    - a2 + a3 - b0 + b1 - b2 - b4 >= -3
    /\
    a1 + a3 - b0 - b1 + b2 - b4 >= -2
    /\
    - a1 + a3 + b0 - b1 + b2 - b4 >= -2
    /\
    - a1 + a3 - b0 + b1 + b2 - b4 >= -2
    /\
    a1 + a3 + b0 + b1 + b2 - b4 >= 0
    /\
    a0 - a1 + a4 + b2 - b3 - b4 >= -2
    /\
    a0 - a1 + a4 - b2 + b3 - b4 >= -2
    /\
    - a1 - a2 + a3 + b0 + b1 + b4 >= -1
    /\
    - a1 + a3 - b0 - b1 + b2 + b4 >= -2
    /\
    a1 + a3 - b0 + b1 + b2 + b4 >= 0
    /\
    a3 - a4 + b0 + b1 + b2 + b4 >= 0
    /\
    a0 - a3 + a4 - b2 - b3 + b4 >= -2
    /\
    a0 + a1 - a2 - a3 - b0 + p0 >= -2
    /\
    - a0 - a1 - a2 - a3 - b1 + p0 >= -4
    /\
    - a0 + a2 + a3 + a4 - b2 + p0 >= -1
    /\
    a0 - a1 - a2 - b1 - b3 + p0 >= -3
    /\
    a1 - a2 - a3 - b0 - b4 + p0 >= -3
    /\
    a0 + a2 - a4 - b0 - b4 + p0 >= -2
    /\
    - a0 - a1 - b2 - b3 - b4 + p0 >= -4
    /\
    - a0 + a1 + a3 - a4 + b4 + p0 >= -1
    /\
    - a1 + a3 - b2 + b3 + b4 + p0 >= -1
    /\
    a0 + a1 - a2 + b0 + b3 - p1 >= -1
    /\
    - a0 + a1 + a3 + a4 + b1 + p1 >= 0
    /\
    a1 + a2 + a4 - b0 - b4 + p1 >= -1
    /\
    a0 + a1 + a4 - b0 + b4 + p2 >= 0
    /\
    a1 - a2 + a4 - b0 + p0 + p2 >= -1
    /\
    a0 - a2 + a3 + b0 + p0 + p2 >= 0
    /\
    a0 - a2 + a4 - b1 + p0 + p2 >= -1
    /\
    a2 - a3 - a4 - b3 + p0 + p2 >= -2
    /\
    a0 + a3 + a4 - b3 + p1 + p2 >= 0
    /\
    a2 + a3 - a4 - b0 - b1 - b2 - b4 >= -4
    /\
    - a1 + a3 - a4 + b2 + b3 - b4 + p0 >= -2
    /\
    - a1 + a3 - a4 + b2 - b3 + b4 - p1 >= -3
    /\
    a2 + a3 + b0 - b1 + b4 + p0 + p1 >= 0
    /\
    - a0 + a1 - a2 - a3 + b0 - b4 + p2 >= -3
    /\
    - a2 + a3 - b0 - b1 - b2 + b4 + p2 >= -3
    /\
    a0 - a1 + a4 + b2 + b3 + b4 + p2 >= 0
    /\
    a1 + a4 + b0 + b4 - p0 + p1 + p2 >= 0
    /\
    - a2 + a3 - a4 + b0 - b1 - b2 - b4 + p2 >= -4
    /\
    a2 + a3 - a4 + b0 + b1 - b2 - b4 + p2 >= -2
    /\
    a2 + a3 - a4 - b0 + b1 - b2 + b4 + p2 >= -2
    /\
    a0 + a1 + a2 + b0 + b4 + p0 - p1 + p2 >= 0
    /\
    - a1 + a2 + a3 - a4 + b0 - b1 - b2 + b4 + p2 >= -3
    /\
    - a0 - a1 + a2 - a3 + a4 + b1 - b3 - b4 - p0 + p1 + p2 >= -5
    /\
    - a0 - a1 + a2 - a3 + a4 - b1 + b3 - b4 - p0 + p1 + p2 >= -5
    /\
    - a0 - a1 + a2 - a3 + a4 - b1 - b3 + b4 - p0 + p1 + p2 >= -5
    /\
    - a0 - a1 - a3 + a4 + b1 + b3 + b4 >= -2
;

predicate sbox_diff_star(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3, var 0..1: a4, 
                        var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3, var 0..1: b4, 
                        var 0..1: d) =
    a1 + a2 + a3 + a4 - b2 >= 0
    /\
    a1 - a2 + a3 + a4 + b2 >= 0
    /\
    - a1 + a2 + a3 + a4 + b2 >= 0
    /\
    a0 + a3 + a4 + b2 - b3 >= 0
    /\
    a0 + a3 + a4 - b2 + b3 >= 0
    /\
    - a0 - a3 - a4 + b2 + b3 >= -2
    /\
    a1 - a3 - a4 - b0 - b4 >= -3
    /\
    a1 + a3 + a4 - b0 - b4 >= -1
    /\
    a1 + a3 - a4 + b0 + b4 >= 0
    /\
    a1 - a3 + a4 + b0 + b4 >= 0
    /\
    - a0 - a1 - a2 - a3 - a4 - b1 >= -5
    /\
    a0 - a1 - a2 - a3 - a4 + b1 >= -3
    /\
    - a0 + a1 + a3 + b0 + b1 + b2 >= 0
    /\
    a0 + a1 + a2 - a3 + b0 - b3 >= -1
    /\
    - a1 - a2 - a3 - a4 - b1 - b3 >= -5
    /\
    a0 + a1 - a2 - a3 - a4 + b3 >= -2
    /\
    - a0 + a1 + a2 + a4 - b0 - b4 >= -2
    /\
    - a0 + a1 - a2 - a3 + b0 - b4 >= -3
    /\
    a0 + a1 - a2 + a3 + b0 - b4 >= -1
    /\
    - a0 + a1 + a2 - a4 + b0 - b4 >= -2
    /\
    a0 + a1 + a2 + a4 + b0 - b4 >= 0
    /\
    a0 + a1 + a2 + a3 + b1 - b4 >= 0
    /\
    a1 - a4 - b0 - b1 + b2 - b4 >= -3
    /\
    - a1 + a3 - b0 + b1 + b2 - b4 >= -2
    /\
    a0 - a1 + a4 - b2 + b3 - b4 >= -2
    /\
    - a0 - a1 - a4 + b2 + b3 - b4 >= -3
    /\
    - a0 + a1 - a2 + a3 - b0 + b4 >= -2
    /\
    - a0 + a1 + a2 - a4 - b0 + b4 >= -2
    /\
    - a1 - a2 + a3 + b0 + b1 + b4 >= -1
    /\
    a1 + a3 + b0 - b1 + b2 + b4 >= 0
    /\
    a1 + a3 - b0 + b1 + b2 + b4 >= 0
    /\
    a3 - a4 + b0 + b1 + b2 + b4 >= 0
    /\
    a0 + a1 - a3 - b0 - b3 + b4 >= -2
    /\
    a0 + a1 + a4 - b0 + b3 + b4 >= 0
    /\
    a0 + a1 - a4 + b0 + b3 + b4 >= 0
    /\
    - a0 + a3 - a4 - b2 + b3 + b4 >= -2
    /\
    - a2 + a3 - a4 + b0 - b1 - b2 - b4 >= -4
    /\
    - a2 + a3 - a4 - b0 + b1 - b2 - b4 >= -4
    /\
    - a0 - a1 + a3 - a4 - b2 - b3 - b4 >= -5
    /\
    - a2 + a3 - a4 - b0 - b1 - b2 + b4 >= -4
    /\
    - a0 - a1 + a3 - a4 + b2 - b3 + b4 >= -3
    /\
    a0 - a1 - a3 + a4 + b2 + b3 + b4 >= -1
    /\
    - a0 - a1 + a2 - a3 + a4 + b1 - b3 - b4 >= -4
    /\
    - a0 - a1 + a2 - a3 + a4 - b1 + b3 - b4 >= -4
    /\
    - a0 - a1 + a2 - a3 + a4 - b1 - b3 + b4 >= -4
    /\
    a2 + a3 - a4 - b0 - b1 - b2 - b4 >= -4
    /\
    a2 + a3 - a4 + b0 + b1 - b2 - b4 >= -2
    /\
    a2 + a3 - a4 + b0 - b1 - b2 + b4 >= -2
    /\
    a2 + a3 - a4 - b0 + b1 - b2 + b4 >= -2
    /\
    a0 - a1 + a4 + b2 - b3 - b4 >= -2
    /\
    - a0 - a1 - a3 + a4 + b1 + b3 + b4 >= -2
    /\
    a0 - a3 + a4 - b2 - b3 + b4 >= -2
    /\
    a0 - a1 + a2 - a3 - a4 + b3 >= -2
    /\
    - a1 + a3 + b0 - b1 + b2 - b4 >= -2
    /\
    - a1 + a3 - b0 - b1 + b2 + b4 >= -2
    /\
    - a0 - a3 - a4 - b2 - b3 >= -4
    /\
    - a0 - a1 - a2 - a3 + a4 + b1 >= -3
    /\
    - a1 - a2 + a3 + a4 - b2 >= -2
    /\
    a0 - a1 - a2 - a3 + a4 - b1 >= -3
    /\
    if exists([a0 == 1, a1 == 1, a2 == 1, a3 == 1, a4 == 1]) then d = 1
    else d = 0 endif
;

predicate sbox_lin(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3, var 0..1: a4, 
                   var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3, var 0..1: b4, 
                   var 0..1: p0, var 0..1: p1) = 
    - a0 - a3 - p1 >= -2
    /\
    - a1 - a3 - p1 >= -2
    /\
    - a2 + p0 + p1 >= 0
    /\
    - b2 + p0 + p1 >= 0
    /\
    - b3 + p0 + p1 >= 0
    /\
    - b4 + p0 + p1 >= 0
    /\
    a1 + a2 + a3 - p0 >= 0
    /\
    a0 + a1 + a4 - p0 >= 0
    /\
    a0 + a1 + b2 - p0 >= 0
    /\
    a2 - a3 + b2 + p0 >= 0
    /\
    - a0 + a4 + b3 + p0 >= 0
    /\
    a0 - a1 - a4 - p1 >= -2
    /\
    - a0 - a1 + a4 - p1 >= -2
    /\
    - a2 + a3 - b0 - p1 >= -2
    /\
    a0 - a4 - b2 - p1 >= -2
    /\
    a0 - a1 + a3 + a4 + b0 >= 0
    /\
    a0 + a1 - a2 - a3 - b2 >= -2
    /\
    - a0 - a1 - a4 - b1 - b2 >= -4
    /\
    a1 - a2 + a3 + b0 + b2 >= 0
    /\
    - a0 - a1 - a4 + b1 + b2 >= -2
    /\
    - a0 + a3 - a4 + b0 - p0 >= -2
    /\
    a0 + a3 - a4 + b3 + p0 >= 0
    /\
    - a0 + a1 + b0 - b2 - p1 >= -2
    /\
    a0 + a1 + b1 - b2 - p1 >= -1
    /\
    - a1 + a2 - a4 - b4 - p1 >= -3
    /\
    - a1 + a3 - a4 - b0 - b2 - b4 >= -4
    /\
    a1 + a3 + a4 - b0 - b2 - b4 >= -2
    /\
    a2 - b0 + b1 + b2 - b3 - b4 >= -2
    /\
    - a0 + a2 + a3 - a4 - b0 + b4 >= -2
    /\
    a0 + a2 + a3 + a4 - b0 + b4 >= 0
    /\
    - a1 + b0 + b1 + b2 + b3 + b4 >= 0
    /\
    - a1 + a2 + a3 - b0 - b4 - p0 >= -3
    /\
    a1 + a2 + b2 + b3 + b4 - p0 >= 0
    /\
    a1 + a3 - b0 - b3 - b4 + p0 >= -2
    /\
    a0 + a1 + a2 + a3 + a4 - p1 >= 0
    /\
    a1 + a2 + a4 - b0 + b2 - p1 >= -1
    /\
    - a0 - a4 + b0 + b2 - b4 - p1 >= -3
    /\
    a1 + a2 + b2 - b3 - b4 + p1 >= -1
    /\
    - a1 + a2 + a3 + b0 + b4 + p1 >= 0
    /\
    a3 + b0 - b1 + b2 + b4 + p1 >= 0
    /\
    - a1 - a2 - a3 - a4 - b1 - b2 + b3 >= -5
    /\
    a0 - a2 - a3 - a4 + b1 + b2 + b3 >= -2
    /\
    a0 + a2 + a4 + b1 - b2 + b3 - b4 >= -1
    /\
    a0 + a4 + b0 + b1 - b2 + b3 - b4 >= -1
    /\
    a0 + a2 - b0 - b1 + b2 + b3 - b4 >= -2
    /\
    - a0 - a1 - a2 + a3 - a4 + b1 + b4 >= -3
    /\
    a1 + a2 - a3 + a4 - b2 - b3 + b4 >= -2
    /\
    a0 + a2 + a4 + b1 - b2 - b3 + b4 >= -1
    /\
    a0 + a4 + b0 + b1 - b2 - b3 + b4 >= -1
    /\
    a1 + a2 + a3 + b0 + b4 + p0 - p1 >= 0
    /\
    - a1 - a2 + a3 + b0 - b2 - b4 + p1 >= -3
    /\
    - a1 - a2 + a3 - a4 + b2 + b4 + p1 >= -2
    /\
    a0 - a1 + a2 - a3 - a4 + b0 - b1 - b3 >= -4
    /\
    a0 - a1 + a2 - a3 - a4 - b0 + b1 - b3 >= -4
    /\
    - a0 - a1 + a2 - a3 + a4 + b0 - b1 - b4 >= -4
    /\
    - a0 - a1 - a3 + a4 + b0 - b1 + b2 - b4 >= -4
    /\
    a0 - a1 + a2 + a4 - b1 - b2 - b3 - b4 >= -4
    /\
    a0 - a3 + a4 + b0 - b1 + b2 + b3 - b4 >= -2
    /\
    - a0 - a1 + a2 - a3 - b0 - b1 - b2 + b4 >= -5
    /\
    - a0 - a2 + a4 - b0 + b1 - b2 - b3 + b4 >= -4
    /\
    a0 + a2 + a4 + b1 + b2 + b3 + b4 - p0 >= 0
    /\
    a0 + a4 + b0 - b1 - b2 - b3 - b4 + p1 >= -3
    /\
    a0 + a2 + a4 - b1 - b2 + b3 + b4 + p1 >= -1
    /\
    a0 + a4 + b0 - b1 - b2 + b3 + b4 + p1 >= -1
    /\
    - a1 + a2 - a3 + a4 - b0 - b1 + b2 - b3 + b4 >= -4
    /\
    a0 - a2 - a3 + a4 - b0 - b1 + b2 + b3 + b4 >= -3
    /\
    - a0 - a2 + a3 + a4 - b0 - b2 + b4 - p0 + p1 >= -4
    /\
    a0 - a1 - a2 - a3 - a4 + b1 - b2 - b3 - p0 + p1 >= -6
    /\
    a0 - a1 - a2 - a3 - a4 - b1 + b2 - b3 - p0 + p1 >= -6
    /\
    a0 - a1 + a2 - a3 - a4 - b0 - b1 + b3 - p0 + p1 >= -5
    /\
    a0 - a1 + a2 - a3 - a4 + b0 + b1 + b3 - p0 + p1 >= -3
    /\
    - a0 - a1 - a3 + a4 - b0 - b1 + b2 + b4 - p0 + p1 >= -5
    /\
    - a0 - a1 - a2 - a3 + a4 + b0 + b1 - b2 - b4 - p0 + p1 >= -6
    /\
    - a1 - a2 - a3 + a4 - b0 + b1 - b2 + b3 + b4 - p0 + p1 >= -5
    /\
    a0 - a2 - a3 + a4 - b0 - b1 - b2 - b3 + b4 >= -5
    /\
    a0 - a2 - a3 + a4 - b0 + b1 + b2 - b3 + b4 >= -3
    /\
    - a0 - a1 - a2 - a3 - b0 - b1 - b2 - b4 >= -7
    /\
    - a0 - a1 - a2 - a3 + b0 - b1 - b2 + b4 >= -5
    /\
    - a0 - a1 + a2 + a4 - b0 + b1 - b4 >= -3
    /\
    a1 - a2 - a3 + b0 + b3 + b4 >= -1
    /\
    - a2 - a3 - b0 + b1 + b2 + b3 - b4 >= -3
    /\
    - a0 - a1 + b0 + b1 + b2 + b4 >= -1
    /\
    a0 - a2 - a3 - b0 + b1 - b2 - b3 - b4 >= -5
    /\
    a0 - a2 - a3 - b0 - b1 + b2 - b3 - b4 >= -5
    /\
    a0 + a4 + b0 + b1 + b2 - b3 - b4 >= -1
    /\
    a1 - a2 - a3 - b0 - b3 + b4 >= -3
    /\
    a1 - a2 - a3 + b0 - b3 - b4 >= -3
    /\
    - a2 - a3 - b0 - b1 - b2 + b3 - b4 >= -5
    /\
    - a0 + a1 + a2 - a4 + b3 + b4 >= -1
    /\
    - a0 - a1 + a2 + a4 + b0 + b1 + b4 >= -1
    /\
    - a0 + a1 + a2 - a4 - b3 - b4 >= -3
    /\
    a0 + b0 - b1 + b2 - b3 + b4 >= -1
    /\
    a1 - a2 - a3 - b0 + b3 - b4 >= -3
    /\
    - a0 - a1 + a3 - b0 + b2 - b4 >= -3
    /\
    - a0 - a1 - b0 + b1 + b2 - b4 >= -3
    /\
    a1 + a2 + a4 - b2 + b3 - b4 >= -1
;

predicate sbox_lin_star(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3, var 0..1: a4, 
                        var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3, var 0..1: b4, 
                        var 0..1: d) = 
    - a0 - a1 + a3 - a4 - b0 >= -3
    /\
    a0 + a1 - a2 - a3 - b2 >= -2
    /\
    - a0 - a1 - a4 - b1 - b2 >= -4
    /\
    a0 + a1 + a4 + b1 - b2 >= 0
    /\
    a0 + a1 + a2 - a3 + b2 >= 0
    /\
    a0 + a1 - a2 + a3 + b2 >= 0
    /\
    a1 - a2 + a3 + b0 + b2 >= 0
    /\
    - a0 - a1 - a4 + b1 + b2 >= -2
    /\
    a1 - a2 + a3 + a4 - b0 - b2 >= -2
    /\
    a0 + a1 + a3 + a4 + b2 - b3 >= 0
    /\
    a0 + a2 + a3 - a4 - b0 - b4 >= -2
    /\
    - a0 + a2 + a3 + a4 - b0 - b4 >= -2
    /\
    - a0 + a2 + a3 - a4 + b0 - b4 >= -2
    /\
    - a0 + a3 - a4 + b0 + b2 - b4 >= -2
    /\
    a1 + a2 + a3 - b2 - b3 - b4 >= -2
    /\
    a1 + a2 + a3 + a4 + b3 - b4 >= 0
    /\
    a1 + a2 + a3 + b0 + b3 - b4 >= 0
    /\
    a0 + a2 + a3 + a4 - b0 + b4 >= 0
    /\
    a0 + a2 + a3 - a4 + b0 + b4 >= 0
    /\
    a1 + a2 + a3 - a4 - b2 + b4 >= -1
    /\
    a2 + a3 + a4 + b0 - b2 + b4 >= 0
    /\
    - a0 + a1 + a2 + a3 + b2 + b4 >= 0
    /\
    a0 - a2 + a3 - a4 + b2 + b4 >= -1
    /\
    a3 + a4 + b0 - b1 + b2 + b4 >= 0
    /\
    a1 + a2 + a3 - b0 + b3 + b4 >= 0
    /\
    - a0 + a1 + a2 + b2 + b3 + b4 >= 0
    /\
    - a3 + b0 + b1 + b2 + b3 + b4 >= 0
    /\
    - a1 - a2 - a3 - a4 - b1 - b2 + b3 >= -5
    /\
    a0 - a2 - a3 - a4 + b1 + b2 + b3 >= -2
    /\
    a0 - a1 - a2 + a3 - a4 - b2 - b4 >= -4
    /\
    - a1 - a2 + a3 + b0 - b1 - b2 - b4 >= -4
    /\
    a0 + a2 + a4 + b1 + b2 - b3 - b4 >= -1
    /\
    a2 + a4 - b0 + b1 - b2 + b3 - b4 >= -2
    /\
    a0 + a2 + a4 - b1 + b2 + b3 - b4 >= -1
    /\
    a0 + a4 + b0 - b1 + b2 + b3 - b4 >= -1
    /\
    - a0 - a2 + a3 + a4 - b0 - b2 + b4 >= -3
    /\
    - a0 - a2 + a3 - a4 + b0 - b2 + b4 >= -3
    /\
    a1 + a2 - a3 + a4 - b2 - b3 + b4 >= -2
    /\
    a0 + a2 + a4 + b1 - b2 - b3 + b4 >= -1
    /\
    a0 + a4 + b0 + b1 - b2 - b3 + b4 >= -1
    /\
    a0 - a1 + a2 - a3 - a4 + b0 - b1 - b3 >= -4
    /\
    a0 - a1 + a2 - a3 - a4 - b0 + b1 - b3 >= -4
    /\
    a0 - a1 - a2 - a3 - a4 + b1 - b2 - b3 >= -5
    /\
    a0 - a1 - a2 - a3 - a4 - b1 + b2 - b3 >= -5
    /\
    a0 - a1 + a2 - a3 - a4 - b0 - b1 + b3 >= -4
    /\
    a0 - a1 + a2 - a3 - a4 + b0 + b1 + b3 >= -2
    /\
    - a0 - a1 + a2 - a3 + b0 - b1 - b2 - b4 >= -5
    /\
    - a0 - a1 - a2 + a4 + b0 + b1 - b2 - b4 >= -4
    /\
    - a0 - a1 - a3 + a4 + b0 - b1 + b2 - b4 >= -4
    /\
    - a0 + a1 + a2 - a3 + a4 + b2 - b3 - b4 >= -3
    /\
    - a0 - a1 + a2 - a3 - b0 - b1 - b2 + b4 >= -5
    /\
    - a0 - a1 - a3 + a4 - b0 - b1 + b2 + b4 >= -4
    /\
    - a0 - a2 + a4 - b0 + b1 - b2 - b3 + b4 >= -4
    /\
    - a1 + a2 - a3 - b0 - b1 - b2 + b3 + b4 >= -4
    /\
    a0 + a2 + a4 - b0 + b1 + b2 + b3 + b4 >= 0
    /\
    a0 - a2 - a3 + a4 + b0 - b1 - b2 - b3 - b4 >= -5
    /\
    a0 - a1 - a3 + a4 + b0 + b1 - b2 + b3 - b4 >= -3
    /\
    - a1 + a2 - a3 + a4 - b0 - b1 + b2 - b3 + b4 >= -4
    /\
    a0 - a1 - a3 + a4 + b0 - b1 - b2 + b3 + b4 >= -3
    /\
    - a1 - a2 - a3 + a4 - b0 + b1 - b2 + b3 + b4 >= -4
    /\
    a0 - a2 - a3 + a4 - b0 - b1 + b2 + b3 + b4 >= -3
    /\
    a0 - a2 - a3 + a4 - b0 + b1 + b2 - b3 + b4 >= -3
    /\
    a0 - a2 - a3 + a4 - b0 - b1 - b2 - b3 + b4 >= -5
    /\
    a0 - a1 + a2 + a4 - b1 - b2 - b3 - b4 >= -4
    /\
    - a0 - a1 - a2 - a3 - b0 - b1 - b2 - b4 >= -7
    /\
    - a0 - a1 - a2 - a3 + b0 - b1 - b2 + b4 >= -5
    /\
    - a0 - a1 + b0 + b1 + b2 + b4 >= -1
    /\
    - a0 - a1 + a2 + a4 - b0 + b1 - b4 >= -3
    /\
    a1 - a2 - a3 + b0 + b3 + b4 >= -1
    /\
    - a2 - a3 - b0 + b1 + b2 + b3 - b4 >= -3
    /\
    a0 - a2 - a3 - b0 + b1 - b2 - b3 - b4 >= -5
    /\
    a0 - a2 - a3 - b0 - b1 + b2 - b3 - b4 >= -5
    /\
    a0 + a4 + b0 + b1 + b2 - b3 - b4 >= -1
    /\
    a1 - a2 - a3 - b0 - b3 + b4 >= -3
    /\
    - a2 - a3 - b0 - b1 - b2 + b3 - b4 >= -5
    /\
    - a0 + a1 + a2 - a4 + b3 + b4 >= -1
    /\
    a1 - a2 - a3 + b0 - b3 - b4 >= -3
    /\
    - a0 - a1 + a2 + a4 + b0 + b1 + b4 >= -1
    /\
    - a0 + a1 + a2 - a4 - b3 - b4 >= -3
    /\
    a0 + b0 - b1 + b2 - b3 + b4 >= -1
    /\
    a1 - a2 - a3 - b0 + b3 - b4 >= -3
    /\
    - a0 - a1 + a3 - b0 + b2 - b4 >= -3
    /\
    - a0 - a1 - b0 + b1 + b2 - b4 >= -3
    /\
    - a0 + a1 + a3 - a4 + b0 >= -1
    /\
    a1 + a2 + a4 - b2 + b3 - b4 >= -1
    /\
    a0 - a1 + a3 + a4 + b0 >= 0
    /\
    if exists([a0 == 1, a1 == 1, a2 == 1, a3 == 1, a4 == 1]) then d = 1
    else d = 0 endif
;

predicate sbox_diff_m(var -1..1: x4, var -1..1: x3, var -1..1: x2, var -1..1: x1, var -1..1: x0, 
                      var -1..1: y4, var -1..1: y3, var -1..1: y2, var -1..1: y1, var -1..1: y0,
                      var 0..1: d) = 
    if (x4 == 0 /\ x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y4 == 0 /\ y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y4 == -1 /\ y3 == 1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 0 /\ x1 == 1 /\ x0 == 0) then (y4 == 1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 0 /\ x1 == 1 /\ x0 == 1) then (y4 == -1 /\ y3 == -1 /\ y2 == -1 /\ y1 == 0 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == 1 /\ y1 == 1 /\ y0 == 0)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 1) then (y4 == 1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 1 /\ x1 == 1 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 1 /\ x1 == 1 /\ x0 == 1) then (y4 == 0 /\ y3 == -1 /\ y2 == -1 /\ y1 == 1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == -1 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 0)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == -1 /\ x1 == 1 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x4 == 0 /\ x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == 1 /\ y1 == 1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 1 /\ x2 == 0 /\ x1 == 1 /\ x0 == 1) then (y4 == -1 /\ y3 == -1 /\ y2 == -1 /\ y1 == 1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 1 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == 0 /\ y1 == 0 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 1 /\ x2 == 1 /\ x1 == 1 /\ x0 == 0) then (y4 == -1 /\ y3 == 0 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 1 /\ x2 == 1 /\ x1 == 1 /\ x0 == 1) then (y4 == -1 /\ y3 == 1 /\ y2 == -1 /\ y1 == 0 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == 1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y4 == 1 /\ y3 == 0 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == 0 /\ x1 == 1 /\ x0 == 1) then (y4 == 0 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 0)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y4 == 0 /\ y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 1) then (y4 == -1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == 1 /\ x1 == 1 /\ x0 == 0) then (y4 == 1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == 1 /\ x1 == 1 /\ x0 == 1) then (y4 == -1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 0)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == -1 /\ x1 == 0 /\ x0 == 1) then (y4 == -1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == -1 /\ x1 == 1 /\ x0 == 1) then (y4 == -1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 0)
    elseif (x4 == 1 /\ x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 1 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 1 /\ x2 == 1 /\ x1 == 1 /\ x0 == 0) then (y4 == -1 /\ y3 == 1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 1 /\ x2 == 1 /\ x1 == 1 /\ x0 == 1) then (y4 == -1 /\ y3 == 0 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == -1 /\ x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == -1 /\ x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == -1 /\ x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == -1 /\ x3 == 1 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    else (y4 = -1 /\ y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0 /\ x4 == 0) then d = 0 else d = 1 endif
;

predicate sbox_lin_m(var -1..1: x4, var -1..1: x3, var -1..1: x2, var -1..1: x1, var -1..1: x0, 
                     var -1..1: y4, var -1..1: y3, var -1..1: y2, var -1..1: y1, var -1..1: y0,
                     var 0..1: d) = 
    if (x4 == 0 /\ x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y4 == 0 /\ y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y4 == -1 /\ y3 == -1 /\ y2 == 0 /\ y1 == 1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == -1) then (y4 == -1 /\ y3 == -1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 0 /\ x1 == 1 /\ x0 == 0) then (y4 == -1 /\ y3 == 1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 0 /\ x1 == 1 /\ x0 == 1) then (y4 == -1 /\ y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 0 /\ x1 == 1 /\ x0 == -1) then (y4 == -1 /\ y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y4 == 0 /\ y3 == 1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 1) then (y4 == -1 /\ y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == -1) then (y4 == -1 /\ y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 1 /\ x1 == 1 /\ x0 == 0) then (y4 == -1 /\ y3 == 0 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 1 /\ x1 == 1 /\ x0 == 1) then (y4 == -1 /\ y3 == -1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == 1 /\ x1 == 1 /\ x0 == -1) then (y4 == -1 /\ y3 == -1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 0 /\ x2 == -1 /\ x1 == 0 /\ x0 == 0) then (y4 == 0 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y4 == 1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x4 == 0 /\ x3 == 1 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y4 == 1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 0 /\ x3 == 1 /\ x2 == -1 /\ x1 == 0 /\ x0 == 0) then (y4 == 1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y4 == -1 /\ y3 == -1 /\ y2 == -1 /\ y1 == 1 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y4 == 1 /\ y3 == -1 /\ y2 == -1 /\ y1 == 0 /\ y0 == 1)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 1) then (y4 == 1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 0 /\ x2 == -1 /\ x1 == 0 /\ x0 == 1) then (y4 == 1 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y4 == 0 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 0)
    elseif (x4 == 1 /\ x3 == 1 /\ x2 == 1 /\ x1 == 0 /\ x0 == 1) then (y4 == 0 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x4 == 1 /\ x3 == 1 /\ x2 == -1 /\ x1 == 0 /\ x0 == 1) then (y4 == 0 /\ y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    else (y4 = -1 /\ y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0 /\ x4 == 0) then d = 0 else d = 1 endif
;