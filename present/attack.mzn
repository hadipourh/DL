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

int: RU;
int: RMU;
int: RM;
int: RML;
int: RL;
int: WU;
int: WM;
int: WL;
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
constraint assert(WU >= 1, "Invalid value for WU: " ++
        "WU must be greater than or equal to 1");
constraint assert(WM >= 1, "Invalid value for WM: " ++
        "WM must be greater than or equal to 1");
constraint assert(WL >= 1, "Invalid value for WL: " ++
        "WL must be greater than or equal to 1");
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

array[0..63] of int: round_permutation = array1d(0..63, [0, 16, 32, 48, 1, 17, 33, 49, 2, 18, 34, 50, 3, 19, 35, 51, 4, 20, 36, 52, 5, 21, 37, 53, 6, 22, 38, 54, 7, 23, 39, 55, 8, 24, 40, 56, 9, 25, 41, 57, 10, 26, 42, 58, 11, 27, 43, 59, 12, 28, 44, 60, 13, 29, 45, 61, 14, 30, 46, 62, 15, 31, 47, 63]);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____|| | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____| \___/ 
% Constraints for EU

array[0..RU, 0..63] of var 0..1: xu;
array[0..(RU - 1), 0..63] of var 0..1: yu;
array[0..(RU - 1), 0..15, 0..1] of var 0..1: pu;

% S-box layer
constraint forall(round in 0..(RU - 1), i in 0..15)
(
    sbox_diff_probabilistic(xu[round, 4*i], xu[round, 4*i + 1], xu[round, 4*i + 2], xu[round, 4*i + 3], 
                            yu[round, 4*i], yu[round, 4*i + 1], yu[round, 4*i + 2], yu[round, 4*i + 3],
                            pu[round, i, 0], pu[round, i, 1])
);

% Diffusion layer
constraint forall(round in 0..(RU - 1), i in 0..63)
(
    xu[round + 1, round_permutation[i]] = yu[round, i]
);

constraint sum(i in 0..63) (xu[0, i]) != 0;



% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  __  __ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____||  \/  |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |\/| |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |  | |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_|  |_|
% Constraints for EM

array[0..RM, 0..63] of var -1..1: xmu;
array[0..(RM - 1), 0..63] of var -1..1: ymu;
array[0..(RM - 1), 0..15] of var 0..1: dmu;
array[0..(RM - 1), 0..15, 0..1] of var 0..1: pmu;

% Connect the end of EU to the beginning of EM
constraint forall(i in 0..63)
(
    xmu[0, i] = xu[RU, i]
);

% S-box layer
constraint forall(round in 0..(RM - 1), i in 0..15)
(
    if round < RMU then
    (
        % sbox_diff_star_ddt(xmu[round, 4*i], xmu[round, 4*i + 1], xmu[round, 4*i + 2], xmu[round, 4*i + 3], 
        %                    ymu[round, 4*i], ymu[round, 4*i + 1], ymu[round, 4*i + 2], ymu[round, 4*i + 3])
        sbox_diff_probabilistic(xmu[round, 4*i], xmu[round, 4*i + 1], xmu[round, 4*i + 2], xmu[round, 4*i + 3], 
                                ymu[round, 4*i], ymu[round, 4*i + 1], ymu[round, 4*i + 2], ymu[round, 4*i + 3],
                                pmu[round, i, 0], pmu[round, i, 1])
    ) else
    (
        sbox_diff_deterministic(xmu[round, 4*i], xmu[round, 4*i + 1], xmu[round, 4*i + 2], xmu[round, 4*i + 3], 
                                ymu[round, 4*i], ymu[round, 4*i + 1], ymu[round, 4*i + 2], ymu[round, 4*i + 3])
    ) endif
);

% Diffusion layer
constraint forall(round in 0..(RM - 1), i in 0..63)
(
    ymu[round, i] = xmu[round + 1, round_permutation[i]]
);

array[0..RM, 0..63] of var -1..1: xml;
array[0..(RM - 1), 0..63] of var -1..1: yml;
array[0..(RM - 1), 0..15] of var 0..1: dml;
array[0..(RM - 1), 0..15, 0..1] of var 0..1: pml;

% S-box layer
constraint forall(round in 0..(RM - 1), i in 0..15)
(
    if round < (RM - RML) then
    (
        sbox_inv_lin_deterministic(yml[round, 4*i], yml[round, 4*i + 1], yml[round, 4*i + 2], yml[round, 4*i + 3], 
                                   xml[round, 4*i], xml[round, 4*i + 1], xml[round, 4*i + 2], xml[round, 4*i + 3])
    ) else
    (
        % sbox_inv_lin_star_lat(yml[round, 4*i], yml[round, 4*i + 1], yml[round, 4*i + 2], yml[round, 4*i + 3], 
        %                       xml[round, 4*i], xml[round, 4*i + 1], xml[round, 4*i + 2], xml[round, 4*i + 3])
        sbox_inv_lin_probabilistic(yml[round, 4*i], yml[round, 4*i + 1], yml[round, 4*i + 2], yml[round, 4*i + 3], 
                                   xml[round, 4*i], xml[round, 4*i + 1], xml[round, 4*i + 2], xml[round, 4*i + 3],
                                   pml[round, i, 0], pml[round, i, 1])
    ) endif
);

% Diffusion layer
constraint forall(round in 0..(RM - 1), i in 0..63)
(
    yml[round, i] = xml[round + 1, round_permutation[i]]
);

% Encode the active S-boxes (in the level of word)
constraint forall(r in 0..(RM - 1), i in 0..15)
(
    if (ymu[r, 4*i] != 0 \/ ymu[r, 4*i + 1] != 0 \/ ymu[r, 4*i + 2] != 0 \/ ymu[r, 4*i + 3] != 0) then
    (
        dmu[r, i] = 1   
    ) else
    (
        dmu[r, i] = 0
    ) endif
);

constraint forall(r in 0..(RM - 1), i in 0..15)
(
    if (yml[r, 4*i] != 0 \/ yml[r, 4*i + 1] != 0 \/ yml[r, 4*i + 2] != 0 \/ yml[r, 4*i + 3] != 0) then
    (
        dml[r, i] = 1   
    ) else
    (
        dml[r, i] = 0
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

array[0..RL, 0..63] of var 0..1: xl;
array[0..(RL - 1), 0..63] of var 0..1: yl;
array[0..(RL - 1), 0..15, 0..1] of var 0..1: pl;

% Connect the end of EM to the beginning of EL
constraint forall(i in 0..63)
(
    xml[RM, i] = xl[0, i]
);

% S-box layer
constraint forall(round in 0..(RL - 1), i in 0..15)
(
    sbox_inv_lin_probabilistic(yl[round, 4*i], yl[round, 4*i + 1], yl[round, 4*i + 2], yl[round, 4*i + 3], 
                               xl[round, 4*i], xl[round, 4*i + 1], xl[round, 4*i + 2], xl[round, 4*i + 3],
                               pl[round, i, 0], pl[round, i, 1])
);

% Diffusion layer
constraint forall(round in 0..(RL - 1), i in 0..63)
(
    yl[round, i] = xl[round + 1, round_permutation[i]]
);

constraint sum(i in 0..63)(xl[RL, i]) != 0;

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
% ann: search_ann;
% search_ann = seq_search([
%     int_search([yl[RL - 1, i] | i in 0..63], indomain_median, complete)]);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ___   _      _              _    _               _____                     _    _               
%  / _ \ | |__  (_)  ___   ___ | |_ (_)__   __ ___  |  ___|_   _  _ __    ___ | |_ (_)  ___   _ __  
% | | | || '_ \ | | / _ \ / __|| __|| |\ \ / // _ \ | |_  | | | || '_ \  / __|| __|| | / _ \ | '_ \ 
% | |_| || |_) || ||  __/| (__ | |_ | | \ V /|  __/ |  _| | |_| || | | || (__ | |_ | || (_) || | | |
%  \___/ |_.__/_/ | \___| \___| \__||_|  \_/  \___| |_|    \__,_||_| |_| \___| \__||_| \___/ |_| |_|
%             |__/                                                                                  
% some artificial constraints to control the input and output masks
var 0..64: PU;
var 0..64: QL;
var 0..(64*RM + 1000): CMB;
var 0..(64 + 1000): CMW;
var 0..64: PMU;
var 0..64: QML;

constraint PU = sum(r in 0..(RU - 1), i in 0..15)(3*pu[r, i, 0] + 2*pu[r, i, 1]);
constraint PMU = sum(r in 0..(RM - 1), i in 0..15)(3*pmu[r, i, 0] + 2*pmu[r, i, 1]);
constraint QL = sum(r in 0..(RL - 1), i in 0..15)(4*pl[r, i, 0] + 2*pl[r, i, 1]);
constraint QML = sum(r in 0..(RM - 1), i in 0..15)(4*pml[r, i, 0] + 2*pml[r, i, 1]);
constraint CMB = 1*sum(r in 0..(RM - 1), i in 0..63)(
                bool2int(ymu[r, i] = -1 /\ yml[r, i] = -1))
                +
                1*sum(r in 0..(RM - 1), i in 0..63)(
                bool2int(ymu[r, i] = -1 /\ yml[r, i] = 1))
                + 
                1*sum(r in 0..(RM - 1), i in 0..63)(
                bool2int(ymu[r, i] = 1 /\ yml[r, i] = -1));
                % +
                % 1*sum(r in 0..(RM - 1), i in 0..63)(
                % bool2int(xmu[r, i] = -1 /\ xml[r, i] = -1))
                % +
                % 1*sum(r in 0..(RM - 1), i in 0..63)(
                % bool2int(xmu[r, i] = -1 /\ xml[r, i] = 1))
                % + 
                % 1*sum(r in 0..(RM - 1), i in 0..63)(
                % bool2int(xmu[r, i] = 1 /\ xml[r, i] = -1));
constraint CMW = sum(r in 0..(RM - 1), i in 0..15)(bool2int(dmu[r, i] + dml[r, i] = 2));

solve minimize WU*PU + WM*PMU + WM*CMW + WM*QML + WL*QL;
% solve :: restart_constant(10000)
%       :: search_ann minimize PU + CMB + QL;

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

predicate sbox_diff_probabilistic(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3,
                                  var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3,
                                  var 0..1: p0, var 0..1: p1) = 
    % Simplifying the MILP/SAT constraints ...
    % Time used to simplify the constraints: 0.03 seconds
    % Number of constraints: 54
    % Input:	a0||a1||a2||a3; a0: msb
    % Output:	b0||b1||b2||b3; b0: msb
    % Weight: 3.0000 p0 + 2.0000 p1
    - p0 - p1 >= -1
    /\
    a3 + b3 - p1 >= 0
    /\
    - a2 + p0 + p1 >= 0
    /\
    a3 + b0 + b2 - p0 >= 0
    /\
    a3 - b0 + b2 + p0 >= 0
    /\
    - a0 - a2 - a3 + b0 - b2 >= -3
    /\
    - a0 - a1 - a3 - b0 + b2 >= -3
    /\
    - a1 - a2 + a3 + b0 + b2 >= -1
    /\
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    a1 + a2 + b0 + b2 - b3 >= 0
    /\
    a0 + a1 + a2 - b1 + b3 >= 0
    /\
    - a3 + b0 + b1 + b2 + b3 >= 0
    /\
    a1 + a2 + a3 + b3 - p0 >= 0
    /\
    - a0 - a1 + a3 - b2 + p0 >= -2
    /\
    - a1 + a2 - b0 - b2 + p0 >= -2
    /\
    a0 - a1 + b0 - b2 + p0 >= -1
    /\
    - a0 + a1 + b0 - b2 + p0 >= -1
    /\
    - a3 - b0 + b1 - b2 + p0 >= -2
    /\
    - a1 + a2 - a3 - b3 + p0 >= -2
    /\
    a1 - a2 + b1 - b3 - p1 >= -2
    /\
    a0 - a2 - b0 + b3 - p1 >= -2
    /\
    a0 - a2 + b2 + b3 - p1 >= -1
    /\
    - a1 - a2 - b0 - b2 + p1 >= -3
    /\
    a1 + a2 - b0 - b2 + p1 >= -1
    /\
    - a1 - a2 + b0 + b2 + p1 >= -1
    /\
    a0 - a1 - a2 - b3 + p1 >= -2
    /\
    a1 + a2 - a3 - b3 + p1 >= -1
    /\
    a3 - b0 - b2 - b3 + p1 >= -2
    /\
    - b0 - b1 - b2 - b3 + p1 >= -3
    /\
    b0 - b1 + b2 - b3 + p1 >= -1
    /\
    - a1 + a2 + a3 + b0 + b1 - b3 >= -1
    /\
    a1 - a2 + a3 - b1 - b2 - b3 >= -3
    /\
    a1 - a3 - b0 - b1 - b2 - b3 >= -4
    /\
    a0 - a2 - a3 - b0 + b2 - b3 >= -3
    /\
    a1 - a2 + a3 + b1 + b2 - b3 >= -1
    /\
    a0 + a1 + b0 + b1 - b2 + b3 >= 0
    /\
    a0 + a2 - b0 + b1 + b2 + b3 >= 0
    /\
    a0 + a2 + b0 + b1 + b3 - p0 >= 0
    /\
    - a1 + a2 + b0 + b2 + b3 + p0 >= 0
    /\
    - a0 + a1 - a2 - a3 - b1 - p1 >= -4
    /\
    - a2 - a3 + b0 - b1 - b3 - p1 >= -4
    /\
    a0 + a2 - a3 + b0 - b2 + p1 >= -1
    /\
    - a0 + a1 + a2 + b0 + b2 + p1 >= 0
    /\
    a0 + a1 - b0 + b1 + b2 + p1 >= 0
    /\
    - a0 - a1 + a2 + b0 - b1 - b2 + b3 >= -3
    /\
    - a0 + a1 - a2 - b0 - b1 + b2 + b3 >= -3
    /\
    - a0 + a1 + a2 - a3 - b0 + b2 + p0 >= -2
    /\
    - a0 + a1 - a2 + a3 + b0 - b1 + p1 >= -2
    /\
    - a0 - a1 + a2 + a3 - b1 + b2 + p1 >= -2
    /\
    - a1 + a2 + a3 - b1 + b2 - b3 + p1 >= -2
    /\
    - a1 + a2 - a3 - b0 - b1 + b3 + p1 >= -3
    /\
    - a1 + a2 - a3 + b0 - b2 + b3 + p1 >= -2
    /\
    a1 - a2 - a3 - b1 - b2 + b3 + p1 >= -3
    /\
    a1 - a2 - a3 - b0 + b2 + b3 + p1 >= -2
;

predicate sbox_diff_deterministic(var -1..1: x3, var -1..1: x2, var -1..1: x1, var -1..1: x0, var -1..1: y3, var -1..1: y2, var -1..1: y1, var -1..1: y0) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 0)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
;

predicate sbox_diff_star_ddt(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3,
                             var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3) = 
    % Simplifying the MILP/SAT constraints ...
    % Time used to simplify the constraints: 0.00 seconds
    % Number of constraints: 38
    % Input:	a0||a1||a2||a3; a0: msb
    % Output:	b0||b1||b2||b3; b0: msb
    a0 + a1 + a2 + a3 - b0 >= 0
    /\
    - a0 - a2 - a3 + b0 - b2 >= -3
    /\
    - a0 - a1 - a3 - b0 + b2 >= -3
    /\
    - a1 - a2 + a3 + b0 + b2 >= -1
    /\
    - a0 + a1 + a2 - a3 - b3 >= -2
    /\
    a1 + a2 + b0 + b2 - b3 >= 0
    /\
    - a3 + b0 - b1 + b2 - b3 >= -2
    /\
    a3 + b0 + b1 + b2 - b3 >= 0
    /\
    - a0 + a1 + a2 + a3 + b3 >= 0
    /\
    a0 + a1 + a2 - b1 + b3 >= 0
    /\
    - a1 + a3 + b0 + b2 + b3 >= 0
    /\
    - a2 + a3 + b0 + b2 + b3 >= 0
    /\
    - a3 + b0 + b1 + b2 + b3 >= 0
    /\
    a0 - a1 + a2 - a3 + b0 - b2 >= -2
    /\
    - a1 + a2 - a3 - b0 - b1 - b2 >= -4
    /\
    - a1 - a2 - a3 - b0 + b1 - b2 >= -4
    /\
    a1 + a2 - a3 - b0 + b1 - b2 >= -2
    /\
    a0 + a2 + a3 + b0 + b1 - b2 >= 0
    /\
    - a1 + a2 + a3 - b0 - b1 - b3 >= -3
    /\
    a1 - a2 + a3 - b0 + b1 - b3 >= -2
    /\
    - a0 - a2 + a3 - b0 - b2 - b3 >= -4
    /\
    a0 - a1 - a2 + b0 - b2 - b3 >= -3
    /\
    a0 + a1 + a3 - b1 - b2 - b3 >= -2
    /\
    a1 - a3 - b0 - b1 - b2 - b3 >= -4
    /\
    - a1 + a2 + a3 + b1 - b2 - b3 >= -2
    /\
    a0 - a1 - a2 - b0 + b2 - b3 >= -3
    /\
    a0 - a2 - a3 - b0 + b2 - b3 >= -3
    /\
    a1 - a2 - a3 - b0 - b1 + b3 >= -3
    /\
    - a1 - a2 + a3 - b0 - b2 + b3 >= -3
    /\
    a0 - a3 - b0 - b1 - b2 + b3 >= -3
    /\
    a0 + a1 + b0 + b1 - b2 + b3 >= 0
    /\
    a0 - a1 - a2 + b0 + b2 + b3 >= -1
    /\
    - a0 + a2 + a3 - b1 + b2 + b3 >= -1
    /\
    a0 + a1 - b0 + b1 + b2 + b3 >= 0
    /\
    a0 + a2 - b0 + b1 + b2 + b3 >= 0
    /\
    - a0 + a1 - a2 + a3 + b0 - b1 - b2 >= -3
    /\
    - a0 - a1 + a2 + b0 - b1 - b2 + b3 >= -3
    /\
    - a0 + a1 - a2 - b0 - b1 + b2 + b3 >= -3
;

predicate sbox_inv_lin_probabilistic(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3,
                                     var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3,
                                     var 0..1: p0, var 0..1: p1) = 
    % Simplifying the MILP/SAT constraints ...
    % Time used to simplify the constraints: 0.01 seconds
    % Number of constraints: 49
    % Input:	a0||a1||a2||a3; a0: msb
    % Output:	b0||b1||b2||b3; b0: msb
    % Weight: 4.0000 p0 + 2.0000 p1
    - p0 - p1 >= -1
    /\
    a0 + a1 + a2 - p0 >= 0
    /\
    b0 + b1 + b2 - p0 >= 0
    /\
    - a1 - b0 - b3 + p0 >= -2
    /\
    - a0 + a1 - a2 + p1 >= -1
    /\
    b0 - b1 - b2 + p1 >= -1
    /\
    a0 + a1 + a2 - a3 + b0 >= 0
    /\
    - a0 + a1 + b0 + b1 + b2 >= 0
    /\
    - a0 - a2 - b1 - b2 + b3 >= -3
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 + a2 - a3 - b0 + p0 >= -2
    /\
    - a0 - a1 - a3 - b1 + p0 >= -3
    /\
    - a1 + a2 + a3 + b1 + p0 >= 0
    /\
    a0 - a1 + b1 + b3 + p0 >= 0
    /\
    - a1 + a2 + b2 + b3 + p0 >= 0
    /\
    a0 - a1 - a2 + b0 + p1 >= -1
    /\
    - a0 - a1 + a2 + b0 + p1 >= -1
    /\
    a1 - b0 + b1 - b2 + p1 >= -1
    /\
    a1 - b0 - b1 + b2 + p1 >= -1
    /\
    - a0 - a2 + b1 + b2 + p1 >= -1
    /\
    b0 + b1 + b2 - b3 + p1 >= 0
    /\
    - a0 - a1 - a2 - b0 - b1 - b2 >= -5
    /\
    - a0 + a1 + a2 + b0 + b1 + p0 >= 0
    /\
    a0 - a1 - a3 + b0 - b2 + p0 >= -2
    /\
    - a0 - a1 - a2 + b1 - b2 + p0 >= -3
    /\
    a0 - a2 - a3 + b1 - b2 + p0 >= -2
    /\
    a0 + a1 - a2 - a3 + b2 + p0 >= -1
    /\
    a2 - a3 + b0 - b1 + b2 + p0 >= -1
    /\
    a0 - a1 - a2 - b1 - b3 + p0 >= -3
    /\
    a2 + a3 - b0 - b1 - b3 + p0 >= -2
    /\
    a0 + a1 + a3 + b1 - b3 + p0 >= 0
    /\
    - a0 - a1 + a2 - b2 - b3 + p0 >= -3
    /\
    a0 - a2 - b0 - b2 - b3 + p0 >= -3
    /\
    - a2 - b0 + b1 - b2 - b3 + p0 >= -3
    /\
    - a0 + a1 + a2 + b2 - b3 + p0 >= -1
    /\
    - a1 - a2 + a3 + b2 - b3 + p0 >= -2
    /\
    a0 + b0 - b1 + b2 - b3 + p0 >= -1
    /\
    a0 - a2 - b0 - b1 + b3 + p0 >= -2
    /\
    - a2 + a3 + b0 + b1 + b3 + p0 >= 0
    /\
    - a0 + a3 - b0 - b2 + b3 + p0 >= -2
    /\
    a0 + a2 + b0 - b2 + b3 + p0 >= 0
    /\
    a3 + b0 - b1 + b2 + b3 + p0 >= 0
    /\
    a1 - b0 + b1 + b2 + b3 + p0 >= 0
    /\
    a0 + a1 + a2 + b3 + p0 - p1 >= 0
    /\
    a0 + a2 - b0 - b1 - b2 + p1 >= -2
    /\
    a0 + a2 - b0 + b1 + b2 + p1 >= 0
    /\
    - a0 - a2 + a3 - b0 - b1 + b2 + p0 >= -3
    /\
    a1 + a3 + b0 - b1 - b2 + p0 - p1 >= -2
    /\
    - a0 + a1 - a2 - a3 - b3 + p0 - p1 >= -4    
;

predicate sbox_inv_lin_deterministic(var -1..1: x3, var -1..1: x2, var -1..1: x1, var -1..1: x0, var -1..1: y3, var -1..1: y2, var -1..1: y1, var -1..1: y0) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y3 == 1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 1 /\ x0 == 0) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 1 /\ x0 == 1) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 0)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
;

predicate sbox_inv_lin_star_lat(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: a3,
                                var 0..1: b0, var 0..1: b1, var 0..1: b2, var 0..1: b3) = 
    % Simplifying the MILP/SAT constraints ...
    % Time used to simplify the constraints: 0.00 seconds
    % Number of constraints: 36
    % Input:	a0||a1||a2||a3; a0: msb
    % Output:	b0||b1||b2||b3; b0: msb
    a0 + a1 + a2 - a3 + b0 >= 0
    /\
    a0 + a1 + a2 + a3 - b1 >= 0
    /\
    a0 + a1 + a2 + a3 - b2 >= 0
    /\
    a1 - a3 + b0 + b1 + b2 >= 0
    /\
    - a0 + a3 + b0 + b1 + b2 >= 0
    /\
    - a2 + a3 + b0 + b1 + b2 >= 0
    /\
    a0 + a1 + a2 - b0 + b3 >= 0
    /\
    - a1 + b0 + b1 + b2 + b3 >= 0
    /\
    - a1 - a3 + b0 - b1 - b2 - b3 >= -4
    /\
    a1 - a2 - b0 + b1 - b2 - b3 >= -3
    /\
    - a0 + a1 - b0 - b1 + b2 - b3 >= -3
    /\
    a0 + a2 + a3 + b1 + b2 - b3 >= 0
    /\
    a1 + a2 - b0 + b1 - b2 + b3 >= -1
    /\
    - a0 + a1 - a2 + b1 + b2 + b3 >= -1
    /\
    - a0 + a1 + a2 - a3 - b0 + b1 - b2 >= -3
    /\
    a0 - a1 + a2 - b0 - b1 - b2 - b3 >= -4
    /\
    a1 - a2 - a3 - b0 - b1 + b2 - b3 >= -4
    /\
    - a0 - a1 - a2 - b0 + b1 + b2 - b3 >= -4
    /\
    a0 + a1 - a3 - b0 + b1 - b2 + b3 >= -2
    /\
    a0 + a1 - a2 - b0 - b1 + b2 + b3 >= -2
    /\
    a1 + a2 - a3 - b0 - b1 + b2 + b3 >= -2
    /\
    a0 + a2 + b0 - b1 - b2 + b3 >= -1
    /\
    a0 - a1 - a2 + b0 + b1 + b3 >= -1
    /\
    - a0 - a1 - a2 - b0 - b1 - b2 >= -5
    /\
    a0 - a1 - a2 + b0 - b1 - b3 >= -3
    /\
    - a0 - a1 + a2 + b0 - b2 - b3 >= -3
    /\
    - a0 - a1 + a2 + b0 + b2 + b3 >= -1
    /\
    - a0 - a1 + a2 - a3 + b0 - b1 >= -3
    /\
    a0 - a1 - a2 - a3 + b0 - b2 >= -3
    /\
    a0 - a1 + a2 - b0 + b1 + b2 >= -1
    /\
    a0 - a1 - a2 + a3 + b0 + b2 >= -1
    /\
    - a0 - a1 + a2 + a3 + b0 + b1 >= -1
    /\
    - a0 - a2 - b1 - b2 + b3 >= -3
    /\
    - a0 + a1 - a2 - a3 - b3 >= -3
    /\
    - a0 + a1 - a2 + a3 + b3 >= -1
    /\
    a1 + a3 + b0 - b1 - b2 >= -1
;