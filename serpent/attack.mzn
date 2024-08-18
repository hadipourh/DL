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

include "sboxes.mzn";

int: NC;
int: RU;
int: RMU;
int: RM;
int: RML;
int: RL;
int: offset;
int: WU;
int: WM;
int: WL;

constraint assert(NC == 32, "Invalid value for NC: " ++
        "NC must be equal to 32");
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
constraint assert(offset in {0, 1, 2, 3, 4, 5, 6, 7}, "Invalid value for offset: " ++
        "offset must be in {0, 1, 2, 3, 4, 5, 6, 7}");
constraint assert(WU >= 1, "Invalid value for WU: " ++
        "WU must be greater than or equal to 1");
constraint assert(WM >= 1, "Invalid value for WM: " ++
        "WM must be greater than or equal to 1");
constraint assert(WL >= 1, "Invalid value for WL: " ++
        "WL must be greater than or equal to 1");

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
% __     __            _         _      _         __  __                       _               
% \ \   / /__ _  _ __ (_)  __ _ | |__  | |  ___  |  \/  |  __ _  _ __   _ __  (_) _ __    __ _ 
%  \ \ / // _` || '__|| | / _` || '_ \ | | / _ \ | |\/| | / _` || '_ \ | '_ \ | || '_ \  / _` |
%   \ V /| (_| || |   | || (_| || |_) || ||  __/ | |  | || (_| || |_) || |_) || || | | || (_| |
%    \_/  \__,_||_|   |_| \__,_||_.__/ |_| \___| |_|  |_| \__,_|| .__/ | .__/ |_||_| |_| \__, |
%                                                               |_|    |_|               |___/ 
% Variable mapping

% State:
% X[0]     Y[0]     X[0]
% X[1]__S__Y[1]__L__X[1]
% X[2]     Y[2]     X[2]
% X[3]     Y[3]     X[3]

% X[r, row, column]
% X[r, row, 0]: the leftmost bit of the column
% X[r, row, 31]: the rightmost bit of the column
% X[round] --S--> Y[round] --S--> X[round + 1]
% T[round, row, column] Array of dummy variables to model the linear layer in differential cryptanalysis

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____|| | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____| \___/ 
% Constraints for EU

array[0..RU, 0..3, 0..<NC] of var 0..1: xu;
array[0..<RU, 0..3, 0..<NC] of var 0..1: yu;
array[0..<RU, 0..3, 0..<NC] of var 0..1: tu;
array[0..<RU, 0..1, 0..<NC] of var 0..1: pu;
array[0..RU, 0..<NC] of var 0..1: du;

% S-box layer
constraint forall(round in 0..<RU, column in 0..<NC)
(
    sbox_diff(xu[round, 0, column], xu[round, 1, column], xu[round, 2, column], xu[round, 3, column],
              yu[round, 0, column], yu[round, 1, column], yu[round, 2, column], yu[round, 3, column],
              pu[round, 0, column], pu[round, 1, column], (offset + round) mod 8)
);

% count the number of active S-boxes 
constraint forall(round in 0..RU, column in 0..<NC)
(
    du[round, column] = bool2int(xu[round, 0, column] + xu[round, 1, column] + xu[round, 2, column] + xu[round, 3, column] >= 1)
);

% Diffusion layer
constraint forall(r in 0..<RU, column in 0..<NC)
(
    xor_3(yu[r, 0, (column + 13) mod NC], yu[r, 1, column], yu[r, 2, (column + 3) mod NC], tu[r, 1, column])
);
constraint forall(r in 0..<RU, column in 0..<NC)
(
    if column > 28 then
    (
        xor_2(yu[r, 2, (column + 3) mod NC], yu[r, 3, column], tu[r, 3, column])
    ) else
    (
        xor_3(yu[r, 0, (column + 13 + 3) mod NC], yu[r, 2, (column + 3) mod NC], yu[r, 3, column], tu[r, 3, column])
    ) endif
);
constraint forall(r in 0..<RU, column in 0..<NC)
(
    xor_3(yu[r, 0, (column + 13) mod NC], tu[r, 1, (column + 1) mod NC], tu[r, 3, (column + 7) mod NC], tu[r, 0, column])
);
constraint forall(r in 0..<RU, column in 0..<NC)
(
    if column > 24 then 
    (
        xor_2(yu[r, 2, (column + 3) mod NC], tu[r, 3, (column + 7) mod NC], tu[r, 2, column])
    ) else
    (
        xor_3(tu[r, 1, (column + 1 + 7) mod NC], yu[r, 2, (column + 3) mod NC], tu[r, 3, (column + 7) mod NC], tu[r, 2, column])
    ) endif
);
constraint forall(r in 0..<RU, column in 0..<NC)
(
    xu[r + 1, 0, column] = tu[r, 0, (column + 5) mod NC]
    /\
    xu[r + 1, 1, column] = tu[r, 1, (column + 1) mod NC]
    /\
    xu[r + 1, 2, column] = tu[r, 2, (column + 22) mod NC]
    /\
    xu[r + 1, 3, column] = tu[r, 3, (column + 7) mod NC]
);

constraint sum(i in 0..3, j in 0..<NC) (xu[0, i, j]) != 0;

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  __  __  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____||  \/  || | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |\/| || | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |  | || |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_|  |_| \___/ 
% Constraints for EMU                                                                                                                                                             

array[0..RM, 0..3, 0..<NC] of var -1..1: xmu;
array[0..<RM, 0..3, 0..<NC] of var -1..1: ymu;
array[0..<RM, 0..3, 0..<NC] of var -1..1: tmu;
array[0..<RM, 0..<NC] of var 0..1: dmu;

% link the input of EM to the output of EU
constraint forall(row in 0..3, column in 0..<NC)(xmu[0, row, column] = xu[RU, row, column]);

% S-box layer
constraint forall(round in 0..<RM, column in 0..<NC)
(
    if round < RMU then 
    (
        sbox_diff_star(xmu[round, 0, column], xmu[round, 1, column], xmu[round, 2, column], xmu[round, 3, column],
                       ymu[round, 0, column], ymu[round, 1, column], ymu[round, 2, column], ymu[round, 3, column],
                       dmu[round, column], (offset + RU + round) mod 8)
    ) else 
    (
        sbox_diff_m(xmu[round, 0, column], xmu[round, 1, column], xmu[round, 2, column], xmu[round, 3, column],
                    ymu[round, 0, column], ymu[round, 1, column], ymu[round, 2, column], ymu[round, 3, column],
                    dmu[round, column], (offset + RU + round) mod 8)
    ) endif
);

% Diffusion layer
constraint forall(r in 0..<RM, column in 0..<NC)
(
    xor_3_m(ymu[r, 0, (column + 13) mod NC], ymu[r, 1, column], ymu[r, 2, (column + 3) mod NC], tmu[r, 1, column])
);
constraint forall(r in 0..<RM, column in 0..<NC)
(
    if column > 28 then
    (
        xor_2_m(ymu[r, 2, (column + 3) mod NC], ymu[r, 3, column], tmu[r, 3, column])
    ) else
    (
        xor_3_m(ymu[r, 0, (column + 13 + 3) mod NC], ymu[r, 2, (column + 3) mod NC], ymu[r, 3, column], tmu[r, 3, column])
    ) endif
);
constraint forall(r in 0..<RM, column in 0..<NC)
(
    xor_3_m(ymu[r, 0, (column + 13) mod NC], tmu[r, 1, (column + 1) mod NC], tmu[r, 3, (column + 7) mod NC], tmu[r, 0, column])
);
constraint forall(r in 0..<RM, column in 0..<NC)
(
    if column > 24 then 
    (
        xor_2_m(ymu[r, 2, (column + 3) mod NC], tmu[r, 3, (column + 7) mod NC], tmu[r, 2, column])
    ) else
    (
        xor_3_m(tmu[r, 1, (column + 1 + 7) mod NC], ymu[r, 2, (column + 3) mod NC], tmu[r, 3, (column + 7) mod NC], tmu[r, 2, column])
    ) endif
);
constraint forall(r in 0..<RM, column in 0..<NC)
(
    xmu[r + 1, 0, column] = tmu[r, 0, (column + 5) mod NC]
    /\
    xmu[r + 1, 1, column] = tmu[r, 1, (column + 1) mod NC]
    /\
    xmu[r + 1, 2, column] = tmu[r, 2, (column + 22) mod NC]
    /\
    xmu[r + 1, 3, column] = tmu[r, 3, (column + 7) mod NC]
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

array[0..RM, 0..3, 0..<NC] of var -1..1: xml;
array[0..<RM, 0..3, 0..<NC] of var -1..1: yml;
array[0..<RM, 0..3, 0..<NC] of var -1..1: tml;
array[0..<RM, 0..<NC] of var 0..1: dml;

% S-box layer
constraint forall(round in 0..<RM, column in 0..<NC)
(
    if round < (RM - RML) then 
    (
        sbox_lin_m(yml[round, 0, column], yml[round, 1, column], yml[round, 2, column], yml[round, 3, column],
                   xml[round, 0, column], xml[round, 1, column], xml[round, 2, column], xml[round, 3, column],
                   dml[round, column], (offset + RU + round) mod 8)
    ) else
    (
        sbox_lin_star(xml[round, 0, column], xml[round, 1, column], xml[round, 2, column], xml[round, 3, column],
                      yml[round, 0, column], yml[round, 1, column], yml[round, 2, column], yml[round, 3, column],                      
                      dml[round, column], (offset + RU + round) mod 8)
    ) endif
);

% Diffusion layer
constraint forall(r in 0..<RM, column in 0..<NC)
(
    if (column < 7) then 
    (
        fork_2_m(xml[r + 1, 1, column], xml[r + 1, 0, (column - 5 + NC) mod NC], tml[r, 1, column])
    ) else 
    (
        fork_3_m(xml[r + 1, 1, column], xml[r + 1, 2, (column - 22 - 7 + NC) mod NC], xml[r + 1, 0, (column - 5 + NC) mod NC], tml[r, 1, column])
    ) endif
);
constraint forall(r in 0..<RM, column in 0..<NC)
(
    fork_3_m(xml[r + 1, 3, column], xml[r + 1, 2, (column - 22 + NC) mod NC], xml[r + 1, 0, (column - 5 + NC) mod NC], tml[r, 3, column])
);
constraint forall(r in 0..<RM, column in 0..<NC)
(
    if (column < 3) then 
    (
        fork_2_m(xml[r + 1, 0, (column - 5 + NC) mod NC], yml[r, 1, column], tml[r, 0, column])
    ) else
    (
        fork_3_m(xml[r + 1, 0, (column - 5 + NC) mod NC], yml[r, 3, (column - 3 + NC) mod NC], yml[r, 1, column], tml[r, 0, column])
    ) endif
);
constraint forall(r in 0..<RM, column in 0..<NC)
(
    fork_3_m(xml[r + 1, 2, (column - 22 + NC) mod NC], yml[r, 3, column], yml[r, 1, column], tml[r, 2, column])
);
constraint forall(r in 0..<RM, column in 0..<NC)
(
    tml[r, 0, column] = yml[r, 0, (column + 13) mod NC]
    /\
    tml[r, 1, column] = yml[r, 1, (column + 1) mod NC]
    /\
    tml[r, 2, column] = yml[r, 2, (column + 3) mod NC]
    /\
    tml[r, 3, column] = yml[r, 3, (column + 7) mod NC]
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

array[0..RL, 0..3, 0..<NC] of var 0..1: xl;
array[0..<RL, 0..3, 0..<NC] of var 0..1: yl;
array[0..<RL, 0..3, 0..<NC] of var 0..1: tl;
array[0..<RL, 0..1, 0..<NC] of var 0..1: pl;
array[0..RL, 0..<NC] of var 0..1: dl;

% link the output of EML to the input of EL
constraint forall(row in 0..3, column in 0..<NC)(xl[0, row, column] = xml[RM, row, column]);

% S-box layer
constraint forall(round in 0..<RL, column in 0..<NC)
(
    sbox_lin(xl[round, 0, column], xl[round, 1, column], xl[round, 2, column], xl[round, 3, column],
             yl[round, 0, column], yl[round, 1, column], yl[round, 2, column], yl[round, 3, column],             
             pl[round, 0, column], pl[round, 1, column], (offset + round + RM + RU) mod 8)
);

% count the number of active S-boxes
constraint forall(round in 0..RL, column in 0..<NC)
(
    dl[round, column] = bool2int(xl[round, 0, column] + xl[round, 1, column] + xl[round, 2, column] + xl[round, 3, column] >= 1)
);

% Diffusion layer
constraint forall(r in 0..<RL, column in 0..<NC)
(
    if (column < 7) then 
    (
        fork_2(xl[r + 1, 1, column], xl[r + 1, 0, (column - 5 + NC) mod NC], tl[r, 1, column])
    ) else 
    (
        fork_3(xl[r + 1, 1, column], xl[r + 1, 2, (column - 22 - 7 + NC) mod NC], xl[r + 1, 0, (column - 5 + NC) mod NC], tl[r, 1, column])
    ) endif
);
constraint forall(r in 0..<RL, column in 0..<NC)
(
    fork_3(xl[r + 1, 3, column], xl[r + 1, 2, (column - 22 + NC) mod NC], xl[r + 1, 0, (column - 5 + NC) mod NC], tl[r, 3, column])
);
constraint forall(r in 0..<RL, column in 0..<NC)
(
    if (column < 3) then 
    (
        fork_2(xl[r + 1, 0, (column - 5 + NC) mod NC], yl[r, 1, column], tl[r, 0, column])
    ) else
    (
        fork_3(xl[r + 1, 0, (column - 5 + NC) mod NC], yl[r, 3, (column - 3 + NC) mod NC], yl[r, 1, column], tl[r, 0, column])
    ) endif
);
constraint forall(r in 0..<RL, column in 0..<NC)
(
    fork_3(xl[r + 1, 2, (column - 22 + NC) mod NC], yl[r, 3, column], yl[r, 1, column], tl[r, 2, column])
);
constraint forall(r in 0..<RL, column in 0..<NC)
(
    tl[r, 0, column] = yl[r, 0, (column + 13) mod NC]
    /\
    tl[r, 1, column] = yl[r, 1, (column + 1) mod NC]
    /\
    tl[r, 2, column] = yl[r, 2, (column + 3) mod NC]
    /\
    tl[r, 3, column] = yl[r, 3, (column + 7) mod NC]
);

constraint sum(row in 0..3, column in 0..<NC)(xl[RL, row, column]) != 0;

% #############################################################################################################################################
% #############################################################################################################################################
var 0..(4*NC + 5000): PU;
var 0..(4*NC + 5000): CMBit;
var 0..(4*NC + 5000): QL;
var 0..(4*NC + 5000): CMSbox;
constraint PU = 3*sum(r in 0..<RU, column in 0..<NC)(pu[r, 0, column]) + 
                2*sum(r in 0..<RU, column in 0..<NC)(pu[r, 1, column]);

constraint CMBit = 
                % 1*sum(r in 0..RM, row in 0..3, column in 0..<NC)(
                % bool2int(xmu[r, row, column] = -1 /\ xml[r, row, column] = -1))
                % +
                % 1*sum(r in 0..RM, row in 0..3, column in 0..<NC)(
                % bool2int(xmu[r, row, column] = -1 /\ xml[r, row, column] = 1))
                % +
                % 1*sum(r in 0..RM, row in 0..3, column in 0..<NC)(
                % bool2int(xmu[r, row, column] = 1 /\ xml[r, row, column] = -1))
                % +
                1*sum(r in 0..<RM, row in 0..3, column in 0..<NC)(
                bool2int(ymu[r, row, column] = -1 /\ yml[r, row, column] = -1))
                +
                1*sum(r in 0..<RM, row in 0..3, column in 0..<NC)(
                bool2int(ymu[r, row, column] = -1 /\ yml[r, row, column] = 1))
                + 
                1*sum(r in 0..<RM, row in 0..3, column in 0..<NC)(
                bool2int(ymu[r, row, column] = 1 /\ yml[r, row, column] = -1));
    
% constraint CMBit = sum(r in 0..<RM, row in 0..3, column in 0..<NC)(
%                 bool2int((ymu[r, row, column] != 0 /\ yml[r, row, column] != 0)));

constraint CMSbox = 2*sum(r in 0..<RM, column in 0..<NC)(bool2int(dmu[r, column] + dml[r, column] == 2));
constraint QL = 4*sum(r in 0..<RL, column in 0..<NC)(pl[r, 0, column]) +               
                2*sum(r in 0..<RL, column in 0..<NC)(pl[r, 1, column]);

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
%     int_search([yl[RL - 1, 0, i] | i in 0..<NC], indomain_min, complete),
%     int_search([yl[RL - 1, 1, i] | i in 0..<NC], indomain_min, complete),
%     int_search([yl[RL - 1, 2, i] | i in 0..<NC], indomain_min, complete),
%     int_search([yl[RL - 1, 3, i] | i in 0..<NC], indomain_min, complete),
%     int_search([yl[RL - 1, 3, i] | i in 0..<NC], indomain_min, complete)]);

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

% constraint forall(row in 2..3, column in 0..<NC)(xu[0, row, column] = 0);
% constraint forall(column in 0..<NC)((xu[0, 0, column] = 1 /\ xu[0, 3, column] = 1) \/ (xu[0, 0, column] = 0 /\ xu[0, 3, column] = 0));
% constraint forall(row in 1..3, column in 0..<NC)(xl[RL, row, column] = 0);

% solve :: restart_constant(10000)
%       :: search_ann minimize (PU + CMSbox + QL);

% #############################################################################################################################################
% constraint forall(row in {0, 1, 3}, column in 0..<NC)(xu[2, row, column] = 0);
% #############################################################################################################################################

% Artificial constraints to reproduce the results of previous works:
% constraint forall(row in 0..3, column in 0..23)(xu[0, row, column] = 0);
% constraint forall(row in 0..3, column in 25..26)(xu[0, row, column] = 0);
% constraint forall(row in 0..3, column in 28..31)(xu[0, row, column] = 0);

% constraint forall(row in 0..3, column in 0..3)(xl[RL, row, column] = 0);
% constraint forall(row in 0..3, column in 5..15)(xl[RL, row, column] = 0);
% constraint forall(row in 0..3, column in 17..20)(xl[RL, row, column] = 0);
% constraint forall(row in 0..3, column in 22..25)(xl[RL, row, column] = 0);
% constraint forall(row in 0..3, column in 27..30)(xl[RL, row, column] = 0);

% constraint sum(i in 0..<NC)(du[0, i]) <= 2;
% constraint sum(i in 0..<NC)(dl[RL, i]) <= 5;

solve minimize (WU*PU + WM*CMBit + WL*QL);
% solve minimize (PU + CMSbox + QL);

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################

% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xu[r, i, j]) = 1 then "1" else show(fix(xu[r, i, j])) endif | j in 0..<NC]) ++ "\"\n"
%     | r in 0..0, i in 0..3];

% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "yu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(yu[r, i, j]) = 1 then "1" else show(fix(yu[r, i, j])) endif | j in 0..<NC]) ++ "\"\n"
%     | r in 0..0 where RU > 0, i in 0..3];

% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xu[r, i, j]) = 1 then "1" else show(fix(xu[r, i, j])) endif | j in 0..<NC]) ++ "\"\n"
%     | r in 1..RU, i in 0..3];

% output ["\n"];
% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xmu_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xmu[r, i, j]) = -1 then "?" else show(fix(xmu[r, i, j])) endif | j in 0..<NC]) ++ "\"\n"
%     | r in 0..RM, i in 0..3];

% output ["\n"];
% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xml_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xml[r, i, j]) = -1 then "?" else show(fix(xml[r, i, j])) endif | j in 0..<NC]) ++ "\"\n"
%     | r in 0..RM, i in 0..3];

% output ["\n"];
% output [
%     if i = 0 then "\n\n" else "" endif ++
%     "xl_" ++ show(r) ++ "_" ++ show(i) ++ ": \"" ++
%     join("", [if fix(xl[r, i, j]) = 1 then "1" else show(fix(xl[r, i, j])) endif | j in 0..<NC]) ++ "\"\n"
%     | r in 0..RL, i in 0..3];

% output ["\n\nCost of EU:                " ++ show(PU)];
% output ["\n\nNo. common active S-boxes: " ++ show(CMSbox)];
% output ["\n\nNo. common active bits:    " ++ show(CMBit)];
% output ["\n\nCost of EL:                " ++ show(QL)];

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

predicate xor_2(var 0..1: a, var 0..1: b, var 0..1: c) = 
    % Generate the constraints for 2-input xor: a + b = c
    - a - b - c >= -2 /\
    a + b - c >= 0 /\
    a - b + c >= 0 /\ 
    - a + b + c >= 0
;

predicate xor_2_m(var -1..1: a, var -1..1: b, var -1..1: c) = 
    % Deterministic version of xor_2
    if (a == -1 \/ b == -1) then c = -1 else xor_2(a, b, c) endif
;

predicate xor_3(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: b) = 
    % Generate the constraints of a 3-input XOR  (b = a0 xor a1 xor a2)
    b - a2 - a1 - a0 >= -2 /\
    - b + a2 - a1 - a0 >= -2 /\
    - b - a2 + a1 - a0 >= -2 /\
    b + a2 + a1 - a0 >= 0 /\
    - b - a2 - a1 + a0 >= -2 /\
    b + a2 - a1 + a0 >= 0 /\
    b - a2 + a1 + a0 >= 0 /\
    - b + a2 + a1 + a0 >= 0
;

predicate xor_3_m(var -1..1: a0, var -1..1: a1, var -1..1: a2, var -1..1: b) = 
    % Deterministic version of xor_3
    if (a0 == -1 \/ a1 == -1 \/ a2 == -1) then b = -1 else xor_3(a0, a1, a2, b) endif
;

predicate fork_2(var 0..1: a, var 0..1: b, var 0..1: c) = 
    % Generate the constraints for 2-way fork: a + b = c
    - a - b - c >= -2 /\
    a + b - c >= 0 /\
    a - b + c >= 0 /\ 
    - a + b + c >= 0
;

predicate fork_2_m(var -1..1: a, var -1..1: b, var -1..1: c) = 
    % Deterministic version of fork_2
    if (a == -1 \/ b == -1) then c = -1 else fork_2(a, b, c) endif
;

predicate fork_3(var 0..1: a0, var 0..1: a1, var 0..1: a2, var 0..1: b) = 
    % Generate the constraints of a 3-input fork  (b = a0 xor a1 xor a2)
    b - a2 - a1 - a0 >= -2 /\
    - b + a2 - a1 - a0 >= -2 /\
    - b - a2 + a1 - a0 >= -2 /\
    b + a2 + a1 - a0 >= 0 /\
    - b - a2 - a1 + a0 >= -2 /\
    b + a2 - a1 + a0 >= 0 /\
    b - a2 + a1 + a0 >= 0 /\
    - b + a2 + a1 + a0 >= 0
;

predicate fork_3_m(var -1..1: a0, var -1..1: a1, var -1..1: a2, var -1..1: b) = 
    % Deterministic version of fork_3
    if (a0 == -1 \/ a1 == -1 \/ a2 == -1) then b = -1 else fork_3(a0, a1, a2, b) endif