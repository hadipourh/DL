/*
MIT License

Copyright (c) 2024 Hosein Hadippour

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

int: RU; % Number of rounds for EU
int: RM; % Number of rounds for EM
int: RL; % Number of rounds for EL
int: RMU; % Number of rounds for EMU
int: RML; % Number of rounds for EML
int: WU; % Weight of active S-boxes in EU
int: WM; % Weight of common active S-boxes in EM
int: WL; % Weight of active S-boxes in EL

constraint assert(RU >= 0, "Invalid value for RU: " ++
        "RU must be greater than or equal to 0");
constraint assert(RM >= 0, "Invalid value for RM: " ++
        "RM must be greater than or equal to 0");
constraint assert(RL >= 0, "Invalid value for RL: " ++
        "RL must be greater than or equal to 0");
constraint assert(RMU >= 0 /\ RMU <= RM, "Invalid value for RMU: " ++
        "RMU must be greater than or equal to 0 and less than or equal to RM");
constraint assert(RML >= 0 /\ RML <= RM, "Invalid value for RML: " ++
        "RML must be greater than or equal to 0 and less than or equal to RM");

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
% X_0 --[SB]--[SR]--> Y_0 --[MC]--> X_1 --[SB]--[SR]--> Y_1 --[MC]--> X_2 --

% X_r: State matrix before the S-box layer in round r
% Y_r: State matrix before the MixColumns layer in round r
% First index: round number
% Second index: row number
% Third index: column number

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  _   _ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____|| | | |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | | | |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |_| |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____| \___/ 
% Constraints for EU                                                                                                

array[0..RU, 0..3, 0..3] of var 0..1: xu;
array[0..<RU, 0..3, 0..3] of var 0..1: yu;

% ShiftRows
constraint forall(r in 0..<RU, row in 0..3, column in 0..3)
(
    yu[r, row, column] = xu[r, row, (column + row) mod 4]
);

constraint forall(r in 0..<RU, column in 0..3)(
    mds_probabilistic(yu[r, 0, column], yu[r, 1, column], yu[r, 2, column], yu[r, 3, column],
                      xu[r + 1, 0, column], xu[r + 1, 1, column], xu[r + 1, 2, column], xu[r + 1, 3, column])
);

% Exclude triavial solutions
constraint sum(row in 0..3, column in 0..3)(xu[0, row, column]) != 0;

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%   ____                    _                _         _           __                _____  __  __ 
%  / ___| ___   _ __   ___ | |_  _ __  __ _ (_) _ __  | |_  ___   / _|  ___   _ __  | ____||  \/  |
% | |    / _ \ | '_ \ / __|| __|| '__|/ _` || || '_ \ | __|/ __| | |_  / _ \ | '__| |  _|  | |\/| |
% | |___| (_) || | | |\__ \| |_ | |  | (_| || || | | || |_ \__ \ |  _|| (_) || |    | |___ | |  | |
%  \____|\___/ |_| |_||___/ \__||_|   \__,_||_||_| |_| \__||___/ |_|   \___/ |_|    |_____||_|  |_|
% Constraints for EM

array[0..RM, 0..3, 0..3] of var 0..1: xmu;
array[0..<RM, 0..3, 0..3] of var 0..1: ymu;

% Connect the end of EU to the beginning of EM
constraint forall(row in 0..3, column in 0..3)(
    xmu[0, row, column] = xu[RU, row, column]
);

% ShiftRows
constraint forall(r in 0..<RM, row in 0..3, column in 0..3)
(
    ymu[r, row, column] = xmu[r, row, (column + row) mod 4]
);

constraint forall(r in 0..<RM, column in 0..3)(
   if r < RMU then
    (
          mds_probabilistic(ymu[r, 0, column], ymu[r, 1, column], ymu[r, 2, column], ymu[r, 3, column],
                            xmu[r + 1, 0, column], xmu[r + 1, 1, column], xmu[r + 1, 2, column], xmu[r + 1, 3, column])
    )
    else
    (
          mds_deterministic(ymu[r, 0, column], ymu[r, 1, column], ymu[r, 2, column], ymu[r, 3, column],
                            xmu[r + 1, 0, column], xmu[r + 1, 1, column], xmu[r + 1, 2, column], xmu[r + 1, 3, column])
    ) endif
);

array[0..RM, 0..3, 0..3] of var 0..1: xml;
array[0..<RM, 0..3, 0..3] of var 0..1: yml;

% ShiftRows
constraint forall(r in 0..<RM, row in 0..3, column in 0..3)
(
    yml[r, row, column] = xml[r, row, (column + row) mod 4]
);

constraint forall(r in 0..<RM, column in 0..3)(
   if r < (RM - RML) then
    (
          mds_deterministic(xml[r + 1, 0, column], xml[r + 1, 1, column], xml[r + 1, 2, column], xml[r + 1, 3, column],
                            yml[r, 0, column], yml[r, 1, column], yml[r, 2, column], yml[r, 3, column])
    )
    else
    (
          mds_probabilistic(xml[r + 1, 0, column], xml[r + 1, 1, column], xml[r + 1, 2, column], xml[r + 1, 3, column],
                            yml[r, 0, column], yml[r, 1, column], yml[r, 2, column], yml[r, 3, column])
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

array[0..RL, 0..3, 0..3] of var 0..1: xl;
array[0..<RL, 0..3, 0..3] of var 0..1: yl;

% Connect the end of EM to the beginning of EL
constraint forall(row in 0..3, column in 0..3)(
    xl[0, row, column] = xml[RM, row, column]
);

% ShiftRows
constraint forall(r in 0..<RL, row in 0..3, column in 0..3)
(
   yl[r, row, column] = xl[r, row, (column + row) mod 4]
);

constraint forall(r in 0..<RL, column in 0..3)(
    mds_probabilistic(yl[r, 0, column], yl[r, 1, column], yl[r, 2, column], yl[r, 3, column],
                      xl[r + 1, 0, column], xl[r + 1, 1, column], xl[r + 1, 2, column], xl[r + 1, 3, column])
);

% Exclude triavial solutions
constraint sum(row in 0..3, column in 0..3)(xl[RL, row, column]) != 0;

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

var int: PU; % Number of active S-boxes in EU
var int: CM; % Number of common active S-boxes in EM
var int: QL; % Number of active S-boxes in EL

constraint PU = sum(r in 0..<RU, row in 0..3, column in 0..3)(xu[r, row, column]);
constraint CM = sum(r in 0..<RM, row in 0..3, column in 0..3)(bool2int(xmu[r, row, column] + xml[r, row, column] == 2));
constraint QL = sum(r in 0..<RL, row in 0..3, column in 0..3)(xl[r, row, column]);

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

% constraint yml[RM - 1, 2, 1] = 1;
% constraint xmu[0, 3, 2] = 1;
solve minimize WU*PU + WM*CM + WL*QL;

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%     _                 _  _  _                       _____                     _    _                    
%    / \   _   _ __  __(_)| |(_)  __ _  _ __  _   _  |  ___|_   _  _ __    ___ | |_ (_)  ___   _ __   ___ 
%   / _ \ | | | |\ \/ /| || || | / _` || '__|| | | | | |_  | | | || '_ \  / __|| __|| | / _ \ | '_ \ / __|
%  / ___ \| |_| | >  < | || || || (_| || |   | |_| | |  _| | |_| || | | || (__ | |_ | || (_) || | | |\__ \
% /_/   \_\\__,_|/_/\_\|_||_||_| \__,_||_|    \__, | |_|    \__,_||_| |_| \___| \__||_| \___/ |_| |_||___/
%                                             |___/       
% Auxiliary functions

predicate xor_cell_wise(var 0..1: a, var 0..1: b, var 0..1: c) = 
    % f(a, b) = a xor b = c
    a + b - c >= 0
    /\
    a - b + c >= 0
    /\
    - a + b + c >= 0
;

predicate mds_probabilistic(var 0..1: x0, var 0..1: x1, var 0..1: x2, var 0..1: x3,
                            var 0..1: y0, var 0..1: y1, var 0..1: y2, var 0..1: y3) = 
    % f(x0, x1, x2, x3) = MDS(y0, y1, y2, y3)
    x0 + x1 + x2 + x3 + y0 + y1 + y2 + y3 in {0, 5, 6, 7, 8}
;

predicate mds_deterministic(var 0..1: x0, var 0..1: x1, var 0..1: x2, var 0..1: x3,
                            var 0..1: y0, var 0..1: y1, var 0..1: y2, var 0..1: y3) = 
    % f(x0, x1, x2, x3) = MDS(y0, y1, y2, y3)
    if (x0 + x1 + x2 + x3 == 0) then
    (
        y0 + y1 + y2 + y3 == 0
    )
    else
    (
        y0 + y1 + y2 + y3 == 4
    ) endif
;
        
