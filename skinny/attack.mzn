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

%########## Distinguisher parameters ############
int: RU;
int: RM;
int: RL;
int: RMU;
int: RML;
int: WU;
int: WM;
int: WL;
int: is_related_tweakey;
%################################################
int: cell_size;
int: block_size;
block_size = 16*cell_size;
int: domain_margine = 100;
int: pr_tolerance = 15;
int: Rd;     % the length of Ed
int: NPT;    % number of parallel tweakey paths in the tweakey schedule
Rd = RU + RM + RL;

constraint assert(cell_size in {4, 8}, "cell_size must be 4 or 8");
constraint assert(NPT in {1, 2, 3, 4}, "NPT must be in {1, 2, 3, 4}");
constraint assert(RU >= 0 /\ RM >= 0 /\ RL >= 0, "RU, RM, and RL must be non-negative");
constraint assert(RMU in 0..RM, "Invalid value for RMU: " ++
                "RMU must be in {0, ..., RM}");
constraint assert(RML in 0..RM, "Invalid value for RML: " ++
                "RML must be in {0, ..., RM}");
constraint assert(is_related_tweakey in {0, 1}, "is_related_tweakey must be 0 or 1");
%##################################################################################
%##################################################################################
%##################################################################################
%  _                   _                    _____       _      _            
% | |     ___    ___  | | __ _   _  _ __   |_   _|__ _ | |__  | |  ___  ___ 
% | |    / _ \  / _ \ | |/ /| | | || '_ \    | | / _` || '_ \ | | / _ \/ __|
% | |___| (_) || (_) ||   < | |_| || |_) |   | || (_| || |_) || ||  __/\__ \
% |_____|\___/  \___/ |_|\_\ \__,_|| .__/    |_| \__,_||_.__/ |_| \___||___/
%                                  |_|                                      


% lower bound for the number of active S-boxes in SKINNY-TK1, SKINNY-TK2, SKINNY-TK3, and SKINNY-TK4
% Reference: https://ia.cr/2020/542
array[0..4, 1..22] of int: lower_bound = array2d(0..4, 1..22, 
                                          [|1, 2, 5, 8, 12, 16, 26, 36, 41, 46, 51, 55, 58, 61, 66, 75, 82, 88, 92, 96, 102, 108|
                                            0, 0, 1, 2, 3 , 6 , 10, 13, 16, 23, 32, 38, 41, 45, 49, 54, 59, 62, 66, 70, 75 , 79 |
                                            0, 0, 0, 0, 1 , 2 , 3 , 6 , 9 , 12, 16, 21, 25, 31, 35, 40, 43, 47, 52, 57, 59 , 64 |
                                            0, 0, 0, 0, 0 , 0 , 1 , 2 , 3 , 6 , 10, 13, 16, 19, 24, 27, 31, 35, 43, 45, 48 , 51 |
                                            0, 0, 0, 0, 0 , 0 , 0 , 0 , 1 , 2 , 3 , 6 , 9 , 12, 16, 19, 21, 24, 30, 35, 39 , 41 |]);

array[0..15] of int: round_permutation = array1d(0..15, [0, 1, 2, 3, 7, 4, 5, 6, 10, 11, 8, 9, 13, 14, 15, 12]);

array[0..15] of int: inv_tweakey_permutation = array1d(0..15, [8, 9, 10, 11, 12, 13, 14, 15, 2, 0, 4, 7, 6, 3, 5, 1]);

%##################################################################################
%##################################################################################
%##################################################################################
%  _____                        _                  ____         _                _         _       
% |_   _|__      __ ___   __ _ | | __ ___  _   _  / ___|   ___ | |__    ___   __| | _   _ | |  ___ 
%   | |  \ \ /\ / // _ \ / _` || |/ // _ \| | | | \___ \  / __|| '_ \  / _ \ / _` || | | || | / _ \
%   | |   \ V  V /|  __/| (_| ||   <|  __/| |_| |  ___) || (__ | | | ||  __/| (_| || |_| || ||  __/
%   |_|    \_/\_/  \___| \__,_||_|\_\\___| \__, | |____/  \___||_| |_| \___| \__,_| \__,_||_| \___|
%                                          |___/                                                   
% Model the difference propagation through the tweakey schedule

array[0..Rd, 0..15] of var int: inv_tkp;
constraint forall(i in 0..15) (inv_tkp[0, i] = i);
constraint forall(r in 1..(Rd - 1), i in 0..15)
(
    inv_tkp[r, i] = inv_tweakey_permutation[inv_tkp[r - 1, i]]
);

array[0..15] of var 0..1: LANEU; % activeness pattern in each lane through the upper trail
array[0..(RU + RM - 1), 0..15] of var 0..1: DSTKU; % activeness pattern in each round tweakey through the upper trail
array[0..(RM + RL - 1), 0..15] of var 0..1: DSTKL; % activeness pattern in each round tweakey through the lower trail
constraint forall(i in 0..15)
(
    forall(r in 0..(RU + RM - 1)) (DSTKU[r, inv_tkp[r, i]] - LANEU[i] <= 0)
    /\
    sum(r in 0..(RU + RM - 1))(DSTKU[r, inv_tkp[r, i]]) - (RU + RM)*LANEU[i] >= (1 - NPT) 
);

constraint if is_related_tweakey == 0 then(
    sum(i in 0..15)(LANEU[i]) = 0
) else(
    true
) endif;
%##################################################################################
%##################################################################################
%##################################################################################
%  ____   _       _    _                       _       _                 
% |  _ \ (_) ___ | |_ (_) _ __    __ _  _   _ (_) ___ | |__    ___  _ __ 
% | | | || |/ __|| __|| || '_ \  / _` || | | || |/ __|| '_ \  / _ \| '__|
% | |_| || |\__ \| |_ | || | | || (_| || |_| || |\__ \| | | ||  __/| |   
% |____/ |_||___/ \__||_||_| |_| \__, | \__,_||_||___/|_| |_| \___||_|   
%                                |___/                                   
% Model the distinguisher
              
array[0..(RU + RM), 0..15] of var 0..1: DXU;
array[0..(RU + RM - 1), 0..15] of var 0..1: DZU;

%##################################################################################
%##################################################################################
%  _____  _   _ 
% | ____|| | | |
% |  _|  | | | |
% | |___ | |_| |
% |_____| \___/ 

constraint forall(r in 0..(RU - 1), i in 0..7) 
(
    xor_probabilistic(DXU[r, i], DSTKU[r, i], DZU[r, i]) 
    /\ 
    DZU[r, i + 8] = DXU[r, i + 8]
);
constraint forall(r in 0..(RU - 1), i in 0..3)
(
    mixcol_diff_probabilistic(DZU[r, round_permutation[i]], DZU[r, round_permutation[i + 4]], DZU[r, round_permutation[i + 8]], DZU[r, round_permutation[i + 12]],
                            DXU[r + 1, i], DXU[r + 1, i + 4], DXU[r + 1, i + 8], DXU[r + 1, i + 12])
);

%##################################################################################
%##################################################################################
%  _____  __  __  _   _ 
% | ____||  \/  || | | |
% |  _|  | |\/| || | | |
% | |___ | |  | || |_| |
% |_____||_|  |_| \___/ 
                      
constraint forall(r in RU..(RU + RM - 1), i in 0..7) 
(
    if r < RU + RMU then
    (
        xor_probabilistic(DXU[r, i], DSTKU[r, i], DZU[r, i])
    ) else
    (
        xor_deterministic(DXU[r, i], DSTKU[r, i], DZU[r, i])
    ) endif
    /\ 
    DZU[r, i + 8] = DXU[r, i + 8]
);
constraint forall(r in RU..(RU + RM - 1), i in 0..3)
(
    if r < RU + RMU then
    (
        mixcol_diff_probabilistic(DZU[r, round_permutation[i]], DZU[r, round_permutation[i + 4]], DZU[r, round_permutation[i + 8]], DZU[r, round_permutation[i + 12]],
                                DXU[r + 1, i], DXU[r + 1, i + 4], DXU[r + 1, i + 8], DXU[r + 1, i + 12])
    ) else
    (
        mixcol_diff_deterministic(DZU[r, round_permutation[i]], DZU[r, round_permutation[i + 4]], DZU[r, round_permutation[i + 8]], DZU[r, round_permutation[i + 12]],
                                DXU[r + 1, i], DXU[r + 1, i + 4], DXU[r + 1, i + 8], DXU[r + 1, i + 12])
    ) endif
);

%##################################################################################
%##################################################################################
%  _____  __  __  _     
% | ____||  \/  || |    
% |  _|  | |\/| || |    
% | |___ | |  | || |___ 
% |_____||_|  |_||_____|

array[0..(RM + RL), 0..15] of var 0..1: DXL;
array[0..(RM + RL - 1), 0..15] of var 0..1: DZL;
constraint forall(r in 0..(RM - 1), i in 0..15) 
(
    DSTKL[r, i] = DXL[r, i]
    /\
    DXL[r, i] = DZL[r, i]
);
constraint forall(r in 0..(RM - 1), i in 0..3)
(
    if r < RM - RML then
    (
        mixcol_lin_deterministic(DXL[r + 1, i], DXL[r + 1, i + 4], DXL[r + 1, i + 8], DXL[r + 1, i + 12],
                                 DZL[r, round_permutation[i]], DZL[r, round_permutation[i + 4]], DZL[r, round_permutation[i + 8]], DZL[r, round_permutation[i + 12]])
    ) else
    (
        mixcol_lin_probabilistic(DXL[r + 1, i], DXL[r + 1, i + 4], DXL[r + 1, i + 8], DXL[r + 1, i + 12],
                                 DZL[r, round_permutation[i]], DZL[r, round_permutation[i + 4]], DZL[r, round_permutation[i + 8]], DZL[r, round_permutation[i + 12]])
    ) endif
);

%##################################################################################
%##################################################################################
%  _____  _     
% | ____|| |    
% |  _|  | |    
% | |___ | |___ 
% |_____||_____|

constraint forall(r in RM..(RM + RL - 1), i in 0..15) 
(
    DSTKL[r, i] = DXL[r, i]
    /\ 
    DXL[r, i] = DZL[r, i]
);
constraint forall(r in RM..(RM + RL - 1), i in 0..3)
(
    mixcol_lin_probabilistic(DXL[r + 1, i], DXL[r + 1, i + 4], DXL[r + 1, i + 8], DXL[r + 1, i + 12],
                             DZL[r, round_permutation[i]], DZL[r, round_permutation[i + 4]], DZL[r, round_permutation[i + 8]], DZL[r, round_permutation[i + 12]])
);
%##################################################################################
constraint sum(i in 0..15)(LANEU[i]) + sum(i in 0..15)(DXU[0, i]) != 0;
constraint sum(r in 0..(RU + RM - 1), i in 0..15)(DXU[r, i]) >= lower_bound[NPT, RU + RM];
constraint sum(i in 0..15)(DXL[RM + RL, i]) != 0;

var 0..(block_size + domain_margine): P0;
var 0..(block_size + domain_margine): Pm;
var 0..(block_size + domain_margine): P1;
constraint P0 = sum(r in 0..(RU - 1), i in 0..15)(DXU[r, i]);
constraint Pm = sum(r in 0..(RM - 1), i in 0..15)(bool2int(DXU[RU + r, i] + DXL[r, i] == 2));
constraint P1 = sum(r in RM..(RM + RL - 1), i in 0..15)(DXL[r, i]);

var 0..(block_size + domain_margine): Pr;
Pr = (WU*1*P0 + WM*1*Pm + WL*1*P1);

%##################################################################################
%##################################################################################
%##################################################################################
%   ___   _      _              _    _               _____                     _    _                    
%  / _ \ | |__  (_)  ___   ___ | |_ (_)__   __ ___  |  ___|_   _  _ __    ___ | |_ (_)  ___   _ __   ___ 
% | | | || '_ \ | | / _ \ / __|| __|| |\ \ / // _ \ | |_  | | | || '_ \  / __|| __|| | / _ \ | '_ \ / __|
% | |_| || |_) || ||  __/| (__ | |_ | | \ V /|  __/ |  _| | |_| || | | || (__ | |_ | || (_) || | | |\__ \
%  \___/ |_.__/_/ | \___| \___| \__||_|  \_/  \___| |_|    \__,_||_| |_| \___| \__||_| \___/ |_| |_||___/
%             |__/                                                                                      


% fix the input/output difference of distinguisher

% array[0..15] of 0..1: alpha = array1d(0..15, [0, 0, 0, 0,
%                                               0, 0, 0, 0,
%                                               0, 0, 0, 0,
%                                               0, 0, 0, 1]);
% constraint forall(i in 0..15)(DXU[0, i] = alpha[i]);

% array[0..15] of 0..1: input_tk_diff = array1d(0..15, [0, 0, 1, 0,
%                                                       0, 0, 0, 0,
%                                                       0, 0, 0, 0,
%                                                       1, 0, 0, 0]);
% constraint forall(i in 0..15)(LANEU[i] = input_tk_diff[i]);

% array[0..15] of 0..1: delta = array1d(0..15, [0, 0, 1, 0,
%                                               0, 0, 1, 0,
%                                               0, 0, 0, 0,
%                                               0, 0, 1, 0]);
% constraint forall(i in 0..15)(DXL[RM + RL, i] = delta[i]);

solve minimize Pr;
%##################################################################################
%##################################################################################
%##################################################################################
%     _                 _  _  _                       _____                     _    _                    
%    / \   _   _ __  __(_)| |(_)  __ _  _ __  _   _  |  ___|_   _  _ __    ___ | |_ (_)  ___   _ __   ___ 
%   / _ \ | | | |\ \/ /| || || | / _` || '__|| | | | | |_  | | | || '_ \  / __|| __|| | / _ \ | '_ \ / __|
%  / ___ \| |_| | >  < | || || || (_| || |   | |_| | |  _| | |_| || | | || (__ | |_ | || (_) || | | |\__ \
% /_/   \_\\__,_|/_/\_\|_||_||_| \__,_||_|    \__, | |_|    \__,_||_| |_| \___| \__||_| \___/ |_| |_||___/
%                                             |___/                                                    
% Auxiliary functions

predicate xor_deterministic(var 0..1: AX, var 0..1: AY, var 0..1: AZ) =     
    % Deterministic propagation
    AZ >= AX
    /\
    AZ >= AY
    /\
    AZ <= AX + AY;

predicate xor_probabilistic(var 0..1: AX, var 0..1: AY, var 0..1: AZ) = 
    % Probabilistic propagation
    AX + AY - AZ >= 0
    /\
    AX - AY + AZ >= 0
    /\
    AY + AZ - AX >= 0;

predicate xor3_deterministic(var 0..1: AX0, var 0..1: AX1, var 0..1: AX2, var 0..1: AY) = 
    AY >= AX0
    /\
    AY >= AX1
    /\
    AY >= AX2
    /\
    AY <= AX0 + AX1 + AX2;

predicate xor3_probabilistic(var 0..1: AX0, var 0..1: AX1, var 0..1: AX2, var 0..1: AY) = 
    let {
        var 0..1: D;        
    } in
    xor_probabilistic(AX0, AX1, D)
    /\
    xor_probabilistic(D, AX2, AY);

predicate mixcol_diff_probabilistic(var 0..1: AX0, var 0..1: AX1, var 0..1: AX2, var 0..1: AX3,
                                    var 0..1: AY0, var 0..1: AY1, var 0..1: AY2, var 0..1: AY3) = 
    % left to right propagation
    AY1 = AX0
    /\
    xor_probabilistic(AX1, AX2, AY2)
    /\
    xor_probabilistic(AX0, AX2, AY3)
    /\
    xor_probabilistic(AY3, AX3, AY0);

predicate mixcol_lin_probabilistic(var 0..1: AX0, var 0..1: AX1, var 0..1: AX2, var 0..1: AX3,
                                   var 0..1: AY0, var 0..1: AY1, var 0..1: AY2, var 0..1: AY3) = 
    % right to left propagation
    xor3_probabilistic(AX0, AX1, AX3, AY0)
    /\
    AY1 = AX2
    /\
    xor3_probabilistic(AX0, AX2, AX3, AY2)
    /\
    AY3 = AX0;

predicate mixcol_diff_deterministic(var 0..1: AX0, var 0..1: AX1, var 0..1: AX2, var 0..1: AX3,
                                    var 0..1: AY0, var 0..1: AY1, var 0..1: AY2, var 0..1: AY3) =
    % left to right propagation
    AY1 = AX0
    /\
    xor_deterministic(AX1, AX2, AY2)
    /\
    xor_deterministic(AX0, AX2, AY3)
    /\
    xor_deterministic(AY3, AX3, AY0);

predicate mixcol_lin_deterministic(var 0..1: AX0, var 0..1: AX1, var 0..1: AX2, var 0..1: AX3,
                                   var 0..1: AY0, var 0..1: AY1, var 0..1: AY2, var 0..1: AY3) =
    % right to left propagation
    xor3_deterministic(AX0, AX1, AX3, AY0)
    /\
    AY1 = AX2
    /\
    xor3_deterministic(AX0, AX2, AX3, AY2)
    /\
    AY3 = AX0;
