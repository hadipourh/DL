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

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
%  __  __             _        _  _                 _    _             ____          _                                     __   ____   _____  ____   ____   _____  _   _  _____ 
% |  \/  |  ___    __| |  ___ | |(_) _ __    __ _  | |_ | |__    ___  / ___|        | |__    ___ __  __ ___  ___    ___   / _| / ___| | ____||  _ \ |  _ \ | ____|| \ | ||_   _|
% | |\/| | / _ \  / _` | / _ \| || || '_ \  / _` | | __|| '_ \  / _ \ \___ \  _____ | '_ \  / _ \\ \/ // _ \/ __|  / _ \ | |_  \___ \ |  _|  | |_) || |_) ||  _|  |  \| |  | |  
% | |  | || (_) || (_| ||  __/| || || | | || (_| | | |_ | | | ||  __/  ___) ||_____|| |_) || (_) |>  <|  __/\__ \ | (_) ||  _|  ___) || |___ |  _ < |  __/ | |___ | |\  |  | |  
% |_|  |_| \___/  \__,_| \___||_||_||_| |_| \__, |  \__||_| |_| \___| |____/        |_.__/  \___//_/\_\\___||___/  \___/ |_|   |____/ |_____||_| \_\|_|    |_____||_| \_|  |_|  
%                                           |___/                                                                                                                               
% Modeling the S-boxes of SERPENT

% We use S-box Analyzer to encode DDT: 
% https://github.com//sboxanalyzer

/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.02 seconds
Number of constraints: 54
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1

*/
predicate s0_diff(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                  var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                  var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - a0 + p0 + p1 >= 0
    /\
    - a3 + p0 + p1 >= 0
    /\
    b0 + b2 - b3 + p0 >= 0
    /\
    a1 + b0 + p0 - p1 >= 0
    /\
    a0 + a1 + a2 + a3 - b0 >= 0
    /\
    a0 + a2 + a3 + b0 - b2 >= 0
    /\
    - a0 - a2 - b0 + b2 - b3 >= -3
    /\
    - a2 - a3 - b0 - b2 + b3 >= -3
    /\
    a0 - a1 + a3 + b2 + b3 >= 0
    /\
    - a1 - b0 - b1 + b2 + b3 >= -2
    /\
    a1 - a2 + b1 + b2 + b3 >= 0
    /\
    a0 - a1 - a3 - b0 + p0 >= -2
    /\
    - a0 - a1 + a3 - b0 + p0 >= -2
    /\
    a1 + a2 - a3 - b3 + p0 >= -1
    /\
    - b0 + b1 - b2 - b3 + p0 >= -2
    /\
    - a0 + a1 + a2 + b3 + p0 >= 0
    /\
    - a0 - a3 + b0 + b3 + p0 >= -1
    /\
    - a1 + b0 + b1 + b3 + p0 >= 0
    /\
    - a0 - a1 - b2 + b3 + p0 >= -2
    /\
    a0 - a2 - b2 + b3 + p0 >= -1
    /\
    - a0 - a1 + a2 - a3 + p1 >= -2
    /\
    a0 - a1 + a3 - b0 + p1 >= -1
    /\
    - a0 - a3 - b2 - b3 + p1 >= -3
    /\
    - a0 - a3 + b2 + b3 + p1 >= -1
    /\
    b0 + b2 + b3 - p0 + p1 >= 0
    /\
    a0 - a1 - a3 + b0 - b1 - b2 >= -3
    /\
    - a0 - a1 + a3 + b0 + b1 - b2 >= -2
    /\
    - a0 + a1 + a2 + a3 + b1 + b2 >= 0
    /\
    a0 + a1 - a3 - b0 - b1 - b3 >= -3
    /\
    - a1 + a2 + a3 + b0 - b1 - b3 >= -2
    /\
    a0 - a1 - a3 + b0 + b1 - b3 >= -2
    /\
    - a0 - a1 - a2 + b0 - b2 - b3 >= -4
    /\
    a0 - a1 - b0 - b1 - b2 - b3 >= -4
    /\
    - a0 + a3 - b0 - b1 - b2 - b3 >= -4
    /\
    - a1 + a2 - a3 - b0 + b2 - b3 >= -3
    /\
    a0 + a1 - a2 + b0 - b1 + b3 >= -1
    /\
    a0 + a1 + a2 - a3 + b1 + b3 >= 0
    /\
    - a0 - a1 + a2 - b0 - b2 + b3 >= -3
    /\
    a1 + a2 + a3 + b1 - b2 + b3 >= 0
    /\
    - a0 + a3 - b0 - b1 + b3 + p0 >= -2
    /\
    a0 - a3 - b0 + b2 + b3 + p0 >= -1
    /\
    a0 - a2 + a3 - b2 - b3 + p1 >= -2
    /\
    a0 + a2 + b1 + b2 - b3 + p1 >= 0
    /\
    a0 + a3 - b1 + b2 + b3 + p1 >= 0
    /\
    a0 + a1 + a3 + b0 - p0 + p1 >= 0
    /\
    a0 + a1 - a2 - a3 - b1 + b2 - b3 >= -3
    /\
    - a0 + a1 - a2 + a3 - b1 - b2 + b3 >= -3
    /\
    a0 - a2 + a3 - b0 + b2 - b3 + p0 >= -2
    /\
    a0 + a1 - a3 - b0 - b1 - b2 + p1 >= -3
    /\
    - a0 - a2 + a3 + b0 - b1 + b2 + p1 >= -2
    /\
    a0 + a1 - a3 - b0 + b2 - b3 + p1 >= -2
    /\
    a1 + a2 + a3 - b0 + b2 - b3 + p1 >= -1
    /\
    - a0 + a1 + a3 - b0 - b2 + b3 + p1 >= -2
;

% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.03 seconds
Number of constraints: 53
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1
*/
predicate s1_diff(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                  var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                  var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    a0 + a2 - a3 + p0 >= 0
    /\
    a0 - a2 + a3 + p0 >= 0
    /\
    a0 - a2 + b1 + p0 >= 0
    /\
    - a0 - a1 + a2 + a3 - b1 >= -2
    /\
    - a0 - a2 + a3 - b0 - b3 >= -3
    /\
    a1 + a2 + a3 + b1 - b3 >= 0
    /\
    - a0 + a2 - a3 - b2 - b3 >= -3
    /\
    a0 + b0 - b1 + b2 + b3 >= 0
    /\
    - a0 + a1 - a2 - b0 + p0 >= -2
    /\
    a0 + a3 - b0 - b2 + p0 >= -1
    /\
    - a2 - a3 + b0 - b2 + p0 >= -2
    /\
    - a0 + b0 - b1 - b2 + p0 >= -2
    /\
    - a0 - b0 - b1 + b2 + p0 >= -2
    /\
    a0 - a2 - b2 - b3 + p0 >= -2
    /\
    - a2 + a3 + b2 - b3 + p0 >= -1
    /\
    - b0 + b1 - b2 + b3 + p0 >= -1
    /\
    - a1 - a3 - b0 + b2 - p1 >= -3
    /\
    - a0 - a1 - a2 - a3 + p1 >= -3
    /\
    a0 - a2 - a3 - b1 + p1 >= -2
    /\
    - a2 - a3 - b0 - b2 + p1 >= -3
    /\
    a2 + a3 - b0 - b2 + p1 >= -1
    /\
    - a0 + b0 - b1 + b2 + p1 >= -1
    /\
    - b0 - b1 - b2 + b3 + p1 >= -2
    /\
    a0 + a2 + a3 - p0 + p1 >= 0
    /\
    a0 + a1 + a3 + b0 - b1 - b2 >= -1
    /\
    - a0 + a1 - a2 - a3 + b0 + b2 >= -2
    /\
    a0 + a1 - a3 - b0 - b1 + b2 >= -2
    /\
    a0 - a1 + a2 + b1 + b2 - b3 >= -1
    /\
    - a0 + a2 - a3 - b0 - b1 + b3 >= -3
    /\
    - a0 - a2 + a3 - b1 - b2 + b3 >= -3
    /\
    a1 + a2 + b0 + b1 - b2 + b3 >= 0
    /\
    a1 + a3 + b0 + b1 - b2 + b3 >= 0
    /\
    a1 + a2 - b0 + b1 + b2 + b3 >= 0
    /\
    a1 + a3 + b1 + b2 + b3 - p0 >= 0
    /\
    - a1 + a3 + b0 + b1 - b2 + p0 >= -1
    /\
    a0 - a1 + a3 + b0 + b2 + p0 >= 0
    /\
    a1 + a2 + a3 - b1 + b2 + p0 >= 0
    /\
    a2 + a3 - b0 + b1 + b2 + p0 >= 0
    /\
    - a0 - a2 - a3 - b1 + b2 - p1 >= -4
    /\
    b0 + b1 + b2 + b3 + p0 - p1 >= 0
    /\
    a0 + b0 + b1 + b2 - p0 + p1 >= 0
    /\
    - a1 - a2 + a3 + b0 + b1 - b2 - b3 >= -3
    /\
    - a0 + a2 - a3 + b0 + b2 - b3 + p0 >= -2
    /\
    - a0 - a2 + a3 + b0 + b1 - b2 + p1 >= -2
    /\
    - a0 + a2 - a3 - b0 + b1 + b2 + p1 >= -2
    /\
    - a0 + a2 + a3 + b0 + b1 + b2 + p1 >= 0
    /\
    a0 - a1 + a2 + b0 - b2 - b3 + p1 >= -2
    /\
    a0 - a1 + a3 - b0 + b2 - b3 + p1 >= -2
    /\
    - a0 + a2 - a3 + b0 + b1 + b3 + p1 >= -1
    /\
    a0 - a1 + a2 - b1 - b2 + b3 + p1 >= -2
    /\
    a0 - a1 + a3 - b1 + b2 + b3 + p1 >= -1
    /\
    - a0 - a2 + a3 + b1 + b2 + b3 + p1 >= -1
;

/*% #############################################################################################################################################
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.02 seconds
Number of constraints: 56
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1
*/
predicate s2_diff(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                  var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                  var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - b3 + p0 + p1 >= 0
    /\
    - a3 - b0 - b3 - p1 >= -3
    /\
    - a0 + a1 - a2 + a3 - b3 >= -2
    /\
    - a1 + a2 + a3 + b0 - b3 >= -1
    /\
    - a0 - a2 + a3 + b1 - b3 >= -2
    /\
    - a0 + a1 + a2 + b0 + b3 >= 0
    /\
    a3 - b0 - b1 + b2 + b3 >= -1
    /\
    - a0 + a3 + b0 + b1 + p0 >= 0
    /\
    a2 - a3 + b0 - b2 + p0 >= -1
    /\
    a1 + a2 + b1 - b2 + p0 >= 0
    /\
    a2 + b0 - b1 + b2 + p0 >= 0
    /\
    a1 - a2 - a3 + b3 + p0 >= -1
    /\
    - a0 - a1 - b0 + b3 + p0 >= -2
    /\
    - a0 - a1 - b1 + b3 + p0 >= -2
    /\
    - a0 - b0 - b1 + b3 + p0 >= -2
    /\
    a0 - a2 - b2 + b3 + p0 >= -1
    /\
    - a2 + a3 + b2 + b3 + p0 >= 0
    /\
    a3 + b0 + b1 + b2 - p1 >= 0
    /\
    a0 - a1 - a2 - b3 - p1 >= -3
    /\
    - a1 - a2 - a3 - b3 - p1 >= -4
    /\
    a0 - a1 + a2 + a3 + p1 >= 0
    /\
    a3 + b0 - b1 - b2 + p1 >= -1
    /\
    a3 - b0 + b1 - b2 + p1 >= -1
    /\
    - b0 - b1 + b2 - b3 + p1 >= -2
    /\
    a0 + a1 + a3 - p0 + p1 >= 0
    /\
    b0 + b1 + b2 - p0 + p1 >= 0
    /\
    - a0 + a1 - a2 - a3 + b1 - b2 >= -3
    /\
    a0 + a1 - a2 - b0 + b1 - b2 >= -2
    /\
    - a1 + a2 - a3 - b0 + b1 - b2 >= -3
    /\
    a1 - a2 - a3 + b0 - b1 - b3 >= -3
    /\
    a0 + a1 + a2 + b0 + b1 - b3 >= 0
    /\
    a0 + a1 + a2 - b0 + b2 - b3 >= -1
    /\
    a0 - a1 - a3 - b0 + b2 - b3 >= -3
    /\
    - a0 - a2 + b0 - b1 + b2 - b3 >= -3
    /\
    - a1 + a2 - b0 + b1 + b2 - b3 >= -2
    /\
    a0 + a2 + b0 + b1 - b2 + b3 >= 0
    /\
    - a0 + a2 - a3 + b0 + b2 + b3 >= -1
    /\
    - a0 - a1 - a3 + b1 + b2 + b3 >= -2
    /\
    a0 - a2 - a3 + b1 + b2 + b3 >= -1
    /\
    a0 + a1 + b0 + b2 + b3 - p0 >= 0
    /\
    a0 - a2 - b0 - b1 + b2 + p0 >= -2
    /\
    a2 + a3 - b0 + b1 + b2 + p0 >= 0
    /\
    a0 + a2 - a3 + b2 + b3 + p0 >= 0
    /\
    a2 + a3 - b0 - b1 - b2 - p1 >= -3
    /\
    a0 + a1 + b0 - b1 - b2 - p1 >= -2
    /\
    - a1 - a2 - a3 + b0 + b3 + p1 >= -2
    /\
    a1 + a2 - b1 - b2 + b3 + p1 >= -1
    /\
    - a0 + a1 + a2 - a3 - b0 - b1 - b2 >= -4
    /\
    - a0 - a1 - a2 - b0 - b1 - b2 - b3 >= -6
    /\
    a0 - a2 - a3 + b0 + b1 - b2 - b3 >= -3
    /\
    - a0 + a1 + a2 - a3 + b0 + b2 - b3 >= -2
    /\
    a0 - a1 - a3 - b0 - b1 - b2 + b3 >= -4
    /\
    a1 + a2 - b0 + b1 + b2 + b3 + p1 >= 0
    /\
    - a0 - a2 - a3 - b0 - b1 - b2 + b3 - p0 + p1 >= -6
    /\
    - a1 + a2 - b1 - b2 - b3 >= -3
;

% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.02 seconds
Number of constraints: 64
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1
*/
predicate s3_diff(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                  var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                  var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    a2 - a3 - b2 + p0 >= -1
    /\
    a2 + b1 - b3 + p0 >= 0
    /\
    - a3 + b1 + b3 + p0 >= 0
    /\
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    a0 + a1 + a3 - b1 + b3 >= 0
    /\
    b0 + b1 + b2 + b3 - p0 >= 0
    /\
    - a0 + a2 + a3 + b0 + p0 >= 0
    /\
    - a0 - a2 + a3 - b1 + p0 >= -2
    /\
    a0 - a1 - b0 - b1 + p0 >= -2
    /\
    a1 + a2 - b0 + b1 + p0 >= 0
    /\
    - a1 + a3 - b0 - b2 + p0 >= -2
    /\
    a3 + b0 + b1 - b2 + p0 >= 0
    /\
    a0 + a1 - b0 + b2 + p0 >= 0
    /\
    - a0 - a2 - a3 - b3 + p0 >= -3
    /\
    - a2 + b0 - b1 - b3 + p0 >= -2
    /\
    a1 + b1 - b2 - b3 + p0 >= -1
    /\
    a0 + b1 + b2 - b3 + p0 >= 0
    /\
    - a2 + a3 + b0 + b3 + p0 >= 0
    /\
    - a0 - b1 + b2 + b3 + p0 >= -1
    /\
    a0 - a2 - a3 - b0 + b1 - b2 >= -3
    /\
    - a0 + a1 + a2 - a3 - b1 + b2 >= -2
    /\
    a0 - a1 + a2 + b0 - b1 + b2 >= -1
    /\
    a0 - a1 + a2 + a3 + b1 + b2 >= 0
    /\
    a0 - a2 - a3 + b0 + b1 + b2 >= -1
    /\
    - a0 - a1 + a2 + a3 - b0 - b3 >= -3
    /\
    a1 - a2 + a3 + b0 + b1 - b3 >= -1
    /\
    - a0 + a2 + a3 - b0 - b2 - b3 >= -3
    /\
    - a0 - a1 + b0 + b1 - b2 - b3 >= -3
    /\
    - a0 + a1 + a3 + b0 + b2 - b3 >= -1
    /\
    - a1 + a3 - b0 + b1 + b2 - b3 >= -2
    /\
    - a0 - a2 + a3 - b0 + b1 + b3 >= -2
    /\
    a1 - a3 - b0 - b1 - b2 + b3 >= -3
    /\
    - a2 - a3 + b0 - b1 - b2 + b3 >= -3
    /\
    a1 + a2 + b0 - b1 + b2 + b3 >= 0
    /\
    a0 + a1 + a2 + b1 + b3 - p0 >= 0
    /\
    a1 + a2 + a3 + b1 + b3 - p0 >= 0
    /\
    a0 + a3 + b0 + b1 + p0 - p1 >= 0
    /\
    - a0 - a1 + a2 - a3 - b0 + p1 >= -3
    /\
    a0 + a1 - a2 - a3 - b1 + p1 >= -2
    /\
    a0 + a1 - a3 + b0 + b2 + p1 >= 0
    /\
    - a0 - a2 + a3 + b1 + b2 + p1 >= -1
    /\
    - a1 - a3 - b0 - b1 - b3 + p1 >= -4
    /\
    a1 + a3 - b0 - b1 - b3 + p1 >= -2
    /\
    a1 + a2 + b0 - b1 - b3 + p1 >= -1
    /\
    - a0 + a2 - b1 + b2 - b3 + p1 >= -2
    /\
    a0 - a1 + a2 + b0 + b3 + p1 >= 0
    /\
    a0 - a2 + a3 + b1 + b3 + p1 >= 0
    /\
    a1 - b0 - b1 - b2 + b3 + p1 >= -2
    /\
    - a1 - a2 + b0 + b2 + b3 + p1 >= -1
    /\
    - a0 - a1 - a2 + a3 + b0 - b1 - b2 >= -4
    /\
    a0 - a1 - a2 + a3 - b0 + b1 + p1 >= -2
    /\
    a0 - a1 - a2 - a3 - b2 - b3 + p1 >= -4
    /\
    a0 - a1 + a2 - b1 - b2 - b3 + p1 >= -3
    /\
    - a0 - a1 - a3 - b1 - b2 + b3 + p1 >= -4
    /\
    a0 - a1 - a2 - a3 + b2 + b3 + p1 >= -2
    /\
    - a0 - a1 + a3 - b0 + b2 + b3 + p1 >= -2
    /\
    - a1 + a2 + a3 + b0 + b1 - b2 - p0 + p1 >= -2
    /\
    - a0 + a1 - a2 - a3 - b2 - b3 - p0 + p1 >= -5
    /\
    - a0 - a1 - a2 - a3 + b2 - b3 - p0 + p1 >= -5
    /\
    a0 - a1 - a2 + a3 - b1 + b2 - b3 - p0 + p1 >= -4
    /\
    a1 + a2 - a3 - b0 + b1 + b2 - b3 - p0 + p1 >= -3
    /\
    a0 - a1 + a3 - b0 - b1 - b2 + b3 - p0 + p1 >= -4
    /\
    - a0 + a1 - a2 + b1 + b2 + b3 >= -1
;

% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 60
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1
*/
predicate s4_diff(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                  var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                  var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - a0 + p0 + p1 >= 0
    /\
    a0 - a3 + b1 + p0 >= 0
    /\
    - b0 + b2 - b3 + p0 >= -1
    /\
    - a2 - b0 + b3 + p0 >= -1
    /\
    - a2 + b1 + b3 + p0 >= 0
    /\
    - a0 - a1 - a2 - a3 - b3 >= -4
    /\
    - a1 - a2 - a3 + b2 - b3 >= -3
    /\
    - a0 + a1 - a3 - b0 + b3 >= -2
    /\
    - a0 + a2 - a3 - b0 + b3 >= -2
    /\
    a2 + a3 + b0 - b2 + b3 >= 0
    /\
    - a0 + a3 + b0 + b2 + b3 >= 0
    /\
    b0 + b1 + b2 + b3 - p0 >= 0
    /\
    a0 - a2 - a3 + b0 + p0 >= -1
    /\
    a1 - a2 - a3 - b1 + p0 >= -2
    /\
    - a1 - a2 + a3 - b2 + p0 >= -2
    /\
    a1 + a3 + b1 - b2 + p0 >= 0
    /\
    a0 + a1 - b1 + b2 + p0 >= 0
    /\
    a0 + a3 - b1 - b3 + p0 >= -1
    /\
    - a0 + a2 - b0 + b1 - p1 >= -2
    /\
    - a0 + a1 - b1 - b2 - p1 >= -3
    /\
    - a0 - a1 + a3 + b2 - p1 >= -2
    /\
    a0 + a2 + a3 + b3 - p1 >= 0
    /\
    a2 + b0 + b2 + p0 - p1 >= 0
    /\
    a0 + a1 + a3 - b0 + b1 - b2 >= -1
    /\
    a0 + a1 + a2 - b0 + b1 + b2 >= 0
    /\
    a0 - a2 + a3 + b0 + b1 + b2 >= 0
    /\
    - a0 - a1 + a2 + b0 - b2 - b3 >= -3
    /\
    - a1 - a2 + a3 - b1 - b2 - b3 >= -4
    /\
    a1 + a2 - b0 - b1 - b2 - b3 >= -3
    /\
    a0 + a1 + a3 - b1 + b2 - b3 >= -1
    /\
    a0 - a3 + b0 - b1 + b2 - b3 >= -2
    /\
    a1 + a3 + b0 + b1 + b2 - b3 >= 0
    /\
    - a0 - a2 + a3 - b0 - b1 + b3 >= -3
    /\
    a0 + a1 - a2 + a3 + b1 + b3 >= 0
    /\
    a0 - a3 + b0 - b1 - b2 + b3 >= -2
    /\
    - a1 + a2 - a3 - b1 - b2 + p0 >= -3
    /\
    a0 - a2 - a3 - b0 - b1 + p1 >= -3
    /\
    a0 - a1 + a2 + a3 + b1 + p1 >= 0
    /\
    a0 + a1 + a2 - b1 - b2 + p1 >= -1
    /\
    - a0 + a1 - a3 + b1 - b2 + p1 >= -2
    /\
    - a0 + a3 - b1 - b2 - b3 + p1 >= -3
    /\
    - a1 - b0 + b1 - b2 - b3 + p1 >= -3
    /\
    - a0 - a2 - a3 + b0 + b3 + p1 >= -2
    /\
    a0 - a2 + a3 + b0 + b3 + p1 >= 0
    /\
    a0 - a1 + a3 - b1 + b3 + p1 >= -1
    /\
    - a1 + a2 + a3 - b2 + b3 + p1 >= -1
    /\
    - a1 + a2 + b1 - b2 + b3 + p1 >= -1
    /\
    a0 + a1 + a2 + b0 - p0 + p1 >= 0
    /\
    - a0 - a1 + a2 - b0 - b1 + b2 - b3 >= -4
    /\
    a1 + a2 + a3 - b0 - b1 + b2 + p1 >= -1
    /\
    - a0 - a2 + a3 + b0 - b1 + b2 + p1 >= -2
    /\
    - a0 - a1 + a2 + a3 - b0 - b3 + p1 >= -3
    /\
    - a0 + a1 - a2 + b0 + b1 - b3 + p1 >= -2
    /\
    a0 + a2 - a3 - b1 - b2 - b3 + p1 >= -3
    /\
    - a0 + a1 - a2 + b1 + b2 - b3 + p1 >= -2
    /\
    a1 + a2 + a3 - b0 + b2 + b3 + p1 >= 0
    /\
    a0 - a1 + a2 - b1 + b2 + b3 + p1 >= -1
    /\
    a1 - a3 + b0 + b1 - b2 - b3 >= -2
    /\
    a1 + a2 - a3 + b0 - b1 + b2 - b3 >= -2
;

% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 60
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1
*/
predicate s5_diff(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                  var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                  var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - a2 + a3 + b1 + p0 >= 0
    /\
    a0 + b0 - b3 + p0 >= 0
    /\
    - a0 - b2 + b3 + p0 >= -1
    /\
    - a0 + b1 + b3 - p1 >= -1
    /\
    a2 + b1 + b3 - p1 >= 0
    /\
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    - a0 - a1 - a2 + b0 - b3 >= -3
    /\
    a0 - a2 - a3 - b2 + b3 >= -2
    /\
    a0 + a2 + a3 - b2 + b3 >= 0
    /\
    a0 + a3 - b0 + b2 + b3 >= 0
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    a1 + a2 - a3 - b0 + p0 >= -1
    /\
    - a0 + a1 - a2 - b1 + p0 >= -2
    /\
    a0 - a2 - a3 - b1 + p0 >= -2
    /\
    - a0 + a2 + a3 - b3 + p0 >= -1
    /\
    - a1 - a3 + b1 - b3 + p0 >= -2
    /\
    a2 - a3 + b2 + b3 + p0 >= 0
    /\
    a0 + a1 - a2 - a3 - b0 - b2 >= -3
    /\
    - a1 + a2 + a3 + b0 + b1 + b2 >= 0
    /\
    - a0 - a1 + a2 - b0 - b1 - b3 >= -4
    /\
    - a0 - a1 - a3 - b1 + b2 - b3 >= -4
    /\
    - a1 - a3 - b0 - b1 + b2 - b3 >= -4
    /\
    a1 + a2 + b0 + b1 + b2 - b3 >= 0
    /\
    - a0 - a1 + a3 - b0 - b1 + b3 >= -3
    /\
    - a0 + a1 + a2 + a3 + b1 + b3 >= 0
    /\
    a0 + a1 + a3 + b0 + b1 - p0 >= 0
    /\
    a0 - a1 + a3 - b0 - b1 + p0 >= -2
    /\
    - a1 - a2 + a3 - b1 + b2 + p0 >= -2
    /\
    - a0 - a3 - b1 - b2 - b3 + p0 >= -4
    /\
    a0 - a1 + a2 + b0 + b3 + p0 >= 0
    /\
    a1 + a3 + b0 - b1 + b3 + p0 >= 0
    /\
    a0 - a1 + a2 - a3 - b0 + p1 >= -2
    /\
    - a0 + a1 - a2 - a3 + b1 + p1 >= -2
    /\
    a0 - a1 + a2 + a3 + b1 + p1 >= 0
    /\
    a0 - a1 - b0 + b1 - b2 + p1 >= -2
    /\
    a0 - a2 + a3 - b1 + b2 + p1 >= -1
    /\
    - a2 - a3 - b0 + b1 + b2 + p1 >= -2
    /\
    a2 - a3 - b0 - b1 - b3 + p1 >= -3
    /\
    - a0 - a1 - b0 - b2 - b3 + p1 >= -4
    /\
    - a0 - a2 - a3 + b2 + b3 + p1 >= -2
    /\
    - a0 + a2 + a3 + b2 + b3 + p1 >= 0
    /\
    - a0 - a3 + b0 + b2 + b3 + p1 >= -1
    /\
    a0 + a2 + b2 + b3 - p0 + p1 >= 0
    /\
    a1 - a2 + a3 - b0 - b1 - b2 - b3 >= -4
    /\
    - a0 + a1 - a3 + b0 + b1 - b2 - b3 >= -3
    /\
    a0 + a1 + a2 + b0 - b1 - b2 + p1 >= -1
    /\
    - a2 + a3 - b0 - b1 - b2 - b3 + p1 >= -4
    /\
    a0 - a1 + a2 + b1 - b2 - b3 + p1 >= -2
    /\
    a2 + a3 - b0 + b1 - b2 - b3 + p1 >= -2
    /\
    - a0 + a1 + a2 - a3 + b2 - b3 + p1 >= -2
    /\
    a1 - a2 + a3 - b0 - b1 + b3 + p1 >= -2
    /\
    a0 + a1 - a3 + b0 - b2 + b3 + p1 >= -1
    /\
    - a1 + a3 + b0 - b1 - b2 + b3 + p1 >= -2
    /\
    - a0 + a1 + a3 + b0 - b1 - b3 - p0 + p1 >= -3
    /\
    - a0 + a1 - a2 - b1 - b2 + b3 - p0 + p1 >= -4
    /\
    - a0 + a2 - a3 - b1 - b2 + b3 - p0 + p1 >= -4
    /\
    a0 - a1 - a3 + b0 - b1 - b2 - b3 - p0 + p1 >= -5
    /\
    a0 + a1 - a2 + b0 - b1 + b2 - b3 >= -2
    /\
    a1 - a2 - b0 + b1 + b2 - b3 >= -2
;

% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 52
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1
*/
predicate s6_diff(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                  var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                  var 0..1: p0, var 0..1: p1) =
    - p0 - p1 >= -1
    /\
    - a0 + a3 + b2 + p0 >= 0
    /\
    a0 + a1 + a2 + a3 - b0 >= 0
    /\
    a0 - a1 - a2 + a3 - b2 >= -2
    /\
    - a0 + a3 - b0 - b1 - b2 >= -3
    /\
    a0 - a3 - b0 - b2 - b3 >= -3
    /\
    - a2 - b0 + b1 - b2 + b3 >= -2
    /\
    a2 - b0 + b1 + b2 + b3 >= 0
    /\
    a0 - a3 + b0 - b1 + p0 >= -1
    /\
    - a0 + a1 - a3 - b2 + p0 >= -2
    /\
    - a3 - b0 - b1 + b2 + p0 >= -2
    /\
    - a0 + a2 + b0 - b3 + p0 >= -1
    /\
    a2 + a3 + b0 - b3 + p0 >= 0
    /\
    - a2 + b1 - b2 - b3 + p0 >= -2
    /\
    - a3 + b1 + b2 - b3 + p0 >= -1
    /\
    - a0 - a3 + b0 + b3 + p0 >= -1
    /\
    - a2 - b1 - b2 + b3 + p0 >= -2
    /\
    a3 - b1 + b2 - b3 - p1 >= -2
    /\
    - a0 - a1 - a3 - b2 + p1 >= -3
    /\
    - a0 - a2 - a3 + b2 + p1 >= -2
    /\
    - a0 - a3 - b1 - b3 + p1 >= -3
    /\
    a0 + a3 - b1 - b3 + p1 >= -1
    /\
    - a2 + b0 - b1 - b3 + p1 >= -2
    /\
    - a0 - a3 + b1 + b3 + p1 >= -1
    /\
    - a1 + b0 + b1 + b3 + p1 >= 0
    /\
    a0 + a3 + b2 - p0 + p1 >= 0
    /\
    a0 + a1 + a2 + b0 + b1 - b3 >= 0
    /\
    - a0 - a2 + a3 + b0 - b2 - b3 >= -3
    /\
    a1 - a2 + a3 + b1 + b2 - b3 >= -1
    /\
    a1 + a2 + a3 + b0 - b1 + b3 >= 0
    /\
    - a1 + a2 - a3 - b1 - b2 + b3 >= -3
    /\
    a0 + a1 - a2 - b1 + b2 + b3 >= -1
    /\
    - a1 + a2 - b0 + b1 - b3 + p0 >= -2
    /\
    a0 + a2 + a3 - b1 + b3 + p0 >= 0
    /\
    a0 + a2 - a3 + b1 + b3 + p0 >= 0
    /\
    a0 - a2 + b1 + b2 + b3 + p0 >= 0
    /\
    - a0 + a1 + a3 + b0 + b1 + p1 >= 0
    /\
    a0 + a1 - a3 + b0 + b3 + p1 >= 0
    /\
    - a1 - a3 + b0 - b2 + b3 + p1 >= -2
    /\
    a0 + a3 + b1 - b2 + b3 + p1 >= 0
    /\
    - a1 + a2 + a3 - b0 + b1 + b2 - b3 >= -2
    /\
    a0 - a1 + a2 - a3 - b0 - b1 + b3 >= -3
    /\
    - a0 + a2 + a3 + b1 - b2 + b3 + p0 >= -1
    /\
    a0 + a3 + b0 + b1 + b3 + p0 - p1 >= 0
    /\
    - a1 + a2 - a3 + b1 - b2 - b3 + p1 >= -3
    /\
    - a0 + a2 + a3 + b1 - b2 - b3 + p1 >= -2
    /\
    a0 - a1 - a2 + b1 + b2 - b3 + p1 >= -2
    /\
    a0 - a1 - b0 + b1 + b2 - b3 + p1 >= -2
    /\
    - a0 - a1 + a2 - b1 - b2 + b3 + p1 >= -3
    /\
    a0 + a1 + a2 - b1 - b2 + b3 + p1 >= -1
    /\
    - a1 - a2 + a3 - b1 + b2 + b3 + p1 >= -2
    /\
    - a1 + a3 - b0 - b1 + b2 + b3 + p1 >= -2
;

% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 63
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1
*/
predicate s7_diff(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                  var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                  var 0..1: p0, var 0..1: p1) =
    - p0 - p1 >= -1
    /\
    a0 + a2 - b0 + p0 >= 0
    /\
    - a0 + a2 + b1 + p0 >= 0
    /\
    - a1 - b0 + b1 + p0 >= -1
    /\
    - a1 + a3 - b3 + p0 >= -1
    /\
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    a0 - a2 + b0 + b2 + b3 >= 0
    /\
    - a0 - a2 + a3 - b1 + p0 >= -2
    /\
    a0 + a1 - a2 + b1 + p0 >= 0
    /\
    - a1 - a3 + b0 - b2 + p0 >= -2
    /\
    - a0 - a3 - b1 - b2 + p0 >= -3
    /\
    - a0 - a1 + a2 + b2 + p0 >= -1
    /\
    a1 - b0 - b1 - b3 + p0 >= -2
    /\
    - a2 - b1 + b2 - b3 + p0 >= -2
    /\
    a0 - a2 - a3 + b3 + p0 >= -1
    /\
    a3 + b0 - b1 + b3 + p0 >= 0
    /\
    - a0 + a1 - b2 + b3 + p0 >= -1
    /\
    a1 - a3 + b2 + b3 + p0 >= 0
    /\
    - a0 - a1 - a2 - a3 + b0 + b1 >= -3
    /\
    - a1 + a2 + a3 + b0 + b1 - b2 >= -1
    /\
    - a0 - a1 - a2 + a3 + b1 + b2 >= -2
    /\
    - a1 + a2 - a3 + b0 + b1 + b2 >= -1
    /\
    a1 - a2 + a3 + b0 + b1 + b2 >= 0
    /\
    a0 + a2 - a3 + b0 - b1 - b3 >= -2
    /\
    - a0 + a1 + a3 + b0 + b2 - b3 >= -1
    /\
    a1 - a2 - b0 + b1 + b2 - b3 >= -2
    /\
    a0 + a2 + b0 + b1 + b2 - b3 >= 0
    /\
    - a1 - a2 + a3 - b0 - b1 + b3 >= -3
    /\
    a1 + a2 + b0 + b1 - b2 + b3 >= 0
    /\
    - a1 + a2 + a3 + b0 + b2 + b3 >= 0
    /\
    a1 + a2 + a3 - b1 + b2 + b3 >= 0
    /\
    a0 + a2 + a3 + b1 + b3 - p0 >= 0
    /\
    a1 + a2 + b1 + b2 + b3 - p0 >= 0
    /\
    a0 + a1 + a3 + b0 + p0 - p1 >= 0
    /\
    - a0 - a2 - a3 - b0 - b2 + p1 >= -4
    /\
    a1 + a2 - a3 + b0 - b2 + p1 >= -1
    /\
    a0 + a1 + a3 - b1 - b2 + p1 >= -1
    /\
    a0 - a1 - a2 - a3 + b2 + p1 >= -2
    /\
    - a0 - a2 - a3 + b0 + b2 + p1 >= -2
    /\
    - a1 - a2 - b0 - b1 + b2 + p1 >= -3
    /\
    - a0 + a1 - a2 + b1 - b3 + p1 >= -2
    /\
    - a3 - b0 + b1 - b2 - b3 + p1 >= -3
    /\
    a0 + a2 + b0 - b1 + b3 + p1 >= 0
    /\
    - a2 + a3 + b1 - b2 + b3 + p1 >= -1
    /\
    a0 - a1 - a3 + b2 + b3 + p1 >= -1
    /\
    a1 + a3 - b0 + b2 + b3 + p1 >= 0
    /\
    - a2 - a3 - b1 + b2 + b3 + p1 >= -2
    /\
    a0 + a1 - a2 - a3 + b0 - b2 - b3 >= -3
    /\
    - a0 + a1 - a2 + b0 + b1 - b2 - b3 >= -3
    /\
    - a0 - a1 - a3 - b0 - b1 + b2 + p1 >= -4
    /\
    - a0 + a1 + a2 + b0 - b1 - b3 + p1 >= -2
    /\
    - a0 + a1 - a3 - b0 + b1 - b3 + p1 >= -3
    /\
    a0 - a1 - a2 - b0 - b2 - b3 + p1 >= -4
    /\
    a0 + a1 - b0 - b1 - b2 + b3 + p1 >= -2
    /\
    - a0 + a2 + a3 - b0 - b1 - b2 + b3 + p1 >= -3
    /\
    - a0 - a1 + a2 + a3 - b1 - b3 - p0 + p1 >= -4
    /\
    - a0 + a1 - b0 - b1 - b2 - b3 - p0 + p1 >= -5
    /\
    - a0 - a1 + a2 + b1 - b2 - b3 - p0 + p1 >= -4
    /\
    - a1 + a2 - a3 - b1 - b2 + b3 - p0 + p1 >= -4
    /\
    - a1 - a2 + a3 + b0 - b1 - b2 - b3 - p0 + p1 >= -5
    /\
    a0 - a1 - a2 - a3 + b0 - b2 + b3 - p0 + p1 >= -4
    /\
    - a0 + a1 - a2 + b0 - b1 - b2 + b3 - p0 + p1 >= -4
    /\
    a0 + a1 + a2 - b1 + b2 - b3 >= -1
;

% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 38
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s0_diff_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                       var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) =
    a0 - a1 - a2 + a3 - b0 >= -2
    /\
    a0 + a1 + a2 + a3 - b0 >= 0
    /\
    a0 - a1 + a2 + a3 + b0 >= 0
    /\
    a0 + a1 + a3 + b0 - b2 >= 0
    /\
    a0 + a1 + a3 + b0 - b3 >= 0
    /\
    - a0 - a2 - b0 + b2 - b3 >= -3
    /\
    - a2 - a3 - b0 - b2 + b3 >= -3
    /\
    a0 - a1 + a3 + b2 + b3 >= 0
    /\
    - a0 - a3 + b0 + b2 + b3 >= -1
    /\
    - a1 - b0 - b1 + b2 + b3 >= -2
    /\
    a1 + b0 - b1 + b2 + b3 >= 0
    /\
    a1 - a2 + b1 + b2 + b3 >= 0
    /\
    - a3 + b0 + b1 + b2 + b3 >= 0
    /\
    a0 - a1 - a3 + b0 - b1 - b2 >= -3
    /\
    a0 + a1 + a2 - a3 + b1 + b2 >= 0
    /\
    - a0 + a1 + a2 + a3 + b1 + b2 >= 0
    /\
    a0 + a1 - a3 - b0 - b1 - b3 >= -3
    /\
    - a0 - a1 + a3 + b0 - b1 - b3 >= -3
    /\
    a0 - a1 - a3 + b0 + b1 - b3 >= -2
    /\
    - a0 + a1 + a2 - a3 - b2 - b3 >= -3
    /\
    - a0 - a2 - a3 + b0 - b2 - b3 >= -4
    /\
    - a0 - a1 + a3 + b0 - b2 - b3 >= -3
    /\
    a0 - a1 - b0 - b1 - b2 - b3 >= -4
    /\
    - a0 + a3 - b0 - b1 - b2 - b3 >= -4
    /\
    - a0 - a3 - b0 + b1 - b2 - b3 >= -4
    /\
    a0 + a3 - b0 + b1 - b2 - b3 >= -2
    /\
    - a0 - a1 + a2 - a3 + b2 - b3 >= -3
    /\
    - a1 + a2 - a3 - b0 + b2 - b3 >= -3
    /\
    - a0 - a1 + a3 + b0 + b1 + b3 >= -1
    /\
    - a0 - a1 + a2 - a3 - b2 + b3 >= -3
    /\
    - a0 + a2 + a3 - b0 - b2 + b3 >= -2
    /\
    a0 - a2 - a3 - b1 - b2 + b3 >= -3
    /\
    a0 + a1 + a2 + b1 - b2 + b3 >= 0
    /\
    a2 + a3 + b0 + b1 - b2 + b3 >= 0
    /\
    - a0 + a1 + a2 - a3 + b2 + b3 >= -1
    /\
    a0 + a1 - a2 - a3 - b1 + b2 - b3 >= -3
    /\
    - a0 - a2 + a3 + b0 - b1 + b2 - b3 >= -3
    /\
    - a0 + a1 - a2 + a3 - b1 - b2 + b3 >= -3
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 39
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s1_diff_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                       var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) =
    - a0 - a1 + a2 + a3 - b1 >= -2
    /\
    a0 - a1 + a2 + a3 + b1 >= 0
    /\
    - a0 + a1 + a2 + a3 + b1 >= 0
    /\
    a0 + a1 + a2 + a3 - b2 >= 0
    /\
    a0 + a2 + a3 - b0 - b2 >= -1
    /\
    a2 + a3 + b0 - b1 + b2 >= 0
    /\
    - a0 - a2 + a3 - b0 - b3 >= -3
    /\
    a0 + a2 + a3 + b1 - b3 >= 0
    /\
    - a0 + a2 - a3 - b2 - b3 >= -3
    /\
    - a0 + b0 - b1 + b2 - b3 >= -2
    /\
    a0 + b0 + b1 + b2 - b3 >= 0
    /\
    a0 - a3 + b0 + b2 + b3 >= 0
    /\
    - a1 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 + a1 - a2 - a3 - b0 - b2 >= -4
    /\
    a0 + a1 - a2 + b0 - b1 - b2 >= -2
    /\
    - a1 - a2 - a3 + b0 - b1 - b2 >= -4
    /\
    - a0 - a1 - a2 + b0 + b1 - b2 >= -3
    /\
    - a0 + a1 - a2 - a3 + b0 + b2 >= -2
    /\
    - a0 - a1 - a2 - a3 - b1 + b2 >= -4
    /\
    a0 + a1 + a2 - b0 - b1 + b2 >= -1
    /\
    a0 - a2 - a3 - b0 - b1 + b2 >= -3
    /\
    - a0 - a1 - a3 - b0 + b1 + b2 >= -3
    /\
    a0 - a1 + a3 + b0 + b1 - b3 >= -1
    /\
    a0 - a2 - a3 - b0 - b2 - b3 >= -4
    /\
    a0 - a1 + a2 + b1 + b2 - b3 >= -1
    /\
    a0 + a1 - a2 + a3 + b0 + b3 >= 0
    /\
    - a0 + a2 - a3 - b0 - b1 + b3 >= -3
    /\
    a0 - a1 - a3 + b0 - b1 + b3 >= -2
    /\
    a1 + a2 - a3 + b0 + b1 + b3 >= 0
    /\
    - a0 - a2 + a3 + b0 - b2 + b3 >= -2
    /\
    a2 - a3 - b0 - b1 - b2 + b3 >= -3
    /\
    - a2 + a3 - b0 - b1 - b2 + b3 >= -3
    /\
    - a2 - a3 - b0 + b1 - b2 + b3 >= -3
    /\
    a2 + a3 - b0 + b1 - b2 + b3 >= -1
    /\
    a1 - a2 + a3 + b1 + b2 + b3 >= 0
    /\
    a1 + a2 - b0 + b1 + b2 + b3 >= 0
    /\
    a0 - a1 + a2 - a3 + b0 - b2 - b3 >= -3
    /\
    a0 - a1 - a2 + a3 - b0 + b2 - b3 >= -3
    /\
    a0 - a1 - a2 + a3 - b1 + b2 + b3 >= -2
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 43
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s2_diff_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                       var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) =
    - a0 + a1 - a2 + a3 - b3 >= -2
    /\
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    - a1 + a2 + a3 + b0 - b3 >= -1
    /\
    - a0 - a2 + a3 + b1 - b3 >= -2
    /\
    - a3 - b0 - b1 + b2 - b3 >= -3
    /\
    a3 + b0 + b1 + b2 - b3 >= 0
    /\
    - a0 + a1 + a2 + a3 + b3 >= 0
    /\
    a0 + a1 + a3 - b2 + b3 >= 0
    /\
    a1 + a2 + b0 - b2 + b3 >= 0
    /\
    - a1 + b0 + b1 + b2 + b3 >= 0
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    a0 - a1 - a2 - a3 + b0 + b1 >= -2
    /\
    a0 + a1 + a2 - a3 + b0 + b1 >= 0
    /\
    a1 - a2 - a3 - b0 + b1 - b2 >= -3
    /\
    - a1 + a2 - a3 - b0 + b1 - b2 >= -3
    /\
    - a1 - a2 + a3 - b0 + b1 - b2 >= -3
    /\
    - a0 + a1 + a2 - a3 + b0 + b2 >= -1
    /\
    a0 - a1 - a2 - a3 + b1 + b2 >= -2
    /\
    a0 + a1 + a2 - b0 + b1 + b2 >= 0
    /\
    a1 - a2 - a3 + b0 - b2 - b3 >= -3
    /\
    a0 + a3 + b0 - b1 - b2 - b3 >= -2
    /\
    a1 + a3 - b0 + b1 - b2 - b3 >= -2
    /\
    - a0 - a1 - a2 + b0 + b2 - b3 >= -3
    /\
    a1 - a2 - a3 - b1 + b2 - b3 >= -3
    /\
    a0 - a2 - b0 - b1 + b2 - b3 >= -3
    /\
    - a1 + a2 - b0 + b1 + b2 - b3 >= -2
    /\
    a0 + a2 + a3 - b0 - b1 + b3 >= -1
    /\
    a0 - a1 + a2 + b0 + b1 + b3 >= 0
    /\
    - a0 + a1 - a3 + b0 + b1 + b3 >= -1
    /\
    - a1 - a2 + b0 - b1 - b2 + b3 >= -3
    /\
    a0 - a1 + a2 + a3 + b2 + b3 >= 0
    /\
    a0 + a1 + b0 - b1 + b2 + b3 >= 0
    /\
    a0 + a1 - b0 + b1 + b2 + b3 >= 0
    /\
    - a0 + a1 + a2 - a3 - b0 - b1 - b2 >= -4
    /\
    - a0 - a1 - a2 - b0 - b1 - b2 - b3 >= -6
    /\
    a0 - a1 - a3 - b0 - b1 - b2 + b3 >= -4
    /\
    - a0 - a2 - a3 - b0 - b1 - b2 + b3 >= -5
    /\
    - a0 + a2 + a3 + b0 - b1 - b2 + b3 >= -2
    /\
    - a0 - a1 + a3 - b0 + b1 - b2 + b3 >= -3
    /\
    - a0 - a1 - a3 + b0 - b1 + b2 + b3 >= -3
    /\
    - a0 - a1 - a3 - b0 + b1 + b2 + b3 >= -3
    /\
    a3 - b0 - b1 + b2 + b3 >= -1
    /\
    - a1 + a2 - b1 - b2 - b3 >= -3
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 48
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s3_diff_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                       var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) =
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    a0 + a1 + a3 - b1 + b3 >= 0
    /\
    - a0 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 - a1 + a2 - a3 - b0 + b1 >= -3
    /\
    a0 - a2 - a3 - b0 + b1 - b2 >= -3
    /\
    - a1 + a2 + a3 + b0 + b1 - b2 >= -1
    /\
    a0 - a2 - a3 - b0 - b1 + b2 >= -3
    /\
    a0 - a1 + a3 + b0 - b1 + b2 >= -1
    /\
    a0 + a1 + a2 - a3 + b1 + b2 >= 0
    /\
    a0 - a1 + a2 + a3 + b1 + b2 >= 0
    /\
    a0 - a2 - a3 + b0 + b1 + b2 >= -1
    /\
    - a0 + a1 + a2 + b0 - b1 - b3 >= -2
    /\
    - a0 - a1 + a2 - b0 + b1 - b3 >= -3
    /\
    a1 - a2 + a3 + b0 + b1 - b3 >= -1
    /\
    - a0 + a1 - a2 - a3 - b2 - b3 >= -4
    /\
    - a0 + a2 + a3 - b0 - b2 - b3 >= -3
    /\
    - a1 - a3 - b0 - b1 - b2 - b3 >= -5
    /\
    a2 + a3 - b0 - b1 - b2 - b3 >= -3
    /\
    a0 - a3 + b0 - b1 - b2 - b3 >= -3
    /\
    - a0 - a1 + b0 + b1 - b2 - b3 >= -3
    /\
    - a0 - a1 - a2 - a3 + b2 - b3 >= -4
    /\
    - a0 + a1 + a3 + b0 + b2 - b3 >= -1
    /\
    a0 - a2 - b0 - b1 + b2 - b3 >= -3
    /\
    a0 - a1 + a2 - a3 + b0 + b3 >= -1
    /\
    - a0 + a2 - a3 - b0 - b1 + b3 >= -3
    /\
    a0 - a2 + a3 + b0 + b1 + b3 >= 0
    /\
    a0 + a1 - a3 - b0 - b2 + b3 >= -2
    /\
    a1 - a2 - a3 - b1 - b2 + b3 >= -3
    /\
    - a1 - a3 + b0 - b1 - b2 + b3 >= -3
    /\
    a1 + a2 + a3 + b1 - b2 + b3 >= 0
    /\
    a0 + a2 + b0 + b1 - b2 + b3 >= 0
    /\
    - a0 - a1 - a2 + b0 + b2 + b3 >= -2
    /\
    a1 + a2 + b0 - b1 + b2 + b3 >= 0
    /\
    a1 + a3 - b0 + b1 + b2 + b3 >= 0
    /\
    - a0 - a1 - a2 + a3 + b0 - b1 - b2 >= -4
    /\
    - a0 - a1 - a2 + a3 - b0 + b1 + b2 >= -3
    /\
    - a0 + a1 - a2 + a3 - b0 - b1 - b3 >= -4
    /\
    a0 - a1 - a2 + a3 - b0 + b1 - b3 >= -3
    /\
    - a0 + a1 + a2 - a3 - b0 + b2 - b3 >= -3
    /\
    a0 + a1 - a2 - a3 - b1 + b2 - b3 >= -3
    /\
    a0 - a1 + a2 - a3 - b1 + b2 - b3 >= -3
    /\
    - a0 - a1 + a2 + a3 - b1 + b2 - b3 >= -3
    /\
    - a0 + a1 - a2 + a3 - b0 - b2 + b3 >= -3
    /\
    a0 - a1 + a3 - b0 - b1 - b2 + b3 >= -3
    /\
    - a1 - a2 + a3 - b0 + b1 - b2 + b3 >= -3
    /\
    - a0 - a1 + a3 - b0 - b1 + b2 + b3 >= -3
    /\
    a0 - a1 - a2 - a3 + b1 + b2 + b3 >= -2
    /\
    - a0 + a1 - a2 + b1 + b2 + b3 >= -1
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 44
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s4_diff_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                       var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) =
    a0 + a1 + a2 + a3 - b2 >= 0
    /\
    - a1 - a2 - a3 + b2 - b3 >= -3
    /\
    - a0 + a1 - a2 - a3 + b3 >= -2
    /\
    a0 - a1 + a2 + a3 + b3 >= 0
    /\
    - a0 + a2 + a3 + b0 + b3 >= 0
    /\
    a0 + a2 + a3 - b1 + b3 >= 0
    /\
    a0 + a2 + b0 - b2 + b3 >= 0
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    - a3 + b0 + b1 + b2 + b3 >= 0
    /\
    a0 + a1 - a2 - a3 - b0 - b1 >= -3
    /\
    a0 - a1 - a2 + a3 - b1 - b2 >= -3
    /\
    - a0 - a2 + a3 - b0 - b1 - b2 >= -4
    /\
    a1 + a2 + a3 + b0 - b1 - b2 >= -1
    /\
    - a0 - a1 + a2 + a3 + b1 - b2 >= -2
    /\
    a0 + a1 + a3 - b0 + b1 - b2 >= -1
    /\
    - a0 - a2 - a3 + b0 + b1 - b2 >= -3
    /\
    a0 + a2 - a3 + b0 - b1 + b2 >= -1
    /\
    a0 + a1 + a3 + b0 - b1 + b2 >= 0
    /\
    a0 - a1 + a3 + b0 + b1 + b2 >= 0
    /\
    a1 + a2 + a3 - b0 - b1 - b3 >= -2
    /\
    - a0 - a1 - a2 + b0 - b1 - b3 >= -4
    /\
    - a0 - a1 - a2 - b1 - b2 - b3 >= -5
    /\
    a1 + a2 - b0 - b1 - b2 - b3 >= -3
    /\
    - a0 - a1 + b0 - b1 - b2 - b3 >= -4
    /\
    - a0 + a2 - a3 + b1 - b2 - b3 >= -3
    /\
    - a1 - a2 - b0 + b1 - b2 - b3 >= -4
    /\
    a0 + a1 - a2 - b1 + b2 - b3 >= -2
    /\
    a0 + a1 + a2 + b1 + b2 - b3 >= 0
    /\
    a1 + a3 + b0 + b1 + b2 - b3 >= 0
    /\
    a0 - a1 - a2 - b0 - b1 + b3 >= -3
    /\
    a0 + a1 + a3 - b0 + b1 + b3 >= 0
    /\
    a0 - a3 + b0 - b1 - b2 + b3 >= -2
    /\
    a0 + a3 + b0 + b1 - b2 + b3 >= 0
    /\
    - a0 - a2 + a3 - b1 + b2 + b3 >= -2
    /\
    a1 + a2 - b0 + b1 + b2 + b3 >= 0
    /\
    a0 - a1 + a2 - a3 - b0 + b1 - b2 >= -3
    /\
    - a0 + a1 - a2 + a3 + b0 - b2 - b3 >= -3
    /\
    a0 - a1 + a2 - a3 - b1 - b2 - b3 >= -4
    /\
    - a0 - a1 + a2 - b0 - b1 + b2 - b3 >= -4
    /\
    - a0 + a1 - a2 - b0 + b1 + b2 - b3 >= -3
    /\
    - a1 + a2 + a3 - b0 + b1 + b2 - b3 >= -2
    /\
    a1 + a2 - a3 + b0 - b1 + b2 - b3 >= -2
    /\
    a1 - a3 + b0 + b1 - b2 - b3 >= -2
    /\
    - a0 + a2 - a3 - b0 + b3 >= -2
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 43
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s5_diff_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                       var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) =
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    - a0 - a1 - a2 + b0 - b3 >= -3
    /\
    - a0 + a1 - a2 - a3 + b3 >= -2
    /\
    a0 + a2 + a3 - b2 + b3 >= 0
    /\
    a0 + a2 - b0 + b2 + b3 >= 0
    /\
    a0 + a3 - b0 + b2 + b3 >= 0
    /\
    a2 - a3 + b0 + b2 + b3 >= 0
    /\
    a0 + a3 - b1 + b2 + b3 >= 0
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 - a2 - a3 - b0 + b1 + b2 >= -3
    /\
    - a1 + a2 + a3 + b0 + b1 + b2 >= 0
    /\
    - a0 - a1 + a2 - b0 - b1 - b3 >= -4
    /\
    a1 + a2 - a3 - b0 - b1 - b3 >= -3
    /\
    - a0 + a1 + a3 + b0 - b1 - b3 >= -2
    /\
    a0 - a1 - a3 - b0 + b1 - b3 >= -3
    /\
    a0 + a1 + a3 + b0 + b1 - b3 >= 0
    /\
    - a0 - a1 - a3 - b0 - b2 - b3 >= -5
    /\
    - a1 - a2 - b0 + b1 - b2 - b3 >= -4
    /\
    - a0 - a1 - a3 - b1 + b2 - b3 >= -4
    /\
    - a1 - a3 - b0 - b1 + b2 - b3 >= -4
    /\
    a1 + a2 + b0 + b1 + b2 - b3 >= 0
    /\
    - a0 + a1 + a2 + a3 + b1 + b3 >= 0
    /\
    - a0 + a2 - a3 - b1 - b2 + b3 >= -3
    /\
    - a0 - a2 + a3 - b1 - b2 + b3 >= -3
    /\
    - a1 + a2 + a3 - b1 - b2 + b3 >= -2
    /\
    a0 - a1 - b0 + b1 - b2 + b3 >= -2
    /\
    a0 + a1 + b0 + b1 - b2 + b3 >= 0
    /\
    - a2 + a3 - b0 - b1 + b2 + b3 >= -2
    /\
    a1 + a2 + b0 - b1 + b2 + b3 >= 0
    /\
    a0 + a1 - a2 - a3 - b0 - b2 - b3 >= -4
    /\
    a0 - a1 - a2 + a3 - b0 - b2 - b3 >= -4
    /\
    a1 - a2 + a3 - b0 - b1 - b2 - b3 >= -4
    /\
    a0 - a1 - a3 + b0 - b1 - b2 - b3 >= -4
    /\
    a0 + a2 - a3 + b0 - b1 - b2 - b3 >= -3
    /\
    - a0 + a2 + a3 - b0 + b1 - b2 - b3 >= -3
    /\
    a0 - a1 + a2 + b0 + b1 - b2 - b3 >= -2
    /\
    - a0 + a1 - a3 + b0 + b1 - b2 - b3 >= -3
    /\
    a0 - a1 - a2 + a3 - b1 + b2 - b3 >= -3
    /\
    - a0 + a1 - a3 - b0 + b1 + b2 - b3 >= -3
    /\
    - a0 - a1 + a2 + a3 - b0 + b2 + b3 >= -2
    /\
    a0 + a1 - a2 + b0 - b1 + b2 - b3 >= -2
    /\
    a0 - a2 - a3 - b2 + b3 >= -2
    /\
    a1 - a2 - b0 + b1 + b2 - b3 >= -2
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 38
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s6_diff_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                       var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) =
    a0 + a1 + a2 + a3 - b0 >= 0
    /\
    a0 - a1 - a2 + a3 - b2 >= -2
    /\
    - a0 + a3 - b0 - b1 - b2 >= -3
    /\
    a0 - a1 + a2 + a3 + b2 >= 0
    /\
    a0 - a3 - b0 - b2 - b3 >= -3
    /\
    a0 + a3 - b1 + b2 - b3 >= -1
    /\
    a0 - a2 + a3 + b1 + b3 >= 0
    /\
    - a2 - b0 + b1 - b2 + b3 >= -2
    /\
    a2 + b0 + b1 - b2 + b3 >= 0
    /\
    - a0 + a2 + b1 + b2 + b3 >= 0
    /\
    a2 - a3 + b1 + b2 + b3 >= 0
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 - a1 - a2 - a3 + b1 - b2 >= -4
    /\
    - a1 - a2 - a3 + b0 + b1 + b2 >= -2
    /\
    a0 - a2 - a3 + b0 - b1 - b3 >= -3
    /\
    - a0 + a2 - a3 + b0 - b1 - b3 >= -3
    /\
    a0 + a2 + a3 + b0 - b1 - b3 >= -1
    /\
    a0 + a1 + a2 + b0 + b1 - b3 >= 0
    /\
    a1 + a2 + a3 + b0 + b1 - b3 >= 0
    /\
    - a0 - a2 + a3 + b0 - b2 - b3 >= -3
    /\
    - a0 + a1 - a3 - b1 - b2 - b3 >= -4
    /\
    - a0 - a1 + a2 + b1 - b2 - b3 >= -3
    /\
    - a0 - a3 - b0 - b1 + b2 - b3 >= -4
    /\
    - a2 + a3 + b0 - b1 + b2 - b3 >= -2
    /\
    - a0 - a2 - a3 + b1 + b2 - b3 >= -3
    /\
    a1 - a2 + a3 + b1 + b2 - b3 >= -1
    /\
    a0 + a1 - a3 + b0 - b1 + b3 >= -1
    /\
    a1 + a2 + a3 + b0 - b1 + b3 >= 0
    /\
    - a0 - a1 - b0 - b1 - b2 + b3 >= -4
    /\
    - a1 - a3 + b0 - b1 - b2 + b3 >= -3
    /\
    - a0 + a1 - a3 + b1 - b2 + b3 >= -2
    /\
    - a0 - a1 - a2 + b0 + b2 + b3 >= -2
    /\
    a0 + a1 - a2 - b1 + b2 + b3 >= -1
    /\
    - a0 - a2 - a3 - b1 + b2 + b3 >= -3
    /\
    - a0 - a1 + a2 + a3 - b0 + b1 - b3 >= -3
    /\
    a0 - a1 - a3 - b0 + b1 + b2 - b3 >= -3
    /\
    a0 - a1 + a2 - a3 - b0 - b1 + b3 >= -3
    /\
    - a0 - a1 + a3 - b0 - b1 + b2 + b3 >= -3
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 49
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s7_diff_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                       var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) =
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    - a3 + b0 + b1 + b2 + b3 >= 0
    /\
    a0 - a1 - a2 + a3 - b0 - b2 >= -3
    /\
    - a0 + a1 + a2 - a3 + b0 - b2 >= -2
    /\
    - a0 - a2 - a3 - b0 - b1 - b2 >= -5
    /\
    - a0 - a1 - a2 - a3 + b1 - b2 >= -4
    /\
    - a1 + a2 + a3 + b0 + b1 - b2 >= -1
    /\
    - a0 - a1 - a2 + a3 - b0 + b2 >= -3
    /\
    - a1 - a2 + a3 - b0 - b1 + b2 >= -3
    /\
    a0 + a1 - a2 + a3 + b1 + b2 >= 0
    /\
    - a0 - a2 + a3 + b0 + b1 + b2 >= -1
    /\
    - a0 - a1 + a2 + a3 - b1 - b3 >= -3
    /\
    a0 + a2 - a3 + b0 - b1 - b3 >= -2
    /\
    a0 - a2 - a3 - b0 + b1 - b3 >= -3
    /\
    - a0 - a1 - a3 + b0 + b1 - b3 >= -3
    /\
    a0 - a2 + a3 - b1 - b2 - b3 >= -3
    /\
    - a0 + a1 - b0 - b1 - b2 - b3 >= -4
    /\
    - a0 - a1 + a2 + b1 - b2 - b3 >= -3
    /\
    a2 - a3 - b0 + b1 - b2 - b3 >= -3
    /\
    a1 + a2 + a3 + b0 + b2 - b3 >= 0
    /\
    - a1 - a2 - a3 - b1 + b2 - b3 >= -4
    /\
    a0 + a2 + b0 + b1 + b2 - b3 >= 0
    /\
    a0 - a1 - a2 - a3 + b0 + b3 >= -2
    /\
    - a0 + a1 - a2 - a3 - b1 + b3 >= -3
    /\
    a1 + a2 + a3 - b0 - b1 + b3 >= -1
    /\
    a0 + a2 + a3 - b0 + b1 + b3 >= 0
    /\
    - a1 - a2 + a3 - b0 - b2 + b3 >= -3
    /\
    a1 - a2 + a3 + b0 - b2 + b3 >= -1
    /\
    a0 + a2 + a3 + b0 - b2 + b3 >= 0
    /\
    - a1 + a2 - a3 - b1 - b2 + b3 >= -3
    /\
    - a0 - a2 - b0 + b1 - b2 + b3 >= -3
    /\
    - a2 + a3 - b0 + b1 - b2 + b3 >= -2
    /\
    a1 + a2 + b0 + b1 - b2 + b3 >= 0
    /\
    a0 - a1 - a3 - b0 + b2 + b3 >= -2
    /\
    a0 - a1 + a3 + b0 + b2 + b3 >= 0
    /\
    - a0 + a2 + a3 + b0 + b2 + b3 >= 0
    /\
    a1 - a2 - a3 - b1 + b2 + b3 >= -2
    /\
    a0 + a1 + b0 - b1 + b2 + b3 >= 0
    /\
    a1 + a2 - b0 + b1 + b2 + b3 >= 0
    /\
    - a0 - a1 + a2 - a3 - b0 - b1 + b2 >= -4
    /\
    - a0 + a1 + a2 - a3 - b0 + b1 - b3 >= -3
    /\
    a0 + a1 - a2 - a3 + b0 - b2 - b3 >= -3
    /\
    - a1 - a2 + a3 + b0 - b1 - b2 - b3 >= -4
    /\
    - a0 + a1 - a2 + b0 + b1 - b2 - b3 >= -3
    /\
    - a0 + a1 - a2 + b0 - b1 + b2 - b3 >= -3
    /\
    - a0 + a1 - a2 - b0 + b1 + b2 - b3 >= -3
    /\
    a0 + a1 - a3 - b0 - b1 - b2 + b3 >= -3
    /\
    - a0 + a1 - a2 - b0 - b1 + b2 + b3 >= -3
    /\
    a0 + a1 + a2 - b1 + b2 - b3 >= -1
;
% #############################################################################################################################################
% We use S-box Analyzer to encode deterministic differential behaviour:
% x0: lsb
predicate s0_diff_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                    var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                    var 0..1: d) =
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 1 /\ x0 == 0) then (y3 == 1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y3 == 1 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 0 /\ x2 == 1 /\ x1 == 1 /\ x0 == 0) then (y3 == 0 /\ y2 == -1 /\ y1 == -1 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;

% #############################################################################################################################################
predicate s1_diff_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                    var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                    var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y3 == -1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s2_diff_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                    var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                    var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 1 /\ x0 == 0) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 1 /\ x0 == 0) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 0)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s3_diff_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                    var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                    var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s4_diff_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                    var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                    var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 1 /\ x0 == 1) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x3 == 1 /\ x2 == 1 /\ x1 == 1 /\ x0 == 1) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 0)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s5_diff_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                    var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                    var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 1 /\ x0 == 1) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 1)
    elseif (x3 == 1 /\ x2 == 1 /\ x1 == 1 /\ x0 == 1) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 0)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s6_diff_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                    var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                    var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 1 /\ x0 == 0) then (y3 == -1 /\ y2 == -1 /\ y1 == 1 /\ y0 == -1)
    elseif (x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y3 == -1 /\ y2 == -1 /\ y1 == 1 /\ y0 == -1)
    elseif (x3 == 0 /\ x2 == 1 /\ x1 == 1 /\ x0 == 0) then (y3 == -1 /\ y2 == -1 /\ y1 == 0 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s7_diff_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                    var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                    var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
% We use S-box Analyzer to encode Squared-LAT: 
% https://github.com//sboxanalyzer

/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 49
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 4.0000 p0 + 2.0000 p1
*/
predicate s0_lin(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                 var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                 var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - a1 + p0 + p1 >= 0
    /\
    b1 + b2 + b3 - p0 >= 0
    /\
    - a2 - b0 - b1 + p0 >= -2
    /\
    - a0 + a2 - a3 + p1 >= -1
    /\
    b1 - b2 - b3 + p1 >= -1
    /\
    - a0 + a1 + a2 - a3 + b1 >= -1
    /\
    - a0 - a3 + b0 - b2 - b3 >= -3
    /\
    - a1 - b0 + b1 - b2 - b3 >= -3
    /\
    a3 + b0 - b1 + b2 + p0 >= 0
    /\
    - a2 + a3 + b0 + b3 + p0 >= 0
    /\
    a1 + a3 - b1 + b3 + p0 >= 0
    /\
    a0 + b0 - b1 + b3 + p0 >= 0
    /\
    a0 + a2 + a3 + b0 - p1 >= 0
    /\
    a0 + a2 + a3 + b1 - p1 >= 0
    /\
    - a1 - a3 - b1 - b2 - p1 >= -4
    /\
    - a1 - a2 + b2 - b3 - p1 >= -3
    /\
    b0 + b1 + b2 + b3 - p1 >= 0
    /\
    a0 - a2 - a3 + b1 + p1 >= -1
    /\
    - a0 - a2 + a3 + b1 + p1 >= -1
    /\
    a0 + a3 - b2 - b3 + p1 >= -1
    /\
    a2 - b1 + b2 - b3 + p1 >= -1
    /\
    a2 - b1 - b2 + b3 + p1 >= -1
    /\
    - a0 - a3 + b2 + b3 + p1 >= -1
    /\
    a0 + a2 + a3 - p0 + p1 >= 0
    /\
    - a0 - a2 - a3 - b1 - b2 - b3 >= -5
    /\
    a1 + a2 + b0 + b1 - b2 + p0 >= 0
    /\
    a0 + a3 + b0 + b1 - b2 + p0 >= 0
    /\
    a0 + a1 + a2 - b0 + b2 + p0 >= 0
    /\
    - a0 + a2 + a3 + b1 + b2 + p0 >= 0
    /\
    - a0 + a1 - b0 + b1 + b2 + p0 >= -1
    /\
    a1 - a2 - a3 + b0 - b3 + p0 >= -2
    /\
    a0 + a1 + a2 - b1 - b3 + p0 >= -1
    /\
    - a0 + a3 - b1 - b2 - b3 + p0 >= -3
    /\
    a0 + a2 + b1 + b2 - b3 + p0 >= 0
    /\
    a0 + a3 + b1 + b2 - b3 + p0 >= 0
    /\
    a0 + a1 - b0 + b1 + b3 + p0 >= 0
    /\
    - a0 - a2 - a3 - b2 + b3 + p0 >= -3
    /\
    a2 - a3 + b0 + b2 + b3 + p0 >= 0
    /\
    - a0 - a1 + a2 + a3 - b0 - p1 >= -3
    /\
    a0 - a2 - a3 - b0 - b2 - p1 >= -4
    /\
    - a0 - a1 + a2 - b1 + b2 - p1 >= -3
    /\
    - a3 - b0 - b1 + b2 - b3 - p1 >= -4
    /\
    a0 - a1 + b1 - b2 + b3 - p1 >= -2
    /\
    - a1 + a3 + b1 - b2 + b3 - p1 >= -2
    /\
    - a1 + a2 - a3 + b2 + b3 - p1 >= -2
    /\
    a0 - a2 + a3 + b2 + b3 + p1 >= 0
    /\
    - a0 + a1 - a2 + a3 - b2 - b3 + p0 >= -3
    /\
    - a0 + a1 + a2 - b0 - b2 + b3 + p0 >= -2
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 49
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 4.0000 p0 + 2.0000 p1
*/
predicate s1_lin(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                 var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                 var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - a0 + p0 + p1 >= 0
    /\
    - a0 - a1 - b3 - p1 >= -3
    /\
    a1 - a2 - a3 + p1 >= -1
    /\
    - b0 - b2 + b3 + p1 >= -1
    /\
    - a2 - a3 - b0 - b1 - b2 >= -4
    /\
    a1 + a2 + a3 - b0 + b3 >= 0
    /\
    a1 + a2 + a3 - b2 + b3 >= 0
    /\
    a1 - b0 + b1 - b2 + b3 >= -1
    /\
    - a1 + a2 + b0 - b3 + p0 >= -1
    /\
    - a2 + a3 - b1 - b3 + p0 >= -2
    /\
    - a1 - b0 + b2 - b3 + p0 >= -2
    /\
    a1 + a2 + a3 + b1 - p1 >= 0
    /\
    a2 + b1 + b2 + b3 - p1 >= 0
    /\
    b0 + b1 + b2 + b3 - p1 >= 0
    /\
    a2 + a3 - b0 - b2 + p1 >= -1
    /\
    - a2 - a3 + b0 + b2 + p1 >= -1
    /\
    a1 + b0 - b2 - b3 + p1 >= -1
    /\
    a1 - b0 + b2 - b3 + p1 >= -1
    /\
    - a1 + a2 - a3 + b3 + p1 >= -1
    /\
    - a1 - a2 + a3 + b3 + p1 >= -1
    /\
    b0 + b2 + b3 - p0 + p1 >= 0
    /\
    - a1 - a2 - a3 - b0 - b2 - b3 >= -5
    /\
    a0 - a1 - a2 - a3 - b2 + p0 >= -3
    /\
    a0 + a1 + a2 - b0 + b2 + p0 >= 0
    /\
    a0 + a1 - a3 + b0 + b2 + p0 >= 0
    /\
    a0 + a1 + a2 - b1 - b3 + p0 >= -1
    /\
    a2 - a3 - b0 - b2 - b3 + p0 >= -3
    /\
    a0 - a2 + b1 - b2 - b3 + p0 >= -2
    /\
    a3 + b0 + b1 + b2 - b3 + p0 >= 0
    /\
    a0 - a1 + a2 + a3 + b3 + p0 >= 0
    /\
    a0 + a2 - a3 + b0 + b3 + p0 >= 0
    /\
    a1 + a3 + b0 - b1 + b3 + p0 >= 0
    /\
    a0 + a1 - a2 + b1 + b3 + p0 >= 0
    /\
    a3 + b0 + b1 - b2 + b3 + p0 >= 0
    /\
    a0 - a1 - a2 + b2 + b3 + p0 >= -1
    /\
    a3 - b0 - b1 + b2 + b3 + p0 >= -1
    /\
    - a0 + a1 + a3 - b0 + b1 - p1 >= -2
    /\
    - a0 + a1 + a2 + b0 + b1 - p1 >= -1
    /\
    - a0 - a1 - a3 - b0 + b2 - p1 >= -4
    /\
    - a0 + a1 - a3 - b1 + b2 - p1 >= -3
    /\
    - a0 - a3 + b0 - b2 - b3 - p1 >= -4
    /\
    - a1 - a2 + b0 + b1 + b3 - p1 >= -2
    /\
    - a0 + b0 - b1 - b2 + b3 - p1 >= -3
    /\
    a2 + a3 + b0 + b2 - b3 + p1 >= 0
    /\
    - a2 - a3 - b0 + b1 + b2 - b3 - p1 >= -4
    /\
    - a0 + a2 - a3 - b0 + b1 + b3 - p1 >= -3
    /\
    - a0 - a2 + a3 + b1 - b2 + b3 - p1 >= -3
    /\
    - a1 - b0 - b1 - b2 + b3 + p0 - p1 >= -4
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.02 seconds
Number of constraints: 44
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 4.0000 p0 + 2.0000 p1
*/
predicate s2_lin(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                 var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                 var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - a3 + p0 + p1 >= 0
    /\
    a0 - a1 - a2 + p1 >= -1
    /\
    a0 - b1 - b2 + p1 >= -1
    /\
    a0 + a1 + a2 - b0 - b1 >= -1
    /\
    - a2 - b0 - b1 - b2 - b3 >= -4
    /\
    a0 - a1 - a2 + a3 + b3 >= -1
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 - a1 - a3 - b0 + p0 >= -3
    /\
    - a0 - a2 + b0 - b2 + p0 >= -2
    /\
    a1 - a2 + b1 - b2 + p0 >= -1
    /\
    - a1 + a2 + b0 + b2 + p0 >= 0
    /\
    a1 + a3 - b1 + b2 + p0 >= 0
    /\
    a1 - b0 - b1 + b2 + p0 >= -1
    /\
    a3 - b0 + b1 + b2 + p0 >= 0
    /\
    - a0 - b0 + b1 - b3 + p0 >= -2
    /\
    - a0 + a3 + b0 + b3 + p0 >= 0
    /\
    a1 - a3 + b1 + b3 + p0 >= 0
    /\
    a0 + a1 + a2 + a3 - p1 >= 0
    /\
    a1 + a2 - b0 - b2 + p1 >= -1
    /\
    a0 + b1 + b2 - p0 + p1 >= 0
    /\
    a0 + a1 + a2 + b0 + b1 - b2 >= 0
    /\
    - a0 - a1 - a2 + a3 - b1 + p0 >= -3
    /\
    a1 - a2 - a3 + b0 - b1 + p0 >= -2
    /\
    - a0 - a1 + a2 + a3 - b2 + p0 >= -2
    /\
    a0 - a1 + a2 + b1 - b2 + p0 >= -1
    /\
    a0 + a2 + a3 - b1 + b2 + p0 >= 0
    /\
    - a0 + a1 + a2 + b0 - b3 + p0 >= -1
    /\
    - a1 + b0 - b1 - b2 - b3 + p0 >= -3
    /\
    - a1 - a3 - b1 + b2 - b3 + p0 >= -3
    /\
    a0 + b0 + b1 + b2 - b3 + p0 >= 0
    /\
    a0 + a2 - a3 - b1 + b3 + p0 >= -1
    /\
    a1 + a3 + b0 - b1 + b3 + p0 >= 0
    /\
    - a0 + a2 - b1 - b2 + b3 + p0 >= -2
    /\
    a0 + b0 - b1 + b2 + b3 + p0 >= 0
    /\
    a0 - a1 - a3 + b1 - b3 - p1 >= -3
    /\
    a0 - a1 - b0 - b2 + b3 - p1 >= -3
    /\
    - a0 + a1 - a2 - b0 - b1 + p1 >= -3
    /\
    - a0 - a1 + a2 - b0 - b1 + p1 >= -3
    /\
    - a0 + a1 - a2 + b0 + b1 + p1 >= -1
    /\
    - a0 - a1 + a2 + b0 + b1 + p1 >= -1
    /\
    - a0 - a1 - a2 - b0 - b2 + p1 >= -4
    /\
    - a0 - a1 - a2 + b0 + b2 + p1 >= -2
    /\
    a1 + a2 + b0 + b2 - p0 + p1 >= 0
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.02 seconds
Number of constraints: 72
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 4.0000 p0 + 2.0000 p1
*/
predicate s3_lin(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                 var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                 var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - b3 + p0 + p1 >= 0
    /\
    - a0 - a1 - a2 - a3 - b0 >= -4
    /\
    - a2 - a3 - b0 - b1 - b2 >= -4
    /\
    - a1 - b0 - b1 - b2 + b3 >= -3
    /\
    a0 + a1 + a2 + a3 - p0 >= 0
    /\
    a1 - a2 - a3 - b1 + p0 >= -2
    /\
    - a1 + a2 + b0 - b2 + p0 >= -1
    /\
    a0 + b0 - b1 + b2 + p0 >= 0
    /\
    - a0 + a1 + b0 + b3 + p0 >= 0
    /\
    a2 + a3 - b1 + b3 + p0 >= 0
    /\
    a0 - b0 + b1 + b3 + p0 >= 0
    /\
    a1 + b1 - b2 + b3 + p0 >= 0
    /\
    a2 - b0 + b2 + b3 + p0 >= 0
    /\
    a0 + a1 + a2 - b3 - p1 >= -1
    /\
    - a0 - a3 + b0 - b1 - b2 - b3 >= -4
    /\
    - a0 - a1 - a3 - b0 + b2 - b3 >= -4
    /\
    a1 + a2 + a3 + b0 + b1 - p0 >= 0
    /\
    a1 + a2 + b0 + b1 + b3 - p0 >= 0
    /\
    a0 + a2 + b1 + b2 + b3 - p0 >= 0
    /\
    - a0 - a1 - a2 + a3 - b1 + p0 >= -3
    /\
    a0 - a1 - a2 - b0 + b1 + p0 >= -2
    /\
    - a0 + a1 - a3 + b0 + b1 + p0 >= -1
    /\
    a0 - a1 - a3 - b1 - b2 + p0 >= -3
    /\
    - a2 + a3 - b0 + b1 - b2 + p0 >= -2
    /\
    a0 - a3 - b0 + b1 + b2 + p0 >= -1
    /\
    - a0 - a3 + b0 + b1 + b3 + p0 >= -1
    /\
    a0 + a3 + b0 - b2 + b3 + p0 >= 0
    /\
    a1 - a2 + a3 + b2 + b3 + p0 >= 0
    /\
    a0 - a3 + b1 + b2 + b3 + p0 >= 0
    /\
    - a0 + a1 + a3 + b0 - b1 - p1 >= -2
    /\
    a2 + a3 - b0 - b1 + b2 - p1 >= -2
    /\
    - a0 + a2 - a3 - b0 - b3 - p1 >= -4
    /\
    - a1 + a2 - a3 + b0 - b3 - p1 >= -3
    /\
    a0 - a2 - b0 - b1 - b3 - p1 >= -4
    /\
    - a0 + a3 - b0 + b1 - b3 - p1 >= -3
    /\
    - a0 - a2 - b0 - b2 - b3 - p1 >= -5
    /\
    - a1 + a3 + b1 - b2 - b3 - p1 >= -3
    /\
    a0 + b0 + b1 - b2 - b3 - p1 >= -2
    /\
    a3 + b0 + b1 + b2 + p0 - p1 >= 0
    /\
    a1 + a2 + a3 - b0 - b1 + p1 >= -1
    /\
    a0 + a2 + a3 - b1 - b2 + p1 >= -1
    /\
    - a0 - a2 - a3 + b1 + b2 + p1 >= -2
    /\
    - a1 - a2 - b0 - b1 + b3 + p1 >= -3
    /\
    a1 + a2 - b0 - b1 + b3 + p1 >= -1
    /\
    - a1 - a2 + b0 + b1 + b3 + p1 >= -1
    /\
    - a0 - a1 - b0 - b2 + b3 + p1 >= -3
    /\
    a0 + a1 - b0 - b2 + b3 + p1 >= -1
    /\
    - a0 - a2 - b1 - b2 + b3 + p1 >= -3
    /\
    a0 + a2 - b1 - b2 + b3 + p1 >= -1
    /\
    - a0 - a1 + b0 + b2 + b3 + p1 >= -1
    /\
    - a0 - a2 + b1 + b2 + b3 + p1 >= -1
    /\
    a0 - a1 - a2 - a3 + b0 + b1 - b3 >= -3
    /\
    - a0 - a1 + a2 + a3 - b0 - b2 - b3 >= -4
    /\
    a0 + a1 - a2 + a3 - b0 + b2 - b3 >= -2
    /\
    a0 + a1 - a2 - a3 - b0 - b2 + p1 >= -3
    /\
    a0 + a1 - a2 - a3 + b0 + b2 + p1 >= -1
    /\
    - a0 - a1 + a2 + a3 + b0 + b2 + p1 >= -1
    /\
    a0 - a1 + a2 + a3 + b1 + b2 + p1 >= 0
    /\
    a0 + a1 + a3 + b0 - b2 - b3 + p1 >= -1
    /\
    - a1 - a3 + b0 + b1 - b2 - b3 + p1 >= -3
    /\
    - a1 - a3 - b0 - b1 + b2 - b3 + p1 >= -4
    /\
    a0 + a1 + b0 + b2 + b3 - p0 + p1 >= 0
    /\
    - a0 - a1 - a2 - a3 + b0 - b1 + b2 - b3 >= -5
    /\
    - a0 + a1 + a2 - a3 + b0 - b1 - b3 + p1 >= -3
    /\
    a0 - a1 - a2 + a3 + b0 - b1 - b3 + p1 >= -3
    /\
    a0 - a1 + a2 - a3 + b1 - b2 - b3 + p1 >= -3
    /\
    a0 - a1 + a2 - a3 - b1 + b2 - b3 + p1 >= -3
    /\
    - a0 + a1 - a2 + a3 - b1 + b2 - b3 + p1 >= -3
    /\
    - a0 + a1 - a2 + a3 + b1 - b2 - b3 >= -3
    /\
    - a0 + a1 + a2 - a3 - b0 + b1 - b3 >= -3
    /\
    a0 - a1 - a2 + a3 - b0 + b1 - b3 >= -3
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.02 seconds
Number of constraints: 51
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 4.0000 p0 + 2.0000 p1
*/
predicate s4_lin(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                 var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                 var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    - a1 + a2 - b0 + b1 + p0 >= -1
    /\
    a0 + a1 + a3 - b2 + p0 >= 0
    /\
    - a0 - a3 - b0 - b2 + p0 >= -3
    /\
    - a2 - a3 - b1 - b2 + p0 >= -3
    /\
    a2 - a3 + b1 - b2 + p0 >= -1
    /\
    - a0 + a1 + a3 + b2 + p0 >= 0
    /\
    - a1 - a2 - b0 + b2 + p0 >= -2
    /\
    - a0 + a2 - b1 - b3 + p0 >= -2
    /\
    - a1 + b0 - b2 - b3 + p0 >= -2
    /\
    a0 + a2 - a3 + b3 + p0 >= 0
    /\
    - a2 + a3 + b0 + b3 + p0 >= 0
    /\
    a1 + a2 + a3 + b2 - p1 >= 0
    /\
    a2 + a3 + b2 + b3 - p1 >= 0
    /\
    a1 - a2 - a3 + b1 + p1 >= -1
    /\
    - a2 - a3 + b0 + b1 + p1 >= -1
    /\
    - a0 + a1 - a3 + b2 + p1 >= -1
    /\
    - a1 + b0 + b1 + b2 + p1 >= 0
    /\
    a1 + b1 + b2 - p0 + p1 >= 0
    /\
    - a0 - a1 - a2 - b0 - b1 - b2 >= -5
    /\
    - a1 + a2 + a3 + b0 + b1 - b2 >= -1
    /\
    a0 - a1 + a2 + a3 + b0 + p0 >= 0
    /\
    - a0 + a1 - a3 + b0 + b1 + p0 >= -1
    /\
    - a0 + a2 - a3 + b0 + b2 + p0 >= -1
    /\
    a1 - a3 - b0 - b1 + b2 + p0 >= -2
    /\
    a0 - a2 - a3 + b1 + b2 + p0 >= -1
    /\
    - a0 - a2 + a3 + b1 - b3 + p0 >= -2
    /\
    a3 - b0 - b1 - b2 - b3 + p0 >= -3
    /\
    - a3 - b0 + b1 - b2 - b3 + p0 >= -3
    /\
    - a1 - a3 - b1 + b2 - b3 + p0 >= -3
    /\
    a0 + b0 - b1 + b2 - b3 + p0 >= -1
    /\
    a0 + a1 + b1 + b2 - b3 + p0 >= 0
    /\
    - a0 + a2 + a3 - b0 + b3 + p0 >= -1
    /\
    - a0 + a1 + a2 + b0 + b3 + p0 >= 0
    /\
    a0 + a1 - a3 + b0 + b3 + p0 >= 0
    /\
    a1 + a3 - b0 + b1 + b3 + p0 >= 0
    /\
    - a0 - a1 + b0 + b2 + b3 + p0 >= -1
    /\
    - a0 + a1 - a2 - b1 - b2 + p1 >= -3
    /\
    a0 + a1 + a2 - b1 - b2 + p1 >= -1
    /\
    a0 + a2 - b0 - b1 - b2 + p1 >= -2
    /\
    a1 + a2 + a3 + b1 - b2 + p1 >= 0
    /\
    a0 + a1 + a3 - b1 + b2 + p1 >= 0
    /\
    a0 + a3 - b0 - b1 + b2 + p1 >= -1
    /\
    a0 - a1 - a2 - b0 - b2 + b3 + p0 >= -3
    /\
    - a0 - a1 + a2 + b0 - b1 - b2 + p1 >= -3
    /\
    - a1 - a2 + a3 - b0 + b1 - b2 + p1 >= -3
    /\
    - a0 - a1 - a3 - b0 - b1 + b2 + p1 >= -4
    /\
    a0 - a1 - a3 + b0 - b1 + b2 + p1 >= -2
    /\
    - a0 - a1 + a3 + b0 - b1 + b2 + p1 >= -2
    /\
    - a1 + a2 - a3 - b0 + b1 - b2 >= -3
    /\
    a0 - a1 - a2 + b0 - b1 - b2 >= -3
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.02 seconds
Number of constraints: 50
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 4.0000 p0 + 2.0000 p1
*/
predicate s5_lin(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                 var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                 var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    a0 + a1 + a2 + a3 - b2 >= 0
    /\
    a1 + a2 - a3 + b0 + p0 >= 0
    /\
    a1 - a2 + a3 + b0 + p0 >= 0
    /\
    - a0 - a2 - b0 - b1 + p0 >= -3
    /\
    - a0 + a3 - b0 - b1 + p0 >= -2
    /\
    - a0 + a1 + a2 + b1 + p0 >= 0
    /\
    a0 - a2 - b0 + b1 + p0 >= -1
    /\
    - a2 - a3 - b0 - b2 + p0 >= -3
    /\
    - a0 - a1 + b0 - b2 + p0 >= -2
    /\
    a0 - a1 + b1 - b2 + p0 >= -1
    /\
    a2 + a3 - b1 + b2 + p0 >= 0
    /\
    a0 - a3 - b1 - b3 + p0 >= -2
    /\
    - a1 - b0 + b2 - b3 + p0 >= -2
    /\
    a0 - a2 + a3 + b3 + p0 >= 0
    /\
    a1 - b0 + b2 + b3 + p0 >= 0
    /\
    a1 - a2 - a3 + b0 + p1 >= -1
    /\
    - a0 + a1 - a2 + b1 + p1 >= -1
    /\
    - a0 - a2 + b1 + b2 + p1 >= -1
    /\
    - a1 + b0 + b1 + b2 + p1 >= 0
    /\
    a1 + b0 + b1 - p0 + p1 >= 0
    /\
    - a0 - a1 - a3 - b0 - b1 - b2 >= -5
    /\
    a0 - a1 + a2 - b0 + b1 + b2 >= -1
    /\
    - a0 - a2 + a3 + b0 + b1 + p0 >= -1
    /\
    a1 - a2 + b0 - b1 - b2 + p0 >= -2
    /\
    a0 - a2 - a3 + b0 + b2 + p0 >= -1
    /\
    - a1 - a2 + b0 - b1 - b3 + p0 >= -3
    /\
    a2 - b0 - b1 - b2 - b3 + p0 >= -3
    /\
    - a2 + a3 + b1 - b2 - b3 + p0 >= -2
    /\
    a0 + a2 + a3 + b2 - b3 + p0 >= 0
    /\
    a0 + a2 - a3 - b2 + b3 + p0 >= -1
    /\
    - a0 - a1 + a3 - b2 + b3 + p0 >= -2
    /\
    - a0 - a1 + a2 + b2 + b3 + p0 >= -1
    /\
    - a1 - a3 + b0 + b2 + b3 + p0 >= -1
    /\
    a1 - a2 - a3 + b1 + b2 - p1 >= -2
    /\
    - a0 + a2 - a3 + b1 - b3 - p1 >= -3
    /\
    a0 + a2 + b0 + b3 + p0 - p1 >= 0
    /\
    - a0 + a1 - a3 - b0 - b1 + p1 >= -3
    /\
    a0 + a1 + a3 - b0 - b1 + p1 >= -1
    /\
    a1 + a2 + a3 + b0 - b1 + p1 >= 0
    /\
    a0 + a1 + a2 - b0 + b1 + p1 >= 0
    /\
    a0 + a3 - b0 - b1 - b2 + p1 >= -2
    /\
    a2 + a3 + b0 - b1 - b2 + p1 >= -1
    /\
    - a1 - a2 - a3 + b0 - b1 - b2 + p1 >= -4
    /\
    - a0 - a1 + a2 - b0 + b1 - b2 + p1 >= -3
    /\
    a0 - a1 - a3 - b0 - b1 + b2 + p1 >= -3
    /\
    - a1 + a2 - a3 + b0 - b1 + b2 + p1 >= -2
    /\
    - a1 - a2 + a3 + b0 - b1 + b2 + p1 >= -2
    /\
    - a0 - a1 + a3 - b0 - b1 + b2 >= -3
    /\
    a0 - a1 - a2 - b0 + b1 - b2 >= -3
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.02 seconds
Number of constraints: 48
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 4.0000 p0 + 2.0000 p1
*/
predicate s6_lin(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                 var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                 var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    a0 + a1 + a3 - p0 >= 0
    /\
    b0 + b1 + b3 - p0 >= 0
    /\
    - a1 - b0 - b2 + p0 >= -2
    /\
    - a0 + a1 - a3 + p1 >= -1
    /\
    b0 - b1 - b3 + p1 >= -1
    /\
    a0 + a1 - a2 + a3 + b0 >= 0
    /\
    - a0 + a1 - a2 - a3 - b2 >= -3
    /\
    - a1 - a2 + b0 - b1 - b3 >= -3
    /\
    a1 + a2 + b0 - b1 - b3 >= -1
    /\
    - a0 - a3 - b1 + b2 - b3 >= -3
    /\
    a0 - a2 - a3 - b0 + p0 >= -2
    /\
    - a0 - a2 + a3 - b0 + p0 >= -2
    /\
    a0 - a1 + b1 + b2 + p0 >= 0
    /\
    - a1 + a3 + b2 + b3 + p0 >= 0
    /\
    a0 + a1 + a3 + b2 - p1 >= 0
    /\
    a0 + a3 + b0 + b2 - p1 >= 0
    /\
    a2 + b0 + b1 + b3 - p1 >= 0
    /\
    a0 - a1 - a3 + b0 + p1 >= -1
    /\
    - a0 - a1 + a3 + b0 + p1 >= -1
    /\
    - a0 - a3 - b1 - b3 + p1 >= -3
    /\
    a0 + a3 - b1 - b3 + p1 >= -1
    /\
    a1 - b0 + b1 - b3 + p1 >= -1
    /\
    a1 - b0 - b1 + b3 + p1 >= -1
    /\
    - a0 - a3 + b1 + b3 + p1 >= -1
    /\
    - a0 + a1 + a3 + b0 + b1 + p0 >= 0
    /\
    a1 + a2 + a3 - b1 - b2 + p0 >= -1
    /\
    a0 + a1 - a3 + b1 - b2 + p0 >= -1
    /\
    a0 + a2 + a3 + b1 - b2 + p0 >= 0
    /\
    a2 - a3 - b0 - b1 + b2 + p0 >= -2
    /\
    - a0 - a1 - a3 + b1 - b3 + p0 >= -3
    /\
    a0 - a2 + b0 + b1 - b3 + p0 >= -1
    /\
    - a0 + a2 - b0 + b2 - b3 + p0 >= -2
    /\
    a2 + b0 + b1 + b2 - b3 + p0 >= 0
    /\
    - a0 - a1 - a3 - b1 + b3 + p0 >= -3
    /\
    - a2 + a3 + b0 - b1 + b3 + p0 >= -1
    /\
    - a0 + a1 + a3 - b2 + b3 + p0 >= -1
    /\
    a0 + b0 - b1 - b2 + b3 + p0 >= -1
    /\
    a0 + a1 - a3 + b2 + b3 + p0 >= 0
    /\
    a2 + b0 - b1 + b2 + b3 + p0 >= 0
    /\
    a1 - b0 + b1 + b2 + b3 + p0 >= 0
    /\
    - a0 + b0 + b1 + b2 + b3 + p0 >= 0
    /\
    - a0 - a1 + a3 - b2 - b3 - p1 >= -4
    /\
    a0 - a2 + a3 + b1 + b3 + p1 >= 0
    /\
    a0 + a3 - b0 + b1 + b3 + p1 >= 0
    /\
    a0 + a2 - a3 - b1 - b2 - b3 + p0 >= -3
    /\
    - a0 + a2 - a3 - b0 + b1 - b3 - p1 >= -4
    /\
    - a0 + a2 - a3 - b0 - b1 + b3 - p1 >= -4
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 70
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 4.0000 p0 + 2.0000 p1
*/
predicate s7_lin(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                 var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                 var 0..1: p0, var 0..1: p1) = 
    - p0 - p1 >= -1
    /\
    b0 + b1 + b2 + b3 - p0 >= 0
    /\
    a0 - a2 + a3 + b0 + p0 >= 0
    /\
    - a0 + a2 + b0 - b1 + p0 >= -1
    /\
    a1 + a3 + b0 - b1 + p0 >= 0
    /\
    - a1 - a3 - b1 - b2 + p0 >= -3
    /\
    - a1 + a3 + b1 - b2 + p0 >= -1
    /\
    a2 + a3 + b1 - b2 + p0 >= 0
    /\
    a0 - a1 + a2 + b2 + p0 >= 0
    /\
    - a2 - b0 - b1 + b2 + p0 >= -2
    /\
    a1 - a3 + b1 - b3 + p0 >= -1
    /\
    - a2 - a3 - b2 - b3 + p0 >= -3
    /\
    a2 + b0 - b1 + b3 + p0 >= 0
    /\
    a1 + a3 - b2 + b3 + p0 >= 0
    /\
    a0 + b1 - b2 + b3 + p0 >= 0
    /\
    a1 - a3 + b2 + b3 + p0 >= 0
    /\
    - a0 + b1 + b2 + b3 + p0 >= 0
    /\
    - b0 - b1 - b2 - b3 + p1 >= -3
    /\
    - a0 + a1 - a3 - b0 - b2 - b3 >= -4
    /\
    a1 + a2 + a3 + b0 + b1 - p0 >= 0
    /\
    a1 + a2 - a3 - b0 - b2 + p0 >= -2
    /\
    a0 + a1 + a3 - b0 + b2 + p0 >= 0
    /\
    a0 + a1 - a3 + b0 + b2 + p0 >= 0
    /\
    - a0 - a1 - a3 - b0 - b3 + p0 >= -4
    /\
    - a0 - a1 + b0 - b1 - b3 + p0 >= -3
    /\
    - a0 - a2 - b0 - b2 - b3 + p0 >= -4
    /\
    a1 + a2 + a3 + b2 - b3 + p0 >= 0
    /\
    - a1 + a3 - b0 + b2 - b3 + p0 >= -2
    /\
    - a2 + a3 + b0 + b2 - b3 + p0 >= -1
    /\
    - a0 - a1 + a2 + a3 + b3 + p0 >= -1
    /\
    - a0 - a2 - a3 - b1 + b3 + p0 >= -3
    /\
    - a1 - a2 - b0 - b2 + b3 + p0 >= -3
    /\
    a0 + a2 + a3 - b0 - b3 - p1 >= -2
    /\
    a0 + b0 + b1 + b2 + p0 - p1 >= 0
    /\
    a0 - a2 - a3 - b1 - b2 + p1 >= -3
    /\
    a0 + a2 + a3 - b1 - b2 + p1 >= -1
    /\
    - a2 - a3 - b0 - b1 - b2 + p1 >= -4
    /\
    a0 - a2 - a3 + b1 + b2 + p1 >= -1
    /\
    a0 - a1 - a2 - b1 - b3 + p1 >= -3
    /\
    a0 + a1 + a2 - b1 - b3 + p1 >= -1
    /\
    - a1 - a2 - b0 - b1 - b3 + p1 >= -4
    /\
    a1 + a2 - b0 - b1 - b3 + p1 >= -2
    /\
    a0 - a1 - a3 - b2 - b3 + p1 >= -3
    /\
    a0 + a1 + a3 - b2 - b3 + p1 >= -1
    /\
    a0 - a1 - a2 + b1 + b3 + p1 >= -1
    /\
    a0 - a1 - a3 + b2 + b3 + p1 >= -1
    /\
    - a0 + a2 + a3 - b0 - b1 - b2 + b3 >= -3
    /\
    - a0 - a1 + a2 - a3 + b0 + b1 + p0 >= -2
    /\
    - a0 + a1 - a2 - a3 + b0 - b2 + p0 >= -3
    /\
    - a0 - a1 - a2 - a3 + b0 + b1 + p1 >= -3
    /\
    - a0 + a1 - a2 - a3 + b0 + b2 + p1 >= -2
    /\
    - a0 - a1 + a2 + a3 + b0 + b2 + p1 >= -1
    /\
    - a0 - a1 + a3 - b0 - b2 - b3 + p1 >= -4
    /\
    - a1 - a3 + b0 + b1 - b2 - b3 + p1 >= -3
    /\
    - a0 + a1 + a2 - a3 + b0 + b3 + p1 >= -1
    /\
    - a0 - a1 - a2 + a3 + b0 + b3 + p1 >= -2
    /\
    a0 + a2 + a3 + b1 + b2 - p0 + p1 >= 0
    /\
    a0 + a1 + a2 + b1 + b3 - p0 + p1 >= 0
    /\
    a0 + a1 + a3 + b2 + b3 - p0 + p1 >= 0
    /\
    - a0 + a1 - a2 + b0 - b1 + b2 - b3 + p1 >= -3
    /\
    - a0 - a2 + a3 - b0 + b1 + b2 - b3 + p1 >= -3
    /\
    - a0 - a2 + a3 + b0 - b1 - b2 + b3 + p1 >= -3
    /\
    - a0 + a1 - a2 - b0 + b1 - b2 + b3 + p1 >= -3
    /\
    - a0 - a1 + a2 - b0 + b1 - b2 + b3 + p1 >= -3
    /\
    - a0 - a1 - a3 - b0 - b1 + b2 + b3 + p1 >= -4
    /\
    - a0 + a1 + a3 - b0 - b1 + b2 + b3 + p1 >= -2
    /\
    - a0 + a1 - a2 + a3 + b0 + b1 - b2 - b3 + p1 >= -3
    /\
    - a0 + a2 - a3 - b0 + b1 + b2 - b3 >= -3
    /\
    - a0 + a2 - a3 + b0 - b1 - b2 + b3 >= -3
    /\
    - a0 - a1 + a2 + b0 - b1 + b2 - b3 >= -3
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 36
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s0_lin_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                      var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) = 
    a0 + a1 + a2 + a3 - b0 >= 0
    /\
    a0 + a2 + a3 + b0 - b1 >= 0
    /\
    a0 + a2 + a3 + b1 - b2 >= 0
    /\
    a0 + a2 + a3 + b1 - b3 >= 0
    /\
    - a0 + a1 + b1 + b2 + b3 >= 0
    /\
    a1 - a2 + b1 + b2 + b3 >= 0
    /\
    - a1 + a2 + b1 + b2 + b3 >= 0
    /\
    a1 - a3 + b1 + b2 + b3 >= 0
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 - a1 + a2 - a3 - b0 - b1 >= -4
    /\
    a0 - a2 - a3 - b0 + b1 - b2 >= -3
    /\
    a0 - a2 - a3 + b0 + b1 + b2 >= -1
    /\
    - a0 - a2 + a3 - b0 + b1 - b3 >= -3
    /\
    a0 + a3 + b0 + b1 - b2 - b3 >= -1
    /\
    - a0 - a2 + a3 + b0 + b1 + b3 >= -1
    /\
    a0 - a2 + a3 - b0 - b1 - b2 - b3 >= -4
    /\
    a0 - a1 - a2 - a3 + b1 + b2 - b3 >= -3
    /\
    - a0 - a1 - a2 + a3 + b1 + b2 - b3 >= -3
    /\
    a0 - a1 - a2 - a3 + b1 - b2 + b3 >= -3
    /\
    - a0 - a1 - a2 + a3 + b1 - b2 + b3 >= -3
    /\
    - a0 - a2 - a3 - b0 - b1 + b2 + b3 >= -4
    /\
    a2 - a3 - b0 - b1 + b2 - b3 >= -3
    /\
    - a0 + a2 - a3 + b0 + b2 + b3 >= -1
    /\
    - a0 - a2 - a3 - b1 - b2 - b3 >= -5
    /\
    - a1 + a2 - a3 - b1 - b2 + b3 >= -3
    /\
    - a0 + a2 - b0 - b1 - b2 + b3 >= -3
    /\
    a0 + a2 + b0 - b1 - b2 + b3 >= -1
    /\
    a2 + a3 + b0 - b1 + b2 - b3 >= -1
    /\
    a0 - a2 + a3 - b1 + b2 + b3 >= -1
    /\
    - a1 - b0 + b1 - b2 - b3 >= -3
    /\
    a1 + a2 + a3 - b1 - b2 + b3 >= -1
    /\
    - a0 - a1 + a2 - b1 + b2 - b3 >= -3
    /\
    a0 + a1 + a2 - b1 + b2 - b3 >= -1
    /\
    a1 + b0 + b1 - b2 - b3 >= -1
    /\
    - a0 + a1 + a2 - a3 + b1 >= -1
    /\
    - a0 - a3 + b0 - b2 - b3 >= -3
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 36
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s1_lin_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                      var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) = 
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    a1 + a2 + a3 + b1 - b3 >= 0
    /\
    a1 + a2 + a3 - b0 + b3 >= 0
    /\
    a1 + a2 + a3 - b2 + b3 >= 0
    /\
    - a0 + a1 + b0 + b2 + b3 >= 0
    /\
    a1 - a2 + b0 + b2 + b3 >= 0
    /\
    a0 - a3 + b0 + b2 + b3 >= 0
    /\
    - a1 + b0 + b1 + b2 + b3 >= 0
    /\
    a1 - a2 - a3 + b0 - b1 + b2 >= -2
    /\
    a0 + a2 + a3 + b0 - b1 + b2 >= 0
    /\
    - a0 + a1 - a3 + b0 - b2 - b3 >= -3
    /\
    a0 + a1 + a3 + b0 - b2 - b3 >= -1
    /\
    - a0 + a1 - a2 - b0 + b2 - b3 >= -3
    /\
    a0 + a1 + a2 - b0 + b2 - b3 >= -1
    /\
    - a0 - a1 - b0 - b1 - b2 + b3 >= -4
    /\
    - a0 - a1 + a2 + a3 - b0 - b2 - b3 >= -4
    /\
    a1 + a2 - a3 + b0 - b1 - b2 - b3 >= -3
    /\
    a1 - a2 + a3 + b0 - b1 - b2 - b3 >= -3
    /\
    - a0 - a1 - a2 - a3 + b0 + b2 - b3 >= -4
    /\
    a1 + a2 - a3 - b0 - b1 + b2 - b3 >= -3
    /\
    a0 + a1 + a3 - b0 - b1 + b2 - b3 >= -2
    /\
    - a1 + a2 - a3 + b1 + b2 + b3 >= -1
    /\
    a0 - a1 + a2 - a3 + b0 + b3 >= -1
    /\
    a0 + a2 + a3 - b0 - b2 + b3 >= -1
    /\
    - a0 + a1 - a2 - a3 - b1 >= -3
    /\
    - a1 - a2 - a3 - b0 - b2 - b3 >= -5
    /\
    - a1 - a2 + a3 + b0 + b1 + b3 >= -1
    /\
    - a1 + a2 + a3 + b0 + b2 - b3 >= -1
    /\
    - a1 - a2 + a3 - b0 - b1 + b3 >= -3
    /\
    - a1 + a2 - a3 - b1 - b2 + b3 >= -3
    /\
    - a0 - a1 - a2 + a3 - b2 + b3 >= -3
    /\
    a0 - a1 - a2 + a3 + b2 + b3 >= -1
    /\
    - a0 - a1 + a2 - a3 - b0 + b3 >= -3
    /\
    a0 + a1 - a2 - a3 + b1 >= -1
    /\
    a1 - b0 + b1 - b2 + b3 >= -1
    /\
    a0 - a2 - a3 - b0 - b2 >= -3
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 32
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s2_lin_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                      var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) = 
    a0 + a1 + a2 + a3 - b2 >= 0
    /\
    a0 + a3 - b0 + b1 + b2 >= 0
    /\
    - a1 + a2 + b0 + b1 + b2 >= 0
    /\
    - a2 - b0 - b1 - b2 - b3 >= -4
    /\
    a0 + b0 + b1 + b2 - b3 >= 0
    /\
    a0 + a1 + a2 - a3 + b3 >= 0
    /\
    - a0 + b0 + b1 + b2 + b3 >= 0
    /\
    a1 + a2 + a3 + b0 - b1 + b2 >= 0
    /\
    a0 + a1 + a2 + b0 + b1 - b3 >= 0
    /\
    a0 - a1 + b0 - b1 - b2 - b3 >= -3
    /\
    a0 + a1 + b0 - b1 - b2 + b3 >= -1
    /\
    a0 - a1 - a2 + b0 + b2 + b3 >= -1
    /\
    a0 + a1 - a2 + b1 + b2 + b3 >= 0
    /\
    a0 + a1 - a2 - a3 + b0 - b1 - b2 >= -3
    /\
    a0 - a1 + a2 - a3 + b0 - b1 - b2 >= -3
    /\
    - a0 + a1 + a2 - a3 - b0 + b1 - b2 >= -3
    /\
    - a0 - a1 - a2 + a3 + b0 - b1 + b2 >= -3
    /\
    a1 + a2 + a3 - b0 + b1 - b2 - b3 >= -2
    /\
    - a1 - a2 - a3 + b0 - b1 + b2 - b3 >= -4
    /\
    - a0 - a1 - a2 - b0 - b2 - b3 >= -5
    /\
    - a0 - a1 + a2 - a3 - b0 - b1 >= -4
    /\
    - a0 + a1 - a2 - b0 - b1 + b2 >= -3
    /\
    - a0 - a1 + a2 + a3 + b0 + b1 >= -1
    /\
    - a0 + a1 - a2 + b0 + b1 - b2 >= -2
    /\
    - a0 + a1 + a2 + b0 + b2 - b3 >= -1
    /\
    a2 - b0 - b1 - b2 + b3 >= -2
    /\
    a0 + a1 + a2 - b0 - b1 >= -1
    /\
    - a0 - a1 - b0 - b1 - b2 >= -4
    /\
    a0 - a1 - a2 - a3 - b3 >= -3
    /\
    - a1 - a2 - a3 - b0 - b2 >= -4
    /\
    a0 - a1 - a2 + a3 + b3 >= -1
    /\
    a0 - a1 + b1 + b2 - b3 >= -1
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 41
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s3_lin_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                      var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) = 
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    - a3 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 + a2 + a3 + b0 + b1 + b2 >= 0
    /\
    - a1 + a2 + a3 + b0 + b1 + b2 >= 0
    /\
    - a0 - a1 - a2 + a3 - b1 + b3 >= -3
    /\
    a2 + a3 - b0 - b1 - b2 + b3 >= -2
    /\
    a1 + a2 + b0 + b1 - b2 + b3 >= 0
    /\
    a0 + a1 + b0 - b1 + b2 + b3 >= 0
    /\
    a0 + a2 - b0 + b1 + b2 + b3 >= 0
    /\
    a0 - a1 - a2 - a3 + b0 + b1 - b3 >= -3
    /\
    - a1 + a2 - a3 + b0 + b1 - b2 - b3 >= -3
    /\
    a1 - a2 + a3 + b0 + b1 - b2 - b3 >= -2
    /\
    - a0 + a1 - a2 + b0 - b1 - b2 + b3 >= -3
    /\
    a0 - a1 + a2 - a3 + b0 - b1 + b2 - b3 >= -3
    /\
    a0 - a1 - a2 + a3 + b0 - b1 + b2 - b3 >= -3
    /\
    - a0 + a1 - a2 + a3 + b0 - b1 + b2 - b3 >= -3
    /\
    - a0 - a1 + a2 + a3 - b0 - b2 - b3 >= -4
    /\
    - a0 + a1 + a2 - a3 - b0 + b1 - b3 >= -3
    /\
    - a0 + a1 - a2 + a3 + b1 - b2 - b3 >= -3
    /\
    - a0 - a3 + b0 - b1 - b2 - b3 >= -4
    /\
    a0 - a1 - a2 + a3 - b0 + b1 - b3 >= -3
    /\
    a0 - a1 - a2 + a3 + b0 + b1 + b3 >= -1
    /\
    - a0 + a1 - a2 + a3 + b1 + b2 + b3 >= -1
    /\
    - a0 + a1 - a2 - a3 + b0 + b1 + b2 >= -2
    /\
    a0 + a1 - a2 - a3 + b0 - b1 + b2 >= -2
    /\
    - a1 - a2 - a3 - b0 - b1 - b3 >= -5
    /\
    - a2 - a3 - b0 - b1 - b2 >= -4
    /\
    - a0 - a1 - a3 - b0 + b2 - b3 >= -4
    /\
    a0 - a1 + a2 - b1 - b2 + b3 >= -2
    /\
    - a0 - a1 - a2 - a3 - b1 - b3 >= -5
    /\
    - a0 - a1 - a2 - b0 - b2 + b3 >= -4
    /\
    a0 + a2 + a3 + b0 - b1 - b2 >= -1
    /\
    a1 + a2 + a3 - b0 - b1 + b2 >= -1
    /\
    a0 + a1 - b0 + b1 - b2 + b3 >= -1
    /\
    - a0 - a1 - a2 - a3 + b1 + b3 >= -3
    /\
    a1 + a2 - b0 - b1 + b2 + b3 >= -1
    /\
    a0 + a1 + a3 - b0 + b2 - b3 >= -1
    /\
    - a0 - a1 + a3 + b0 + b2 + b3 >= -1
    /\
    - a1 - b0 - b1 - b2 + b3 >= -3
    /\
    - a0 - a1 - a2 - a3 - b0 >= -4
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 36
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s4_lin_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                      var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) = 
    a0 + a1 + a2 + a3 - b1 >= 0
    /\
    a0 + a1 - a3 + b1 + b2 >= 0
    /\
    a1 + a3 - b0 + b1 + b2 >= 0
    /\
    - a0 + a1 + b0 + b1 + b2 >= 0
    /\
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    - a1 + b0 + b1 + b2 + b3 >= 0
    /\
    - a1 + a2 + a3 + b0 + b1 - b2 >= -1
    /\
    a0 + a1 - a2 + a3 + b0 + b2 >= 0
    /\
    - a0 + a1 + a2 - a3 - b1 + b2 >= -2
    /\
    a0 - a1 - a2 - a3 + b0 - b3 >= -3
    /\
    a0 + a2 + a3 + b0 + b1 - b3 >= 0
    /\
    a1 + a2 + a3 + b1 - b2 + b3 >= 0
    /\
    a0 - a1 + a2 - a3 + b0 - b1 + b2 >= -2
    /\
    - a0 - a1 + a2 + a3 + b0 - b1 - b3 >= -3
    /\
    a1 - a2 - a3 - b0 + b1 - b2 - b3 >= -4
    /\
    - a1 - a2 - a3 + b0 + b1 - b2 - b3 >= -4
    /\
    a0 + a2 - a3 - b0 - b1 - b2 + b3 >= -3
    /\
    a0 + a2 + a3 - b0 - b1 + b2 + b3 >= -1
    /\
    a0 - a1 - a2 + a3 - b0 + b1 - b2 + b3 >= -3
    /\
    - a0 - a1 - a2 - b0 - b1 - b2 >= -5
    /\
    - a0 - a2 - b0 - b1 - b2 - b3 >= -5
    /\
    - a0 + a1 - a2 + b0 - b1 - b2 + b3 >= -3
    /\
    - a0 - a1 - a2 + a3 - b0 - b2 - b3 >= -5
    /\
    a0 + a2 + a3 - b0 - b1 - b2 - b3 >= -3
    /\
    - a0 - a1 + a2 + b0 - b1 - b2 - b3 >= -4
    /\
    a0 - a1 - a2 + a3 - b0 - b1 + b2 >= -3
    /\
    a1 - a2 - a3 + b0 + b1 + b3 >= -1
    /\
    - a0 - a2 + a3 + b0 + b1 + b2 >= -1
    /\
    - a1 + a2 - a3 - b0 + b1 - b2 >= -3
    /\
    - a0 - a2 - a3 - b0 - b1 >= -4
    /\
    - a0 - a3 - b0 - b1 + b2 - b3 >= -4
    /\
    a0 + a1 + a2 - b1 - b2 + b3 >= -1
    /\
    a0 - a1 - a2 + b0 - b1 - b2 >= -3
    /\
    - a0 - a1 + a3 + b0 + b2 + b3 >= -1
    /\
    - a0 + a2 - a3 + b0 + b1 + b2 >= -1
    /\
    - a0 + a1 - a2 - a3 - b2 >= -3
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 36
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s5_lin_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                      var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) = 
    a0 + a1 + a2 + a3 - b1 >= 0
    /\
    a1 - a2 + a3 + b0 + b1 >= 0
    /\
    a1 + a2 + b0 + b1 - b2 >= 0
    /\
    a1 - a3 + b0 + b1 + b2 >= 0
    /\
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    - a1 + b0 + b1 + b2 + b3 >= 0
    /\
    a0 + a1 - a2 - a3 + b0 - b1 >= -2
    /\
    - a0 + a1 + a2 + a3 + b0 + b2 >= 0
    /\
    a0 - a1 + a2 - b0 + b1 + b2 >= -1
    /\
    - a0 - a1 - a2 + a3 + b2 - b3 >= -3
    /\
    a0 + a2 + a3 + b1 + b2 - b3 >= 0
    /\
    a0 + a1 + a2 - b0 + b1 + b3 >= 0
    /\
    a0 - a1 - a2 + a3 + b0 - b1 + b2 >= -2
    /\
    - a0 + a1 - a2 - b0 + b1 - b2 - b3 >= -4
    /\
    a0 - a1 + a2 - a3 - b1 + b2 - b3 >= -3
    /\
    - a0 - a1 - a2 - b0 + b1 + b2 - b3 >= -4
    /\
    a0 - a2 + a3 - b0 - b1 - b2 + b3 >= -3
    /\
    - a1 + a2 + a3 + b0 - b1 - b2 + b3 >= -2
    /\
    - a0 - a1 + a2 + a3 - b0 + b1 - b2 + b3 >= -3
    /\
    - a0 - a1 + a2 - a3 - b0 - b2 - b3 >= -5
    /\
    a0 + a2 + a3 - b0 - b1 - b2 - b3 >= -3
    /\
    - a0 - a1 + a2 + a3 + b0 - b1 - b2 >= -3
    /\
    - a0 + a1 - a3 - b0 - b1 + b2 + b3 >= -3
    /\
    a0 - a1 - a3 - b0 - b1 + b2 - b3 >= -4
    /\
    a0 - a2 - a3 + b0 + b1 + b2 >= -1
    /\
    - a2 - a3 + b0 - b1 - b2 - b3 >= -4
    /\
    - a0 - a3 - b0 - b1 - b2 - b3 >= -5
    /\
    - a0 - a1 - a3 - b0 - b1 - b2 >= -5
    /\
    a0 - a1 - a2 - b0 + b1 - b2 >= -3
    /\
    a0 + a1 + a3 - b0 - b1 + b3 >= -1
    /\
    - a1 + a2 - a3 + b0 + b2 + b3 >= -1
    /\
    - a0 - a1 + a3 - b0 - b1 + b2 >= -3
    /\
    - a0 + a1 - a2 + b1 + b2 + b3 >= -1
    /\
    - a0 + a2 - a3 + b0 + b1 + b2 >= -1
    /\
    - a0 - a2 - a3 - b1 - b2 >= -4
    /\
    - a0 + a1 - a2 - a3 - b0 >= -3
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 36
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s6_lin_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                      var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) = 
    a0 + a1 - a2 + a3 + b0 >= 0
    /\
    a0 + a1 + a2 + a3 - b1 >= 0
    /\
    a0 + a1 + a3 - b0 + b2 >= 0
    /\
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    - a0 + a1 + b0 + b1 + b3 >= 0
    /\
    a1 - a3 + b0 + b1 + b3 >= 0
    /\
    a2 - a3 + b0 + b1 + b3 >= 0
    /\
    - a1 + b0 + b1 + b2 + b3 >= 0
    /\
    - a1 - a2 + b0 - b1 - b2 - b3 >= -4
    /\
    a1 - a3 - b0 + b1 - b2 - b3 >= -3
    /\
    a1 + a3 - b0 + b1 + b2 - b3 >= -1
    /\
    - a0 + a1 - b0 - b1 - b2 + b3 >= -3
    /\
    a0 + a2 + a3 + b1 - b2 + b3 >= 0
    /\
    - a0 + a1 - a3 + b1 + b2 + b3 >= -1
    /\
    - a0 + a1 - a2 + a3 - b0 + b1 - b3 >= -3
    /\
    a0 - a1 + a3 - b0 - b1 - b2 - b3 >= -4
    /\
    a0 + a1 - a2 - b0 + b1 + b2 - b3 >= -2
    /\
    a1 - a2 - a3 - b0 - b1 - b2 + b3 >= -4
    /\
    - a0 - a1 - a3 - b0 + b1 - b2 + b3 >= -4
    /\
    a0 + a1 - a3 - b0 - b1 + b2 + b3 >= -2
    /\
    a1 - a2 + a3 - b0 - b1 + b2 + b3 >= -2
    /\
    a0 + a3 + b0 - b1 + b2 - b3 >= -1
    /\
    a0 - a1 - a3 + b0 + b1 + b2 >= -1
    /\
    - a0 - a1 - a3 - b0 - b1 - b3 >= -5
    /\
    a0 - a1 - a3 + b0 - b1 - b2 >= -3
    /\
    - a0 - a1 + a3 + b0 - b2 - b3 >= -3
    /\
    - a0 - a1 + a3 + b0 + b2 + b3 >= -1
    /\
    - a0 - a1 - a2 + a3 + b0 - b1 >= -3
    /\
    a0 - a1 - a2 - a3 + b0 - b3 >= -3
    /\
    a0 - a1 + a3 - b0 + b1 + b3 >= -1
    /\
    a0 - a1 + a2 - a3 + b0 + b3 >= -1
    /\
    - a0 - a1 + a2 + a3 + b0 + b1 >= -1
    /\
    - a0 - a3 - b1 + b2 - b3 >= -3
    /\
    - a0 + a1 - a2 - a3 - b2 >= -3
    /\
    - a0 + a1 + a2 - a3 + b2 >= -1
    /\
    a1 + a2 + b0 - b1 - b3 >= -1
;
% #############################################################################################################################################
/*
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 42
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
*/
predicate s7_lin_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                      var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0) = 
    a0 + a1 + a2 + a3 - b2 >= 0
    /\
    a0 + a1 + a2 + a3 - b3 >= 0
    /\
    - a2 + b0 + b1 + b2 + b3 >= 0
    /\
    - a3 + b0 + b1 + b2 + b3 >= 0
    /\
    - a0 + a1 + a2 + a3 + b0 + b1 >= 0
    /\
    a0 - a1 + a2 + a3 + b1 + b2 >= 0
    /\
    a0 - a2 - a3 + b0 + b1 + b2 >= -1
    /\
    a0 - a2 + b0 - b1 - b2 - b3 >= -3
    /\
    a0 + a1 + a2 + b1 - b2 + b3 >= 0
    /\
    a0 + a1 + a3 - b1 + b2 + b3 >= 0
    /\
    a0 + a1 - b0 + b1 + b2 + b3 >= 0
    /\
    - a0 - a1 + a2 + a3 - b0 + b1 - b2 >= -3
    /\
    a0 - a1 - a2 + a3 - b1 + b2 - b3 >= -3
    /\
    - a0 + a2 + a3 - b0 - b1 - b2 + b3 >= -3
    /\
    - a0 - a1 + a2 + a3 + b0 + b2 + b3 >= -1
    /\
    - a0 + a1 - a2 + a3 - b0 + b1 - b2 + b3 >= -3
    /\
    - a0 + a1 - a2 + a3 + b0 - b1 + b2 - b3 >= -3
    /\
    - a0 + a1 - a2 + a3 + b0 - b1 - b2 + b3 >= -3
    /\
    a0 - a1 - a2 + b1 - b2 + b3 >= -2
    /\
    - a0 + a2 - a3 - b0 + b1 + b2 - b3 >= -3
    /\
    - a0 - a1 - a3 + b0 + b1 - b2 - b3 >= -4
    /\
    a1 - a2 + a3 - b0 - b1 + b2 + b3 >= -2
    /\
    a0 - a1 - a3 - b1 - b2 - b3 >= -4
    /\
    - a0 + a2 - a3 + b0 - b1 - b2 + b3 >= -3
    /\
    - a0 - a1 + a2 + b0 - b1 + b2 - b3 >= -3
    /\
    a0 + a2 - b0 - b1 - b2 - b3 >= -3
    /\
    a0 + a2 + a3 + b0 - b1 - b2 + b3 >= -1
    /\
    - a0 - a2 - b0 - b1 - b2 - b3 >= -5
    /\
    a0 + a1 + a2 + b0 - b1 + b2 - b3 >= -1
    /\
    a1 + a2 + a3 - b0 - b1 + b2 - b3 >= -2
    /\
    a0 + a1 + a3 + b0 - b2 - b3 >= -1
    /\
    a0 - a1 + a2 - a3 + b2 + b3 >= -1
    /\
    a0 + a1 - a2 - a3 + b1 + b2 >= -1
    /\
    - a1 - a2 + a3 + b0 + b1 + b3 >= -1
    /\
    a1 - a2 - a3 + b0 + b1 + b2 >= -1
    /\
    - a0 - a2 - a3 - b0 - b1 - b2 >= -5
    /\
    - a0 - a1 - a2 + a3 - b0 - b3 >= -4
    /\
    - a0 + a1 - a3 + b0 + b2 + b3 >= -1
    /\
    - a3 - b0 - b1 - b2 - b3 >= -4
    /\
    - a0 + a1 - a3 - b0 - b2 - b3 >= -4
    /\
    a0 - a1 - a2 - a3 - b2 >= -3
    /\
    - a1 - a2 - a3 - b0 - b1 >= -4
;

% #############################################################################################################################################
% We use S-box Analyzer to encode deterministic linear behaviour:
% x0: lsb (Note: we model the inverse of S-box for the linear trail through the middle part)
predicate s0_lin_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                   var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                   var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 1 /\ x0 == 1) then (y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == -1 /\ y2 == 1 /\ y1 == 1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 1 /\ x0 == 1) then (y3 == -1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;

% #############################################################################################################################################
predicate s1_lin_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                   var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                   var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 1 /\ x1 == 0 /\ x0 == 0) then (y3 == 1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 0 /\ x1 == 1 /\ x0 == 0) then (y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 1 /\ x1 == 1 /\ x0 == 0) then (y3 == -1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s2_lin_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                   var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                   var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y3 == 1 /\ y2 == -1 /\ y1 == 1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 1 /\ x1 == 1 /\ x0 == 0) then (y3 == -1 /\ y2 == 0 /\ y1 == 1 /\ y0 == -1)
    elseif (x3 == 1 /\ x2 == 1 /\ x1 == 1 /\ x0 == 1) then (y3 == -1 /\ y2 == -1 /\ y1 == 0 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s3_lin_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                   var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                   var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 1 /\ x2 == 1 /\ x1 == 1 /\ x0 == 0) then (y3 == -1 /\ y2 == 0 /\ y1 == -1 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s4_lin_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                   var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                   var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s5_lin_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                   var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                   var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 1) then (y3 == -1 /\ y2 == 1 /\ y1 == -1 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s6_lin_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                   var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                   var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 0 /\ x2 == 0 /\ x1 == 1 /\ x0 == 0) then (y3 == -1 /\ y2 == 1 /\ y1 == 1 /\ y0 == -1)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate s7_lin_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                   var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                   var 0..1: d) = 
    if (x3 == 0 /\ x2 == 0 /\ x1 == 0 /\ x0 == 0) then (y3 == 0 /\ y2 == 0 /\ y1 == 0 /\ y0 == 0)
    elseif (x3 == 1 /\ x2 == 1 /\ x1 == 1 /\ x0 == 1) then (y3 == -1 /\ y2 == -1 /\ y1 == -1 /\ y0 == 0)
    else (y3 = -1 /\ y2 = -1 /\ y1 = -1 /\ y0 = -1)
    endif
    /\
    if (x0 == 0 /\ x1 == 0 /\ x2 == 0 /\ x3 == 0) then d = 0 else d = 1 endif
;

% #############################################################################################################################################
% #############################################################################################################################################
% #############################################################################################################################################
predicate sbox_diff(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                    var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                    var 0..1: p0, var 0..1: p1, var int: which) = 
    if (which == 0) then s0_diff(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 1) then s1_diff(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 2) then s2_diff(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 3) then s3_diff(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 4) then s4_diff(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 5) then s5_diff(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 6) then s6_diff(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 7) then s7_diff(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    else false
    endif
;
% #############################################################################################################################################
predicate sbox_diff_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                         var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                         var 0..1: d, var int: which) =
    if (which == 0) then s0_diff_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 1) then s1_diff_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 2) then s2_diff_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 3) then s3_diff_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 4) then s4_diff_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 5) then s5_diff_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 6) then s6_diff_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 7) then s7_diff_star(a3, a2, a1, a0, b3, b2, b1, b0)
    else false
    endif
    /\
    if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate sbox_lin(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                   var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                   var 0..1: p0, var 0..1: p1, var int: which) = 
    if (which == 0) then s0_lin(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 1) then s1_lin(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 2) then s2_lin(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 3) then s3_lin(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 4) then s4_lin(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 5) then s5_lin(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 6) then s6_lin(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    elseif (which == 7) then s7_lin(a3, a2, a1, a0, b3, b2, b1, b0, p0, p1)
    else false
    endif
;
% #############################################################################################################################################
predicate sbox_lin_star(var 0..1: a3, var 0..1: a2, var 0..1: a1, var 0..1: a0,
                        var 0..1: b3, var 0..1: b2, var 0..1: b1, var 0..1: b0,
                        var 0..1: d, var int: which) = 
    if (which == 0) then s0_lin_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 1) then s1_lin_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 2) then s2_lin_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 3) then s3_lin_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 4) then s4_lin_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 5) then s5_lin_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 6) then s6_lin_star(a3, a2, a1, a0, b3, b2, b1, b0)
    elseif (which == 7) then s7_lin_star(a3, a2, a1, a0, b3, b2, b1, b0)
    else false
    endif
    /\
    if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0) then d = 0 else d = 1 endif
;
% #############################################################################################################################################
predicate sbox_diff_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                      var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                      var 0..1: d, var int: which) = 
    if (which == 0) then s0_diff_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 1) then s1_diff_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 2) then s2_diff_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 3) then s3_diff_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 4) then s4_diff_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 5) then s5_diff_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 6) then s6_diff_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 7) then s7_diff_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    else false
    endif
;
% #############################################################################################################################################
predicate sbox_lin_m(var -1..1: x0, var -1..1: x1, var -1..1: x2, var -1..1: x3, 
                     var -1..1: y0, var -1..1: y1, var -1..1: y2, var -1..1: y3,
                     var 0..1: d, var int: which) = 
    if (which == 0) then s0_lin_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 1) then s1_lin_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 2) then s2_lin_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 3) then s3_lin_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 4) then s4_lin_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 5) then s5_lin_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 6) then s6_lin_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    elseif (which == 7) then s7_lin_m(x0, x1, x2, x3, y0, y1, y2, y3, d)
    else false
    endif
;
