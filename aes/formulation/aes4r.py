"""
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
"""

import pickle
import itertools
from diff import Diff
import math

if __name__ == '__main__':
    with open('aes3r.pkl', 'rb') as file:
        big_dlct = pickle.load(file)
    params = {"nrounds" : 1,
        "variant": 1,
        "is_related_key": 0,
        "mode" : 0,
        "startweight" : 0,
        "endweight" : 128,
        "timelimit" : 10000,
        "numberoftrails" : 1,
        "fixedVariables" : {}}
    params["fixedVariables"][f"x_0"] = "00002c00000000af3f000000006b0000"
    for row, column in itertools.product(range(3), range(4)):
        params["fixedVariables"][f"x_1_{row}_{column}"] = "00"
    params["fixedVariables"][f"x_1_3_0"] = "00"
    params["fixedVariables"][f"x_1_3_1"] = "00"    
    params["fixedVariables"][f"x_1_3_3"] = "00"
    
    total_correlation = 0
    for di in range(256):
        params["fixedVariables"][f"x_1_3_2"] = hex(di)[2:].zfill(2)
        diff = Diff(params)
        diff.make_model()
        diff_trail = diff.solve()        
        if diff_trail != None:
            print(diff_trail["total_weight"])
            total_correlation += big_dlct[di][0x9] * (2**(-1*float(diff_trail["total_weight"])))
    total_correlation = math.log2(abs(total_correlation))
    print("-2^({:0.02})".format(total_correlation))


