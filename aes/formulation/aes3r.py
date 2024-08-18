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

from math import log
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm



import os
if os.path.exists('dlct.pkl'):
    with open('dlct.pkl', 'rb') as file:
        dlct = pickle.load(file)
    print("dlct.pkl exists")
else:
    print("dlct.pkl does not exist")

if os.path.exists('ddt.pkl'):
    with open('ddt.pkl', 'rb') as file:
        ddt = pickle.load(file)
    print("ddt.pkl exists")
else:
    print("ddt.pkl does not exist")

def compute_table(di_range=2, lo_range=2**8):
    m = 8
    n = 8
    big_dlct = [[0 for _ in range(2**n)] for _ in range(2**n)]
    for di in range(di_range):
        for lo in range(lo_range):
            s = 0
            for a in range(2**n):
                for b in range(2**n):
                    s += ddt[di][a] * ddt[a][b] * dlct[b][lo]            
            big_dlct[di][lo] = s            
            str_output = f"({di}, {lo}): "
            str_output += "{}".format(s)
            print(str_output)
    with open(f"aes3r.pkl", 'wb') as file:
        pickle.dump(big_dlct, file)
    return big_dlct

if __name__ == '__main__':
    di_range = 2**8
    lo_range = 2**8
    if os.path.exists('aes3r.pkl'):
        with open('aes3r.pkl', 'rb') as file:
            aes3r = pickle.load(file)
        print("aes3r.pkl exists")
    else:
        print("aes3r.pkl does not exist")
        print("Generating aes3r.pkl")
        aes3r = compute_table(di_range=di_range, lo_range=lo_range)
    
    aes3r = np.array(aes3r)

    positive_sign_flag = False
    for i in range(0, di_range):
        for j in range(0, lo_range):            
            if i != 0 and j != 0 and aes3r[i][j] > 0:
                positive_sign_flag = True
    if positive_sign_flag:
        print("Positive sign exists")
    min_values = None
    if min_values is not None:
        min_indices = np.where(aes3r == min_values)
    else:        
        min_values = np.min([aes3r[i][j] for i in range(1, di_range) for j in range(1, lo_range)])
        min_indices = np.where(aes3r == min_values)
    sn = "-" if min_values < 0 else "+"
    str_output = round(log(abs(2**(-3*8) * min_values), 2), 2)
    str_output = "{}2^({:0.2f})".format(sn, str_output)
    print("Maximum values:", str_output)
    print("Indices of maximum values:", list(zip(min_indices[0], min_indices[1])))

    aes3r_image = [[abs(2**(-3*8) * aes3r[i][j]) for j in range(1, lo_range)] for i in range(1, di_range)]
    plt.imshow(aes3r_image, cmap='jet', norm=LogNorm())
    # plt.imshow(aes3r_image, cmap='viridis', norm=LogNorm())
    plt.tick_params(axis='x', top=True, labeltop=True, bottom=False, labelbottom=False)
    plt.colorbar(label='$\mathbb{C}$ (log scale)')
    plt.title('3-Round AES')
    plt.xlabel('$\lambda{o}$')
    plt.ylabel('$\Delta_{i}$')
    plt.savefig('aes3r.svg', format='svg', dpi=800)
    plt.show()
