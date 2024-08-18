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

from sage.all import *
from sboxanalyzer import *
import pickle

from sage.crypto.sboxes import AES as sb
sa = SboxAnalyzer(sb)

try:
    with open('ddt.pkl', 'rb') as file:
        ddt = pickle.load(file)
    print("ddt.pkl exists")
except:
    print("ddt.pkl does not exist")
    print("Generating ddt.pkl")
    ddt = sa.difference_distribution_table()
    ddt = [[ddt[i][j] for j in range(256)] for i in range(256)]
    with open('ddt.pkl', 'wb') as file:
        pickle.dump(ddt, file)

try:
    with open('dlct.pkl', 'rb') as file:
        dlct = pickle.load(file)
    print("dlct.pkl exists")
except:
    print("dlct.pkl does not exist")
    print("Generating dlct.pkl")
    dlct = sa.differential_linear_connectivity_table()
    with open('dlct.pkl', 'wb') as file:
        pickle.dump(dlct, file)