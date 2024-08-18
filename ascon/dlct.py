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
import itertools

def binlist_to_int(L, n):
    """
    Convert a list of n bits to an integer.
    """
    return sum([L[i]*2**(n-1-i) for i in range(n)])

def compute_correlation_1(sa, dx, ly):
    """
    Compute the correlation for the middle part of DL distinguishers of Ascon
    when it involves 1 S-box with truncated input difference dx and fixed
    output linear mask ly.
    """

    dlct = sa.differential_linear_connectivity_table()
    DX = sa.truncated_to_binvectors(dx)
    LY = sa.truncated_to_binvectors(ly)
    summation = 0
    for d in DX:
        for l in LY:
            summation += dlct[binlist_to_int(d, sa.m)][binlist_to_int(l, sa.n)]
    return summation

def compute_correlation_3(sa, dx0, dx1, dx2, ly0, ly1, ly2):
    """
    Compute the correlation for the middle part of DL distinguishers of Ascon
    when it involves 3 S-boxes with truncated input differences (dx0, dx1, dx2)
    and fixed output linear masks (ly0, ly1, ly2).
    """

    dlct = sa.differential_linear_connectivity_table()
    DX0 = sa.truncated_to_binvectors(dx0)
    DX1 = sa.truncated_to_binvectors(dx1)
    DX2 = sa.truncated_to_binvectors(dx2)
    l0 = binlist_to_int(ly0, sa.n)
    l1 = binlist_to_int(ly1, sa.n)
    l2 = binlist_to_int(ly2, sa.n)
    summation = 0
    for DX in itertools.product(*[DX0, DX1, DX2]):
        d0 = binlist_to_int(DX[0], sa.m)
        d1 = binlist_to_int(DX[1], sa.m)
        d2 = binlist_to_int(DX[2], sa.m)
        temp = dlct[d0][l0]*dlct[d1][l1]*dlct[d2][l2]
        summation += temp
    return summation

def compute_correlation_6(sa, dx0, dx1, dx2, dx3, dx4, dx5, ly0, ly1, ly2, ly3, ly4, ly5):
    """
    Compute the correlation for the middle part of DL distinguishers of Ascon
    when it involves 6 S-boxes with truncated input differences (dx0, dx1, dx2, dx3, dx4, dx5)
    and fixed output linear masks (ly0, ly1, ly2, ly3, ly4, ly5).
    """

    dlct = sa.differential_linear_connectivity_table()
    DX0 = sa.truncated_to_binvectors(dx0)
    DX1 = sa.truncated_to_binvectors(dx1)
    DX2 = sa.truncated_to_binvectors(dx2)
    DX3 = sa.truncated_to_binvectors(dx3)
    DX4 = sa.truncated_to_binvectors(dx4)
    DX5 = sa.truncated_to_binvectors(dx5)
    l0 = binlist_to_int(ly0, sa.n)
    l1 = binlist_to_int(ly1, sa.n)
    l2 = binlist_to_int(ly2, sa.n)
    l3 = binlist_to_int(ly3, sa.n)
    l4 = binlist_to_int(ly4, sa.n)
    l5 = binlist_to_int(ly5, sa.n)
    summation = 0
    for DX in itertools.product(*[DX0, DX1, DX2, DX3, DX4, DX5]):
        d0 = binlist_to_int(DX[0], sa.m)
        d1 = binlist_to_int(DX[1], sa.m)
        d2 = binlist_to_int(DX[2], sa.m)
        d3 = binlist_to_int(DX[3], sa.m)
        d4 = binlist_to_int(DX[4], sa.m)
        d5 = binlist_to_int(DX[5], sa.m)
        temp = dlct[d0][l0]*dlct[d1][l1]*dlct[d2][l2]*dlct[d3][l3]*dlct[d4][l4]*dlct[d5][l5]
        summation += temp
    return summation
    


if __name__ == '__main__':
    from sage.crypto.sboxes import Ascon as sb
    sa = SboxAnalyzer(sb)
    dlct = sa.differential_linear_connectivity_table()
    dx = [-1, -1, 0, -1, -1]
    ly = [1, 0, 0, 0, 0]
    output = compute_correlation_1(sa, dx, ly)
    print(output)

    # dx0 = [0, 1, -1, -1, -1]
    # dx1 = [1, 1, 1, 1, 1]
    # dx2 = [0, 0, 0, 0, 0]
    # ly0 = [0, 0, 1, 0, 0]
    # ly1 = ly2 = ly0
    # output = compute_correlation_3(sa, dx0, dx1, dx2, ly0, ly1, ly2)
    # print(output)

    # dx0 = [0, 0, 0, 0, -1]
    # dx1 = [0, 0, 1, 0, 0]
    # dx2 = [-1, 0, 0, 0, 0]
    # dx3 = [-1, 1, 0, 0, 0]
    # dx4 = [0, 0, 0, 0, -1]
    # dx5 = [0, 0, 0, 0, -1]
    # ly0 = [0, 0, 0, 1, 0]
    # ly1 = [1, 0, 0, 1, 1]
    # ly2 = [0, 1, 0, 1, 1]
    # ly3 = [0, 0, 0, 1, 0]
    # ly4 = [0, 0, 0, 1, 0]
    # ly5 = [0, 0, 0, 0, 1]
    # output = compute_correlation_6(sa, dx0, dx1, dx2, dx3, dx4, dx5, ly0, ly1, ly2, ly3, ly4, ly5)
    # print(output)

    # dx = binlist_to_int([1, 0, 0, 0, 1], sa.m)
    # ly = binlist_to_int([0, 0, 1, 0, 0], sa.m)
    # output = dlct[dx][ly]
    # print(output)
    