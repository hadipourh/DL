# !/usr/bin/env sage
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
from argparse import ArgumentParser, RawTextHelpFormatter

def doct_product(a, b, n=4):
    """
    Compute the dot product of two inetgers as a list of bits
    """

    output = 0
    for i in range(n):
        output ^= ((a >> i) & 1) * ((b >> i) & 1)
    return output

def compute_ddt4(sb):
    """
    compute DDT4
    """

    n = sb.input_size()
    ddt = sb.difference_distribution_table()
    ddt4 = [[0 for _ in range(2**n)] for _ in range(2**n)]
    for D1 in range(2**n):
        for D2 in range(2**n):
            for D3 in range(2**n):
                for D4 in range(2**n):
                    for D5 in range(2**n):
                        ddt4[D1][D5] += ddt[D1][D2] * ddt[D2][D3] * ddt[D3][D4] * ddt[D4][D5]
    return ddt4

def compute_correlation_13r_v0(sb, D=None, L=None):
    """
    Compute the correlation for 13-round DLD-v0 
    """
    
    from sboxanalyzer import SboxAnalyzer
    sa = SboxAnalyzer(sb)
    n = sa.input_size()
    ddlct = sa.double_differential_linear_connectivity_table()
    if D is not None and L is not None:
        corr = ddlct[D][L]/2**(2*n)
        return corr
    else: 
        for D in range(2**4):
            for L in range(2**4):
                corr = ddlct[D][L]/2**(2*n)
                if corr > 0:
                    sign = "+"
                elif corr < 0:
                    sign = "-"
                else:
                    sign = ""
                if corr != 0:
                    log2_abscorr = float(log(abs(corr), 2))
                    print("({}, {}): Corr = {}2^({:1.2f})".format(hex(D), hex(L), sign, log2_abscorr))
                else:
                    print("({}, {}): Corr = 0".format(hex(D), hex(L)))

def compute_correlation_13r_v1(sb, D=None, L=None):
    """
    Compute the correlation for 13-round DLD-v1 
    """
    
    from sboxanalyzer import SboxAnalyzer
    sa = SboxAnalyzer(sb)
    n = sa.input_size()
    tdlct = sa.triple_differential_linear_connectivity_table()
    if D is not None and L is not None:
        corr = tdlct[D][L]/2**(3*n)
        return corr**2
    else: 
        for D in range(2**4):
            for L in range(2**4):
                corr = tdlct[D][L]/2**(3*n)
                corr = corr**2
                if corr > 0:
                    sign = "+"
                elif corr < 0:
                    sign = "-"
                else:
                    sign = ""
                if corr != 0:
                    log2_abscorr = float(log(abs(corr), 2))
                    print("({}, {}): Corr = {}2^({:1.2f})".format(hex(D), hex(L), sign, log2_abscorr))
                else:
                    print("({}, {}): Corr = 0".format(hex(D), hex(L)))


def compute_correlation_13r_v2(sb, D=None, L=None):
    """
    Compute the correlation for 13-round DLD-v2 
    """

    from sboxanalyzer import SboxAnalyzer
    sa = SboxAnalyzer(sb)
    n = sa.input_size()

    corr = 0
    ddt4 = compute_ddt4(sa)
    ddt = sa.difference_distribution_table()
    dlct = sa.differential_linear_connectivity_table()
    if D is not None and L is not None:
        for D2 in range(2**n):
            for D6 in range(2**n):
                for D7 in range(2**n):
                    corr += ddt[D][D2] * ddt4[D2][D6] * ddt[D][D7] * dlct[D6 ^ D7][L] * dlct[D2][L]
        corr = corr/(2**(8*n))
        return corr
    else:
        for D in range(2**n):
            for L in range(2**n):
                for D2 in range(2**n):
                    for D6 in range(2**n):
                        for D7 in range(2**n):
                            corr += ddt[D][D2] * ddt4[D2][D6] * ddt[D][D7] * dlct[D6 ^ D7][L] * dlct[D2][L]
                corr = corr/(2**(8*n))            
                if corr > 0:
                    sign = "+"
                elif corr < 0:
                    sign = "-"
                else:
                    sign = ""
                if corr != 0:
                    log2_abscorr = float(log(abs(corr), 2))
                    print("({}, {}): Corr = {}2^({:1.2f})".format(hex(D), hex(L), sign, log2_abscorr))
                else:
                    print("({}, {}): Corr = 0".format(hex(D), hex(L)))

def loadparameters(args):
        """
        Get parameters from the argument list and inputfile.
        """
        # Load default values
        params = {"version": 0}

        # Override parameters if they are set on commandline

        if args.version:
            params["version"] = args.version[0]

        return params

def main():
    """
    Parse the arguments and start the request functionality with the provided
    parameters.
    """

    parser = ArgumentParser(description="This script computes the correlation of the 13-round DL distinguishers of TWINE in our paper: https://ia.cr/2024/255",
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('--version', nargs=1, type=int, default=[0],
                        help="Version of the 13-round DL distinguisher. Default is 0.")

    # Parse command line arguments and construct parameter list.
    args = parser.parse_args()
    params = loadparameters(args)
    from sage.crypto.sboxes import TWINE as sb
    if params["version"] == 0:
        compute_correlation_13r_v0(sb)
    elif params["version"] == 1:
        compute_correlation_13r_v1(sb)
    elif params["version"] == 2:
        compute_correlation_13r_v2(sb)
    else:
        raise ValueError("Invalid version number")

if __name__ == "__main__":
    main()