#!/usr/bin/env python3

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

Disclaimer: We acknowledge that the WARP block cipher doesn't adhere to statistical assumptions 
in differential analysis, such as the random sub-key assumption
or Markov cipher assumption. The tool's primary function is to find some bounds
for the security of TWINE against differential and differential-linear cryptanalysis.
"""

pi = [31, 6, 29, 14, 1, 12, 21, 8, 27, 2, 3, 0, 25, 4, 23, 10, 15, 22, 13, 30, 17, 28, 5, 24, 11, 18, 19, 16, 9, 20, 7, 26]

def tex_init(options=""):
    tex_content = r"% generated by plotdistinguisher.py" + options + "\n" + \
             r"\documentclass[preview]{standalone}" + "\n" + \
            r"\usepackage{comment}" + "\n" + \
             r"\usepackage{tugcolors}" + "\n" + \
             r"\usepackage{warp}" + "\n" + \
             r"\usepackage[margin=1in]{geometry}" + "\n" + \
             r"%\tikzset{warpfig/.append style={white}}" + "\n" + \
             r"\begin{document}" + "\n" + \
             r"%\begin{figure}[htp!]" + "\n"
    return tex_content

def tikz_mark_input_bits(bits, color="tugred"):
    tex_content = r"\begin{tikzpicture}[warpfig]" + "\n" + \
                  r"\foreach \z in {""" + ",".join([str(bit) for bit in bits]) + \
                  r"} { \fill[" + color + r"] (.25*\z,.75) circle[radius=3pt]; }" + "\n" + \
                  r"   \foreach \z[evaluate=\z as \zf using int(4*\z)] in {0,...,31} {" + "\n" + \
                  r"       \draw[gray] (\z,0) node[above] {\tiny\zf};" + "\n" + \
                  r"       \foreach \zb in {0,...,3} { \draw[gray] (\z+.25*\zb,0) -- +(0,-3pt); }" + "\n" + \
                  r"}" + "\n" + \
                  r"\end{tikzpicture}" + "\n"
    return tex_content

def tikz_mark_output_bits(bits, color="tugblue"):
    tex_content = r"\begin{tikzpicture}[warpfig]" + "\n" + \
                  r"\foreach \z[evaluate=\z as \zf using int(4*\z)] in {0,...,31} {" + "\n" + \
                  r"\draw[gray] (\z,-3pt) node[below] {\tiny\zf};" + "\n" + \
                  r"\foreach \zb in {0,...,3} { \draw[gray] (\z+.25*\zb,0) -- +(0,-3pt); }" + "\n" + \
                  r"}" + "\n" + \
                  r"\foreach \z in {" + ",".join([str(bit) for bit in bits]) + r"} { \fill[" + color + \
                  r"] (.25*\z,-.75) circle[radius=3pt]; }" + "\n" + \
                  r"\end{tikzpicture}" + "\n"
    return tex_content

def tex_diff_trail(trail, markpattern, direction="->"):
    tex_content = ""
    nrounds  = trail["nrounds"]
    for r in range(nrounds):
        tex_content += r"\warproundwokey{%s" + "\n"
        active_input_branches = [n for n in range(32) if trail[f"x_{r}"][n] != "0"]
        active_even_branches = [n for n in active_input_branches if n%2 == 0]
        active_odd_branches = [n for n in active_input_branches if n%2 != 0]
        tex_content += r"\markevenbranches{" + \
                        ",".join([f"{n}/{pi[n]}" for n in active_even_branches])+ \
                       r"}{" + markpattern + r"}{" + direction + r"}" + "\n"
        tex_content += r"\markoddbranchbeforexor{" + \
                        ",".join([f"{n}" for n in active_odd_branches]) + \
                        r"}{" + markpattern + r"}" + "\n"
        active_after_xor = [n for n in range(32) if trail[f"x_{r+1}"][pi[n]] != "0" and n%2 == 1]
        tex_content += r"\markoddbranchafterxor{" + \
                       ",".join([f"{n}/{pi[n]}" for n in active_after_xor]) + \
                       r"}{" + markpattern + r"}{" + direction + r"}" + "\n"
        tex_content += "}\n"
    return tex_content

def tex_diff_lower_trail(trail, upper_crossing_difference, markpattern, direction="->"):
    tex_content = ""
    nrounds  = trail["nrounds"]
    for r in range(nrounds):
        tex_content += r"\warproundwokey{%s" + "\n"
        active_input_branches = [n for n in range(32) if trail[f"x_{r}"][n] != "0"]
        active_even_branches = [n for n in active_input_branches if n%2 == 0]
        active_odd_branches = [n for n in active_input_branches if n%2 != 0]
        if r == 0:
            tex_content += r"\markuppercrossingdifferences{" + \
                           ",".join(upper_crossing_difference) + \
                           "}\n"
        tex_content += r"\markevenbranches{" + \
                        ",".join([f"{n}/{pi[n]}" for n in active_even_branches])+ \
                       r"}{" + markpattern + r"}{" + direction + r"}" + "\n"
        tex_content += r"\markoddbranchbeforexor{" + \
                        ",".join([f"{n}" for n in active_odd_branches]) + \
                        r"}{" + markpattern + r"}" + "\n"
        active_after_xor = [n for n in range(32) if trail[f"x_{r+1}"][pi[n]] != "0" and n%2 == 1]
        tex_content += r"\markoddbranchafterxor{" + \
                       ",".join([f"{n}/{pi[n]}" for n in active_after_xor]) + \
                       r"}{" + markpattern + r"}{" + direction + r"}" + "\n"
        if r == nrounds - 1:
            tex_content += r"\markoutputdiff{" + \
                           ",".join([str(i) for i in range(32) if trail[f"x_{nrounds}"][i] != "0"]) + \
                           "}\n"
        tex_content += "}\n"
    return tex_content

def tex_lin_lower_trail(trail, upper_crossing_difference, markpattern, direction="->"):
    tex_content = ""
    nrounds  = trail["nrounds"]
    for r in range(nrounds):
        tex_content += r"\warproundwokey{%s" + "\n"
        active_input_branches = [n for n in range(32) if trail[f"x_{r}"][n] != "0"]
        active_even_branches = [n for n in active_input_branches if n%2 == 0]
        active_odd_branches = [n for n in active_input_branches if n%2 != 0]
        if r == 0:
            tex_content += r"\markuppercrossingdifferences{" + \
                           ",".join(upper_crossing_difference) + \
                           "}\n"
        tex_content += r"\markevenbranchbeforetee{" + \
                        ",".join([f"{n}" for n in active_even_branches])+ \
                       r"}{" + markpattern + r"}{" + direction + r"}" + "\n"
        tex_content += r"\markoddbranch{" + \
                        ",".join([f"{n}/{pi[n]}" for n in active_odd_branches]) + \
                        r"}{" + markpattern + r"}{" + direction + r"}" + "\n"
        active_after_tee = [n for n in range(32) if trail[f"x_{r+1}"][pi[n]] != "0" and n%2 == 0]
        tex_content += r"\markevenbranchaftertee{" + \
                       ",".join([f"{n}/{pi[n]}" for n in active_after_tee]) + \
                       r"}{" + markpattern + r"}{" + direction + r"}" + "\n"
        if r == nrounds - 1:
            tex_content += r"\markoutputdiff{" + \
                           ",".join([str(i) for i in range(32) if trail[f"x_{nrounds}"][i] != "0"]) + \
                           "}\n"
        tex_content += "}\n"
    return tex_content

def tex_middle(upper_trail, midd_trail, lower_trail, RU, RM, RL):
    tex_content = ""
    for r in range(RM):
        ur = r + RU
        tex_content += r"\warproundwokey{%s" + "\n"
        active_branches_upper = [n for n in range(32) if upper_trail[f"x_{ur}"][n] != "0"]
        active_branches_lower = [n for n in range(32) if lower_trail[f"x_{r}"][n] != "0"]
        active_even_branches_upper = [n for n in active_branches_upper if n%2 == 0]
        active_even_branches_lower = [n for n in active_branches_lower if n%2 == 0]
        active_odd_branches_upper = [n for n in active_branches_upper if n%2 != 0]
        active_odd_branches_lower = [n for n in active_branches_lower if n%2 != 0]
        tex_content += r"\markevenbranches{" + \
                       ",".join([f"{n}/{pi[n]}" for n in active_even_branches_upper])+ \
                       r"}{markmidupperpath}{->}" + "\n"
        tex_content += r"\markoddbranchbeforexor{" + \
                        ",".join([f"{n}" for n in active_odd_branches_upper]) + \
                        r"}{markmidupperpath}" + "\n"
        active_after_xor = [n for n in range(32) if upper_trail[f"x_{ur+1}"][pi[n]] != "0" and n%2 == 1]
        tex_content += r"\markoddbranchafterxor{" + \
                       ",".join([f"{n}/{pi[n]}" for n in active_after_xor]) + \
                       r"}{markmidupperpath}{->}" + "\n"

        tex_content += r"\markevenbranchbeforetee{" + \
                        ",".join([f"{n}" for n in active_even_branches_lower])+ \
                       r"}{markmidlowerpath}{<-}" + "\n"
        tex_content += r"\markoddbranch{" + \
                        ",".join([f"{n}/{pi[n]}" for n in active_odd_branches_lower]) + \
                        r"}{markmidlowerpath}{<-}" + "\n"
        active_after_tee = [n for n in range(32) if lower_trail[f"x_{r+1}"][pi[n]] != "0" and n%2 == 0]
        tex_content += r"\markevenbranchaftertee{" + \
                       ",".join([f"{n}/{pi[n]}" for n in active_after_tee]) + \
                       r"}{markmidlowerpath}{<-}" + "\n"
        tex_content += r"\markcommonactivesboxes{" + \
                       ",".join([str(i//2) for i in range(32) if midd_trail[f"s_{r}"][i] == "1"]) + \
                       "}\n"
        if RL == 0 and r == RM - 1:
            tex_content += r"""	\foreach \z in {""" + \
            ",".join([str(i) for i in range(32) if lower_trail[f"x_{RM + RL}"][i] != "0"]) + \
            r"""} {""" + "\n" + \
            r"""\node[tugblue,fill,fill opacity=.2, below] at (\z,{-\diffusion_depth - 0.6}) {\phantom{$X_{\!00}$}};""" + \
	         "}\n" + \
            r"""	\foreach \z in {""" + \
            ",".join([str(i) for i in range(32) if upper_trail[f"x_{RU + RM}"][i] != "0"]) + \
            r"""} {""" + "\n" + \
            r"""\node[tugred,fill,fill opacity=.2, below] at (\z,{-\diffusion_depth - 0.6}) {\phantom{$X_{\!00}$}};""" + \
	         "}\n" + \
            r"""\foreach \z in {0, ..., 31} {""" + \
            r"""\draw (\z,{-\diffusion_depth - 0.6}) node[below] {$X_{\!\z}$};""" + "\n" + \
            "}"
        tex_content += "}\n"
    return tex_content

def tex_fin(rounds):
    tex_content = r"%\caption{Differential-linear distinguisher for " + str(rounds) + r" rounds of \texttt{WARP}.}" + "\n" + \
             r"%\label{fig:difflin_distinguisher}" + "\n" + \
             r"%\end{figure}" + "\n" + \
             r"\end{document}" + "\n"
    return tex_content