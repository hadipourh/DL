#!/usr/env/bin python3
#-*- coding: UTF-8 -*-

"""
MIT License

Copyright (c) 2024 

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

import textwrap

class Draw(object):
    def __init__(self, dlobject, output_file_name="output.tex", attack_summary=""):
        self.result = dlobject.result
        self.RU = dlobject.RU
        self.RM = dlobject.RM
        self.RL = dlobject.RL
        self.RD = dlobject.RD
        self.halfblocksize = dlobject.halfblocksize
        self.blocksize = dlobject.blocksize
        self.upper_trail = dlobject.upper_trail
        self.lower_trail = dlobject.lower_trail
        self.attack_summary = attack_summary
        self.output_file_name = output_file_name
        self.fillcolor_up = {0: "zero", 1: "upperfix", -1: "upperunknown"}
        self.fillcolor_lo = {0: "zero", 1: "lowerfix", -1: "lowerunknown"}

    def draw_eu(self, r):
        """
        Paint EU
        """
        output = dict()
        output["rot0"] = ""
        output["rot5"] = ""
        output["rot1"] = ""
        output["left"] = ""
        output["right"] = ""
        output["key"] = ""
        output["xor"] = ""
        output["and"] = ""
        for i in range(self.halfblocksize):
            output["left"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_up[self.upper_trail["xl"][r][i]], i) if self.upper_trail["xl"][r][i] != 0 else ""
            output["right"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_up[self.upper_trail["xr"][r][i]], i) if self.upper_trail["xr"][r][i] != 0 else ""
            output["rot0"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_up[self.upper_trail["xl"][r][i]], i) if self.upper_trail["xl"][r][i] != 0 else ""
            output["rot5"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_up[self.upper_trail["xl"][r][i]], (i - 5)%self.halfblocksize) if self.upper_trail["xl"][r][i] != 0 else ""
            output["rot1"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_up[self.upper_trail["xl"][r][i]], (i - 1)%self.halfblocksize) if self.upper_trail["xl"][r][i] != 0 else ""
            output["and"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_up[self.result["and_output"][r][i]], i) if self.result["and_output"][r][i] != 0 else ""
            output["xor"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_up[self.result["xor_output"][r][i]], i) if self.result["xor_output"][r][i] != 0 else ""
        return output

    def draw_em(self, r):
        """
        Paint EM
        """

        output = dict()
        output["left"] = ""
        output["right"] = ""
        for i in range(self.halfblocksize):
            output["left"] += "\TFill[{0}]{{s{1}}}".format(self.fillcolor_up[self.upper_trail["xl"][self.RU + r][i]], i) if self.upper_trail["xl"][self.RU + r][i] != 0 else ""
            output["right"] += "\TFill[{0}]{{s{1}}}".format(self.fillcolor_up[self.upper_trail["xr"][self.RU + r][i]], i) if self.upper_trail["xr"][self.RU + r][i] != 0 else ""

            output["left"] += "\BFill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.lower_trail["xl"][r][i]], i) if self.lower_trail["xl"][r][i] != 0 else ""
            output["right"] += "\BFill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.lower_trail["xr"][r][i]], i) if self.lower_trail["xr"][r][i] != 0 else ""
        return output

    def draw_el(self, r):
        """
        Paint EL
        """
        output = dict()
        output["rot0"] = ""
        output["rot5"] = ""
        output["rot1"] = ""
        output["left"] = ""
        output["right"] = ""
        output["key"] = ""
        output["xor"] = ""
        output["and"] = ""
        for i in range(self.halfblocksize):
            output["left"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.lower_trail["xl"][self.RM + r][i]], i) if self.lower_trail["xl"][self.RM + r][i] != 0 else ""
            output["right"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.lower_trail["xr"][self.RM + r][i]], i) if self.lower_trail["xr"][self.RM + r][i] != 0 else ""
            output["rot1"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.result["m0"][r][i]], (i - 1)%self.halfblocksize) if self.result["m0"][r][i] != 0 else ""
            output["rot5"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.result["m1"][r][i]], (i - 5)%self.halfblocksize) if self.result["m1"][r][i] != 0 else ""
            output["rot0"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.result["m2"][r][i]], i) if self.result["m2"][r][i] != 0 else ""
            output["and"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.lower_trail["xr"][self.RM + r][i]], i) if self.lower_trail["xr"][self.RM + r][i] != 0 else ""
            output["xor"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.lower_trail["xr"][self.RM + r][i]], i) if self.lower_trail["xr"][self.RM + r][i] != 0 else ""
        return output

    def draw_final_ed(self, r):
        """
        Paint the final round
        """
        output = dict()

        output["rot0"] = ""
        output["rot5"] = ""
        output["rot1"] = ""
        output["left"] = ""
        output["right"] = ""
        output["key"] = ""
        output["xor"] = ""
        output["and"] = ""

        for i in range(self.halfblocksize):
            if self.RL == 0:
                output["left"] += "\BFill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.lower_trail["xl"][self.RM + r][i]], i) if self.lower_trail["xl"][self.RM + r][i] != 0 else ""
                output["right"] += "\BFill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.lower_trail["xr"][self.RM + r][i]], i) if self.lower_trail["xr"][self.RM + r][i] != 0 else ""

                output["left"] += "\TFill[{0}]{{s{1}}}".format(self.fillcolor_up[self.upper_trail["xl"][self.RU + self.RM + r][i]], i) if self.upper_trail["xl"][self.RU + self.RM + r][i] != 0 else ""
                output["right"] += "\TFill[{0}]{{s{1}}}".format(self.fillcolor_up[self.upper_trail["xr"][self.RU + self.RM + r][i]], i) if self.upper_trail["xr"][self.RU + self.RM + r][i] != 0 else ""
            else:
                output["left"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.lower_trail["xl"][self.RM + r][i]], i) if self.lower_trail["xl"][self.RM + r][i] != 0 else ""
                output["right"] += "\Fill[{0}]{{s{1}}}".format(self.fillcolor_lo[self.lower_trail["xr"][self.RM + r][i]], i) if self.lower_trail["xr"][self.RM + r][i] != 0 else ""

        return output

    def generate_distinguisher_shape(self):
        """
        Draw the figure of the Intergal distinguisher
        """
        contents = ""
        contents += textwrap.dedent(r"""
            \documentclass[varwidth=100cm]{standalone}
            \usepackage{tikz}
            %%% macros %%%
            \usepackage{simeck}
            \usepackage{comment}
            \simeckshowkeyfalse
            \simeckcompacttrue
            \begin{document}
            \begin{figure}
                \begin{tikzpicture}[>=latex,fillopts/.style={black},raster/.style={gray!50}]
                    \simeckcompactfalse
                    \SimeckInit["""+str(self.blocksize)+ """]""" + "\n")
        # draw ED
        for r in range(0, self.RU):
            state = self.draw_eu(r)
            contents += r"""
            \SimeckRound{""" + str(r+1) + """}
            {""" + state["left"] + """}
            {""" + state["right"] + """}
            {""" + state["rot0"] + """}
            {""" + state["rot5"] + """}
            {""" + state["rot1"] + """}
            {}
            {""" + state["and"] + """}
            {""" + state["xor"] + "}""" + "\n"

        for r in range(0, self.RM):
            state = self.draw_em(r)
            contents += r"""
            \SimeckRoundShort{""" + str(self.RU + r + 1) + """}
            {""" + state["left"] + """}
            {""" + state["right"] + """}""" + "\n"
        for r in range(0, self.RL):
            state = self.draw_el(r)
            contents += r"""
            \SimeckRound{""" + str(self.RU + self.RM + r + 1) + """}
            {""" + state["left"] + """}
            {""" + state["right"] + """}
            {""" + state["rot0"] + """}
            {""" + state["rot5"] + """}
            {""" + state["rot1"] + """}
            {}
            {""" + state["and"] + """}
            {""" + state["xor"] + "}""" + "\n"
        for r in range(self.RL, self.RL + 1):
            state = self.draw_final_ed(r)
            contents += r"""
            \SimeckFinal{""" + str(self.RD + 1) + """}
            {""" + state["left"] + """}
            {""" + state["right"] + "}\n" + "\n"
        contents += "\n"
        contents += r"""  \end{tikzpicture}""" + "\n"
        contents += r"""\caption{""" + str(self.RD) + r""" rounds of {\tt Simeck}-""" + str(self.blocksize) + "}\n"
        contents += r"""\end{figure}""" + "\n"
        contents += r"""\begin{comment}""" + "\n"
        contents += self.attack_summary + "\n"
        contents += r"""\end{comment}""" + "\n"
        contents += r"""\end{document}"""
        with open(self.output_file_name, "w") as output_file:
            output_file.write(contents)
