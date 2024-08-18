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


import sys
PRINT_DIFFERENCE_VALUES = False

def trim(docstring):
    if not docstring:
        return ''
    # Convert tabs to spaces (following the normal Python rules)
    # and split into a list of lines:
    lines = docstring.expandtabs().splitlines()
    # Determine minimum indentation (first line doesn't count):
    indent = sys.maxsize
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))
    # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < sys.maxsize:
        for line in lines[1:]:
            trimmed.append(line[indent:].rstrip())
    # Strip off trailing and leading blank lines:
    while trimmed and not trimmed[-1]:
        trimmed.pop()
    while trimmed and not trimmed[0]:
        trimmed.pop(0)
    # Return a single string:
    return '\n'.join(trimmed)

class DrawDL():
    """
    Draw the shape of DL attack
    """

    def __init__(self, dl_object, output_file_name="output.tex"):
        self.result = dl_object.result
        self.RU = dl_object.RU
        self.RM = dl_object.RM
        self.RL = dl_object.RL
        self.RD = dl_object.RD
        self.inv_permutation = [0, 1, 2, 3, 5, 6, 7, 4, 10, 11, 8, 9, 15, 12, 13, 14]
        self.tweakey_permutation = [9, 15, 8, 13, 10, 14, 12, 11, 0, 1, 2, 3, 4, 5, 6, 7]
        self.output_file_name = output_file_name
        self.upper_trail = dl_object.upper_trail
        self.lower_trail = dl_object.lower_trail
        self.attack_summary = dl_object.attack_summary

    def gen_round_tweakey_labels(self, round_number):
        """
        Generate the round tweakey labels
        """
        
        round_tweakey_state = list(range(16))
        for r in range(round_number):
            round_tweakey_state = [self.tweakey_permutation[i] for i in round_tweakey_state]
        text = ""
        if not PRINT_DIFFERENCE_VALUES:
            for i in range(8):
                text += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, hex(round_tweakey_state[i])[2:])
        return text
    
    def draw_e0(self, r):
        """
        Paint E1 or E2
        """
        
        output = dict()
        output["before_sb"] = ""              
        output["after_sb"] = ""
        output["after_addtk"] = ""
        output["after_sr"] = ""
        output["subtweakey"] = ""
        output["after_mix_columns"] = ""
        for i in range(16):
            if self.result["DXU"][r][i] == 1:
                output["before_sb"] += "\FillCell[upperactive]{{s{0}}}".format(i)                
                output["after_sb"] += "\FillCell[upperactive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["before_sb"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.upper_trail[f"x_{r}"][i])
                    output["after_sb"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.upper_trail[f"y_{r}"][i])
            if self.result["DZU"][r][i] == 1:
                output["after_addtk"] += "\FillCell[upperactive]{{s{0}}}".format(i)        
                output["after_sr"] += "\FillCell[upperactive]{{s{0}}}".format(self.inv_permutation[i])
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["after_addtk"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.upper_trail[f"z_{r}"][i])
                    output["after_sr"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(self.inv_permutation[i], self.upper_trail[f"z_{r}"][i])                
            if self.result["DXU"][r + 1][i] == 1:
                if r + 1 == self.RU:
                    output["after_mix_columns"] += "\TFill[upperactive]{{s{0}}}".format(i)                    
                else:
                    output["after_mix_columns"] += "\FillCell[upperactive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["after_mix_columns"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.upper_trail[f"x_{r + 1}"][i])
            if self.result["DSTKU"][r][i] == 1 and i <= 7:
                output["subtweakey"] += "\FillCell[upperactive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["subtweakey"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.upper_trail[f"tk_{r}"][i])
        return output

    def draw_em(self, r):
        """
        Paint Em
        """
        
        output = dict()
        output["before_sb"] = ""              
        output["after_sb"] = ""
        output["after_addtk"] = ""
        output["after_sr"] = ""
        output["subtweakey"] = ""
        output["after_mix_columns"] = ""
        for i in range(16):
            if self.result["DXU"][r + self.RU][i] == 1:
                output["before_sb"] += "\TFill[upperactive]{{s{0}}}".format(i)
                output["after_sb"] += "\TFill[upperactive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["before_sb"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.upper_trail[f"x_{r + self.RU}"][i])
                    output["after_sb"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.upper_trail[f"y_{r + self.RU}"][i])
            if self.result["DZU"][r + self.RU][i] == 1:
                output["after_addtk"] += "\TFill[upperactive]{{s{0}}}".format(i)                
                output["after_sr"] += "\TFill[upperactive]{{s{0}}}".format(self.inv_permutation[i])
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["after_addtk"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.upper_trail[f"z_{r + self.RU}"][i])
                    output["after_sr"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(self.inv_permutation[i], self.upper_trail[f"z_{r + self.RU}"][i])
            if self.result["DXU"][r + 1 + self.RU][i] == 1:
                output["after_mix_columns"] += "\TFill[upperactive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["after_mix_columns"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.upper_trail[f"x_{r + 1 + self.RU}"][i])
            if self.result["DSTKU"][r + self.RU][i] == 1 and i <= 7:
                output["subtweakey"] += "\TFill[upperactive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["subtweakey"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.upper_trail[f"tk_{r + self.RU}"][i])

        for i in range(16):
            if self.result["DXL"][r][i] == 1:
                output["before_sb"] += "\BFill[loweractive]{{s{0}}}".format(i)
                output["after_sb"] += "\BFill[loweractive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["before_sb"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.lower_trail[f"x_{r}"][i])
                    output["after_sb"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.lower_trail[f"y_{r}"][i])
            if self.result["DZL"][r][i] == 1:
                output["after_addtk"] += "\BFill[loweractive]{{s{0}}}".format(i)                
                output["after_sr"] += "\BFill[loweractive]{{s{0}}}".format(self.inv_permutation[i])
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["after_addtk"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.lower_trail[f"z_{r}"][i])
                    output["after_sr"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(self.inv_permutation[i], self.lower_trail[f"z_{r}"][i])
            if self.result["DXL"][r + 1][i] == 1:
                output["after_mix_columns"] += "\BFill[loweractive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["after_mix_columns"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.lower_trail[f"x_{r + 1}"][i])
            if self.result["DSTKL"][r][i] == 1 and i <= 7:
                output["subtweakey"] += "\BFill[loweractive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["subtweakey"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.lower_trail[f"tk_{r}"][i])
        return output

    def draw_e1(self, r):
        """
        Paint E1
        """
        
        output = dict()
        output["before_sb"] = ""              
        output["after_sb"] = ""
        output["after_addtk"] = ""
        output["after_sr"] = ""
        output["subtweakey"] = ""
        output["after_mix_columns"] = ""
        for i in range(16):
            if self.result["DXL"][r + self.RM][i] == 1:
                output["before_sb"] += "\FillCell[loweractive]{{s{0}}}".format(i)
                output["after_sb"] += "\FillCell[loweractive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["before_sb"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.lower_trail[f"x_{r + self.RM}"][i])
                    output["after_sb"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.lower_trail[f"y_{r + self.RM}"][i])
            if self.result["DZL"][r + self.RM][i] == 1:
                output["after_addtk"] += "\FillCell[loweractive]{{s{0}}}".format(i)                
                output["after_sr"] += "\FillCell[loweractive]{{s{0}}}".format(self.inv_permutation[i])
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["after_addtk"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.lower_trail[f"z_{r + self.RM}"][i])
                    output["after_sr"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(self.inv_permutation[i], self.lower_trail[f"z_{r + self.RM}"][i])
            if self.result["DXL"][r + 1 + self.RM][i] == 1:
                output["after_mix_columns"] += "\FillCell[loweractive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["after_mix_columns"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.lower_trail[f"x_{r + 1 + self.RM}"][i])
            if self.result["DSTKL"][r + self.RM][i] == 1 and i <= 7:
                output["subtweakey"] += "\FillCell[loweractive]{{s{0}}}".format(i)
                if PRINT_DIFFERENCE_VALUES:
                    # difference values
                    output["subtweakey"] += "\Cell{{s{0}}}{{\\texttt{{{1}}}}}".format(i, self.lower_trail[f"tk_{r + self.RM}"][i])
        return output


    def generate_attack_shape(self):
        """
        Draw the figure of the Rectangle distinguisher
        """

        contents = ""
        # head lines
        contents += trim(r"""
                    \documentclass[varwidth=50cm]{standalone}
                    \usepackage{skinnyzero}
                    \usepackage{comment}
                    \begin{document}
                    \begin{tikzpicture}
                    \SkinnyInit{}{}{}{} % init coordinates, print labels""") + "\n\n"
        # draw E0
        for r in range(self.RU):
            state = self.draw_e0(r)
            state["subtweakey"] += self.gen_round_tweakey_labels(r)
            contents += trim(r"""
            \SkinnyRoundTK[""" + str(r) + """]
                          {""" + state["before_sb"] + r"""} % state (input)
                          {""" + state["subtweakey"] + r"""} % tk[1]
                          {""" + r"""} % tk[2]
                          {""" + r"""} % tk[3]
                          {""" + state["after_sb"] + """} % state (after subcells)
                          {""" + state["after_addtk"] + r"""} % state (after addtweakey)
                          {""" + state["after_sr"] + r"""} % state (after shiftrows)""") + "\n\n"
            if (r) % 2 == 1:
                contents += trim(r"""\SkinnyNewLine[""" + str(r + 1) + r"""]{""") + state["after_mix_columns"] + r"""} % state (after mixcols)""" + "\n"
            
        
        # draw Em
        for r in range(self.RM):
            state = self.draw_em(r)
            state["subtweakey"] += self.gen_round_tweakey_labels(r + self.RU)
            contents += trim(r"""
            \SkinnyRoundTK[""" + str(r + self.RU) + """]
                          {""" + state["before_sb"] + r"""} % state (input)
                          {""" + state["subtweakey"] + r"""} % tk[1]
                          {""" + r"""} % tk[2]
                          {""" + r"""} % tk[3]
                          {""" + state["after_sb"] + """} % state (after subcells)
                          {""" + state["after_addtk"] + r"""} % state (after addtweakey)
                          {""" + state["after_sr"] + r"""} % state (after shiftrows)""") + "\n\n"
            if self.RL == 0 and r == (self.RM - 1):
                contents += trim(r"""\SkinnyFin[""" + str(r + self.RU + 1) + r"""]{""") + state["after_mix_columns"] + r"""} % state (after mixcols)""" + "\n"
            elif (r + self.RU) % 2 == 1:
                contents += trim(r"""\SkinnyNewLine[""" + str(r + self.RU + 1) + r"""]{""") + state["after_mix_columns"] + r"""} % state (after mixcols)""" + "\n"
        
        # draw E1
        for r in range(self.RL):
            state = self.draw_e1(r)
            state["subtweakey"] += self.gen_round_tweakey_labels(r + self.RU + self.RM)
            contents += trim(r"""
            \SkinnyRoundTK[""" + str(r + self.RU + self.RM) + """]
                          {""" + state["before_sb"] + r"""} % state (input)
                          {""" + state["subtweakey"] + r"""} % tk[1]
                          {""" + r"""} % tk[2]
                          {""" + r"""} % tk[3]
                          {""" + state["after_sb"] + """} % state (after subcells)
                          {""" + state["after_addtk"] + r"""} % state (after addtweakey)
                          {""" + state["after_sr"] + r"""} % state (after shiftrows)""") + "\n\n"
            if r == self.RL - 1:
                contents += trim(r"""\SkinnyFin[""" + str(r + self.RU + self.RM + 1) + r"""]{""") + state["after_mix_columns"] + r"""} % state (after mixcols)""" + "\n"
            elif (r + self.RU + self.RM) % 2 == 1:
                contents += trim(r"""\SkinnyNewLine[""" + str(r + self.RU + self.RM + 1) + r"""]{""") + state["after_mix_columns"] + r"""} % state (after mixcols)""" + "\n"
        contents += r"""\begin{comment}""" + "\n"
        contents += self.attack_summary
        contents += r"""\end{comment}""" + "\n"
        contents += r"""\end{tikzpicture}""" + "\n"
        contents += trim(r"""\end{document}""")
        with open(self.output_file_name, "w") as output_file:
            output_file.write(contents)