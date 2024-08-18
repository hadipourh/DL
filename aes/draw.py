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
    Draw the shape of a given differential-linear distinguisher
    """

    def __init__(self, dlobject, output_file_name="output.tex"):
        self.result = dlobject.result
        self.RU = dlobject.RU
        self.RM = dlobject.RM
        self.RL = dlobject.RL
        self.RD = dlobject.RD
        self.output_file_name = output_file_name
        self.upper_trail = dlobject.upper_trail
        self.lower_trail = dlobject.lower_trail
        self.attack_summary = dlobject.attack_summary        

    def draw_ed(self, r):
        """
        Paint the distinguisher in the r-th round
        """
        
        state = dict()
        state["input"] = ""
        state["round_key"] = ""
        state["after_addrk"] = ""
        state["after_sb"] = ""
        state["after_sr"] = ""
        state["output"] = ""
        if r < self.RU:
            for row in range(4):
                for col in range(4):
                    if self.upper_trail["x"][r][row][col] == 1:
                        state["input"] += "\TFill[upperactive]{{ss{0}{1}}}".format(row, col)
                        state["after_addrk"] += "\TFill[upperactive]{{ss{0}{1}}}".format(row, col)
                        state["after_sb"] += "\TFill[upperactive]{{ss{0}{1}}}".format(row, col)
                    if self.upper_trail["y"][r][row][col] == 1:
                        state["after_sr"] += "\TFill[upperactive]{{ss{0}{1}}}".format(row, col)
                    if self.upper_trail["x"][r + 1][row][col] == 1:
                        state["output"] += "\TFill[upperactive]{{ss{0}{1}}}".format(row, col)
        elif r < self.RU + self.RM:
            for row in range(4):
                for col in range(4):
                    if self.upper_trail["x"][r][row][col] == 1:
                        state["input"] += "\TFill[upperactive]{{ss{0}{1}}}".format(row, col)
                        state["after_addrk"] += "\TFill[upperactive]{{ss{0}{1}}}".format(row, col)
                        state["after_sb"] += "\TFill[upperactive]{{ss{0}{1}}}".format(row, col)
                    if self.upper_trail["y"][r][row][col] == 1:
                        state["after_sr"] += "\TFill[upperactive]{{ss{0}{1}}}".format(row, col)
                    if self.upper_trail["x"][r + 1][row][col] == 1:
                        state["output"] += "\TFill[upperactive]{{ss{0}{1}}}".format(row, col)
                    if self.lower_trail["x"][r - self.RU][row][col] == 1:
                        state["input"] += "\BFill[loweractive]{{ss{0}{1}}}".format(row, col)
                        state["after_addrk"] += "\BFill[loweractive]{{ss{0}{1}}}".format(row, col)
                        state["after_sb"] += "\BFill[loweractive]{{ss{0}{1}}}".format(row, col)
                    if self.lower_trail["y"][r - self.RU][row][col] == 1:
                        state["after_sr"] += "\BFill[loweractive]{{ss{0}{1}}}".format(row, col)
                    if self.lower_trail["x"][r - self.RU + 1][row][col] == 1:
                        state["output"] += "\BFill[loweractive]{{ss{0}{1}}}".format(row, col)
        else:
            for row in range(4):
                for col in range(4):
                    if self.lower_trail["x"][r - self.RU][row][col] == 1:
                        state["input"] += "\BFill[loweractive]{{ss{0}{1}}}".format(row, col)
                        state["after_addrk"] += "\BFill[loweractive]{{ss{0}{1}}}".format(row, col)
                        state["after_sb"] += "\BFill[loweractive]{{ss{0}{1}}}".format(row, col)
                    if self.lower_trail["y"][r - self.RU][row][col] == 1:
                        state["after_sr"] += "\BFill[loweractive]{{ss{0}{1}}}".format(row, col)
                    if self.lower_trail["x"][r - self.RU + 1][row][col] == 1:
                        state["output"] += "\BFill[loweractive]{{ss{0}{1}}}".format(row, col)
        return state

    def generate_distinguisher_shape(self):
        """
        Draw the figure of the distinguisher
        """
        
        contents = ""
        # Head lines

        contents += trim(r"""
                    \documentclass{standalone}
                    \usepackage{aes}
                    \usepackage{comment}
                    \begin{document}
                    %\begin{figure}
                    %\centering
                    \begin{tikzpicture}
                    \AesInit""") + "\n\n"
        for r in range(self.RD):
            state = self.draw_ed(r)
            contents += trim(r"""
                    \AesRound[""" + str(r) + """]
                    {""" + state["input"] + """} % state input
                    {""" + state["round_key"]+ """} % round key
                    {""" + state["after_addrk"] + """} % state after AK
                    {""" + state["after_sb"] + """} % state after SB
                    {""" + state["after_sr"] + """} % state after SR""") + "\n\n"
            if r == self.RD - 1:
                contents += trim(r"""\AesFin[""" + str(r + 1) + r"""]{""") + state["output"] + r"""}""" + "\n"
            elif r % 2 == 1:
                contents += trim(r"""\AesNewLine[""" + str(r + 1) + r"""]{""") + state["output"] + r"""}""" + "\n"
        contents += trim(r"""\end{tikzpicture}""") + "\n"
        contents += trim(r"""%\end{figure}""") + "\n"
        contents += trim(r"""\begin{comment}""") + "\n"
        contents += self.attack_summary
        contents += trim(r"""\end{comment}""") + "\n"
        contents += trim(r"""\end{document}""")
        with open(self.output_file_name, "w") as output_file:
            output_file.write(contents)