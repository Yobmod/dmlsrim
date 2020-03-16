import re

tdata = """===============Target material=======================

           Layer 12     Ceria 1.33 kg/mol """
reg = re.compile(r'''(=+[a-zA-Z]+\ [a-zA-Z]+=+\n\s*)    # 1 header
                 (\b[Ll]ayer\b\s+\d+)\s+                # 2 layer and number
                 (.+)\s+                                   # 3 layer name
                 ([0-9]\.[0-9]+)\                          # 4 layer percentage
                 (?:kg/mol)''',                         # ?:  non captured
                 re.DOTALL | re.VERBOSE)
# \b\w+(?<!s)\b  any word not ending in s

match_target = reg.search(tdata)
match_target = re.search(R'\=+\n*\t*\r*Layer\s+\d+\s', tdata, re.DOTALL)
#match_target = re.search(r'(?<=====\r\n)Layer\s+\d+\s+:.*?(?=====)', f.read(), re.DOTALL)
# ## ?<=  positive lookbehind, ?= positive look ahead

if match_target:
    print(match_target.group(3))
    print(f"span = {match_target.span()}")
else:
    print("target not found")

# Valid double_regex 4, 4.0, 4.0e100
double_regex = r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?'
symbol_regex = r'[A-Z][a-z]?'
int_regex = r'[+-]?\d+'


def _read_target(output: bytes) -> None:
    match_target = re.search(br'(?<=====\r\n)Layer\s+\d+\s+:.*?(?=====)', output, re.DOTALL)
    if match_target:
        print(match_target.group(0))
        layer_regex = (
            r'Layer\s+(?P<i>\d+)\s+:\s+(.+)\r\n'
            r'Layer Width\s+=\s+({0})\s+A\s+;\r\n'
            r'\s+Layer #\s+(?P=i)- Density = ({0}) atoms/cm3 = ({0}) g/cm3\r\n'
            r'((?:\s+Layer #\s+(?P=i)-\s+{1}\s+=\s+{0}\s+Atomic Percent = {0}\s+Mass Percent\r\n)+)'
        ).format(double_regex, symbol_regex)
        layers = re.findall(layer_regex.encode('utf-8'), match_target.group(0))
        if layers:
            element_regex = (
                r'\s+Layer #\s+(\d+)-\s+({1})\s+=\s+({0})\s+Atomic Percent = ({0})\s+Mass Percent\r\n'
            ).format(double_regex, symbol_regex)
            element_regex_str = element_regex.encode()

            layers_elements = []
            for layer in layers:
                # We know that elements will match
                layers_elements.append(re.findall(element_regex_str, layer[5]))


R"""

[]      define charachter class
.       Matches Any char except \n (or including \n using re.DOTALL arg)
?       Match the preceeding char or class one or zero times
+       Match the preceeding char or class repeated any number of times, excluding zero
*       Match the preceeding char or class repeated any number of times, including zero
^       only match following word if beginnning of line (ignoring whitespace unless re.MULTILINE)
$       only match preceeding word if end of a line (eg at \n)
\b      WORD boundary
\B      Not a word boundary
\t      Tab
\n      Newline
\r      Return (windows)
\f
\v
\d      Matches any decimal digit; this is equivalent to the class [0-9].
\D      Matches any non-digit character; this is equivalent to the class [^0-9].
\s      Matches any whitespace character; this is equivalent to the class [ \t\n\r\f\v].
\S     Matches any non-whitespace character; this is equivalent to the class [^ \t\n\r\f\v].
\w     Matches any alphanumeric character; this is equivalent to the class [a-zA-Z0-9_].
\W     Matches any non-alphanumeric character; this is equivalent to the class [^a-zA-Z0-9_].

layer\s
"""
