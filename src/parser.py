#!/usr/bin/env python

import sys
from termcolor import colored
import colorama
import re
from collections import Counter
import warnings

from typing import Dict, Tuple, List, Iterable, Mapping


FormulaType = str

colorama.init()

ATOM_REGEX = r'([A-Z][a-z]*)(\d*)'
OPENERS = r'({['
CLOSERS = r')}]'


def is_balanced(formula: FormulaType) -> bool:
    """Check if all sort of brackets come in pairs."""
    # Very naive check, just here because you always need some input checking
    c = Counter(formula)
    return c['['] == c[']'] and c['{'] == c['}'] and c['('] == c[')']


def is_real_elements(formula: FormulaType) -> bool:
    """Some basic string chacks to be sure there are some elements present"""

    elem_first_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    elem_sec_letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'l', 'm', 'n', 'o', 'r', 's', 't', 'u', 'v', 'y']
    elem_letters = elem_first_letters + elem_sec_letters
    # non_letter = ['J', 'Q', 'j', 'k', 'p', 'q', 'w', 'x', 'z']
    allowed_symbols = ['[', ']', '(', ')', '{', '}', '-', '–']

    # Check if formula string contains only elemental letters, or allowed symbols.
    for letter in formula:
        if letter.isalpha():
            if letter in elem_letters:
                continue
            else:
                return False
        elif not letter.isnumeric() and letter not in allowed_symbols:
            return False

    # check there is at least one element
    if not any(letter.isupper() for letter in formula):
        return False

    # check first character (if not a bracket) is an element
    if formula[0].isalpha() and not formula[0].isupper() and formula[0] not in elem_first_letters:
        return False

    # check first character is not a number
    if formula[0].isnumeric():
        return False

    # if all checks pass, return true
    return True


def _dictify(tuples: Iterable[Tuple[str, str]]) -> Dict[str, int]:
    """Transform tuples of tuples to a dict of atoms."""
    res: Dict[str, int] = {}
    for atom, n in tuples:
        try:
            res[atom] += int(n or 1)
        except KeyError:
            res[atom] = int(n or 1)
    return res


def _fuse(mol1: Mapping[str, int], mol2: Mapping[str, int], w: int = 1) -> Dict[str, int]:
    """
    Fuse 2 dicts representing molecules. Return a new dict.
    This fusion does not follow the laws of physics.
    """
    return {atom: (mol1.get(atom, 0) + mol2.get(atom, 0)) * w for atom in set(mol1) | set(mol2)}


def _parse(formula: FormulaType) -> Tuple[Dict[str, int], int]:
    """
    Return the molecule dict and length of parsed part.
    Recurse on opening brackets to parse the subpart and
    return on closing ones because it is the end of said subpart.
    """
    q: List[str] = []
    mol: Dict[str, int] = {}
    i = 0

    while i < len(formula):
        # Using a classic loop allow for manipulating the cursor
        token = formula[i]

        if token in CLOSERS:
            # Check for an index for this part
            m = re.match(r'\d+', formula[i+1:])
            if m:
                weight = int(m.group(0))
                i += len(m.group(0))
            else:
                weight = 1

            submol = _dictify(re.findall(ATOM_REGEX, ''.join(q)))
            return _fuse(mol, submol, weight), i

        elif token in OPENERS:
            submol, sublength = _parse(formula[i+1:])
            mol = _fuse(mol, submol)
            # skip the already read submol
            i += sublength + 1
        else:
            q.append(token)

        i += 1

    # Fuse in all that's left at base level
    return _fuse(mol, _dictify(re.findall(ATOM_REGEX, ''.join(q)))), i


def parse_formula(formula: FormulaType) -> Dict[str, int]:
    """Parse the formula and return a dict with occurences of each atom."""

    warn_text = (colored('Check your formula [', 'yellow', attrs=['bold']) +
                 colored(formula, 'red', attrs=['bold']) +
                 colored('] - are all the element symbols correct?\n', 'yellow', attrs=['bold']))

    if not is_real_elements(formula):
        warnings.warn(warn_text, category=UserWarning, stacklevel=4)

    try:
        if not is_balanced(formula):
            raise ValueError(colored("Watch your brackets ![{]$[&?)]}!]", 'red', attrs=['bold']))
    except ValueError as v_err:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        #fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        #print(exc_type, fname, exc_tb.tb_lineno)
        if exc_type:
            print(f'Error: {exc_type.__name__}: {v_err}\n')
        sys.exit()

    parsed = _parse(formula)[0]

    for element in parsed.keys():
        if len(element) > 2:
            warnings.warn(warn_text, category=UserWarning, stacklevel=4)
            # raise ValueError(colored("Some elemental symbols don't make sense", 'red'))

    return parsed


def parse_print_from_input() -> None:
    formula = str(input("\nEnter chemical formula: \n"))
    print("")
    res = parse_formula(formula)
    colored_res = colored(str(res), 'green', attrs=['bold'])
    print(f'Dict of elements from given formula: {colored_res}\n')


if __name__ == "__main__":
    parse_print_from_input()
