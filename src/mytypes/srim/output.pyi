""" Read output files of SRIM simulation

TODO: Read header information
"""
import os
import re
from io import BytesIO

import numpy as np

from .core.ion import Ion
from .core.element import Element

from typing import Union, NoReturn, List, Dict, Optional, Tuple
from mytypes import floatArray
from pathlib import Path

double_regex: str = ...
symbol_regex: str = ...
int_regex: str = ...


class SRIMOutputParseError(Exception):
    ...


class SRIM_Output(object):
    def _read_name(self, output: BytesIO) -> NoReturn: ...

    def _read_ion(self, output: bytes) -> Ion: ...

    def _read_target(self, output: bytes) -> NoReturn: ...

    def _read_num_ions(self, output: bytes) -> int: ...

    def _read_table(self, output: bytes) -> floatArray: ...


class Results(object):

    def __init__(self, directory: Union[str, Path]) -> None:
        self.ioniz: Ioniz = ...
        self.vacancy: Vacancy = ...
        self.novac: Optional[NoVacancy] = ...
        self.etorecoils: EnergyToRecoils = ...
        self.phonons: Phonons = ...
        self.range: Range = ...


class Ioniz(SRIM_Output):

    def __init__(self, directory: str, filename:  Union[str, Path] = ...) -> None: ...

    @property
    def ion(self) -> Ion: ...

    @property
    def num_ions(self) -> int: ...

    @property
    def depth(self) -> 'floatArray': ...
    #  return self._depth

    @property
    def ions(self): ...
    # return self._ions

    @property
    def recoils(self): ...
    # return self._recoils


class Vacancy(SRIM_Output):

    def __init__(self, directory: str, filename: str = 'VACANCY.txt') -> None: ...

    @property
    def ion(self) -> Ion: ...

    @property
    def num_ions(self) -> int: ...

    @property
    def depth(self) -> 'floatArray': ...
    # return self._depth

    @property
    def knock_ons(self): ...
    # return self._ion_knock_ons

    @property
    def vacancies(self): ...
    # return self._vacancies


class NoVacancy(SRIM_Output):

    def __init__(self, directory: str, filename: str = ...) -> None: ...

    @property
    def ion(self) -> Ion: ...

    @property
    def num_ions(self) -> int: ...

    @property
    def depth(self) -> 'floatArray': ...

    @property
    def number(self): ...


class EnergyToRecoils(SRIM_Output):

    def __init__(self, directory: str, filename: str = ...) -> None: ...

    @property
    def ion(self) -> Ion: ...

    @property
    def num_ions(self) -> int: ...

    @property
    def depth(self) -> 'floatArray': ...

    @property
    def ions(self): ...

    @property
    def absorbed(self): ...


class Phonons(SRIM_Output):

    def __init__(self, directory: str, filename: str = ...) -> None: ...

    @property
    def ion(self) -> Ion: ...

    @property
    def num_ions(self) -> int: ...

    @property
    def depth(self) -> 'floatArray': ...

    @property
    def ions(self): ...

    @property
    def recoils(self): ...


class Range(SRIM_Output):

    def __init__(self, directory: str, filename: str = ...) -> None: ...

    @property
    def ion(self) -> Ion: ...

    @property
    def num_ions(self) -> int: ...

    @property
    def depth(self) -> 'floatArray': ...

    @property
    def ions(self): ...

    @property
    def elements(self) -> Element: ...


class Backscat(object):
    ...


class Transmit(object):
    ...


class Sputter(object):
    ...


class Collision:

    def __init__(self, directory: str, filename: str = ...) -> None: ...

    def _read_header(self, f) -> List[str]: ...

    def _read_ion(self, ion_str: str) -> Dict[str, Union[float, int]]: ...

    def _read_cascade(self, lines) -> Tuple[Optional[float], Optional[float], Optional[float],
                                            Optional[float], List[Dict[str, Union[int, float]]]]: ...

    def __getitem__(self, i: int) -> Ion: ...

    def __len__(self) -> int: ...


def buffered_findall(filename: str, string: str, start: int = ...) -> List: ...


class SRResults(object):
    """Read SR_OUTPUT.txt file generated by pysrim SR.run()"""

    def __init__(self, directory: str, filename: str = ...) -> None: ...

    def _read_stopping_units(self, output: bytes) -> str: ...

    def _read_ion_info(self, output: bytes) -> Dict[str, Union[str, int, float]]:
        '''Example line to read from the file:
        Ion = Nickel       [28] , Mass = 58.6934 amu'''
        projectile_rexep = r'Ion\s+=\s+(.*?)\s+\[({})\]\s+, Mass\s+=\s({})\s+amu+\r\n'.format(int_regex, double_regex)
        match = re.findall(projectile_rexep.encode('utf-8'), output, re.DOTALL)
        out_dict = {
            'name': match[0][0].decode('utf-8'),
            'Z1': int(match[0][1]),
            'A1': float(match[0][2])
        }
        return out_dict

    def _read_target_info(self, output: bytes):
        '''lines to find from the file:
        Density =  2.3210E+00 g/cm3 = 4.9766E+22 atoms/cm3
        ======= Target  Composition ========
           Atom   Atom   Atomic    Mass
           Name   Numb   Percent   Percent
           ----   ----   -------   -------
            Si     14    100.00    100.00
        ====================================
        '''

        # first read the density info from the file
        density_reexp = r'Density\s+=\s+({})\s+g/cm3\s+=\s({})\s+atoms/cm3'.format(double_regex, double_regex)

        density_match = re.search(density_reexp.encode('utf-8'), output)

        density = np.array([density_match.group(1), density_match.group(2)], dtype='float')

        # find the target composition table
        # .format(symbol_regex, int_regex, double_regex, double_regex)#(=*)\r\n'
        table_regexp = r'=*\s+Target\s+Composition\s+=*\r\n(.*\r\n){3}((?:\s*.+\s\r\n)+)\s=*\r\n\s+Bragg Correction'
        table_match = re.search(table_regexp.encode('utf-8'), output)

        # rearrange the match into list of layer elements
        target_comp = table_match.groups()[-1].decode('utf-8').strip().split('\r\n')

        # create a dict object for target layers
        elements_dict = {}

        for line in target_comp:
            element = line.strip().split()
            Z = int(element[1])
            stoich_percent = float(element[2])
            mass_percent = float(element[3])
            elements_dict[element[0]] = [Z, stoich_percent, mass_percent]
            # print()

        # create a output dict
        target_dict = {'density g/cm3': density[0],
                       'density atoms/cm3': density[1],
                       'target composition': elements_dict
                       }

        return target_dict

    def _read_stopping_table(self, output: bytes) -> np.ndarray[List[float]]: ...

    @property
    def units(self) -> str: ...
    # return self._units

    @property
    def data(self) -> List[float]: ...

    @property
    def ion(self) -> Ion: ...

    @property
    def target(self) -> Dict[str, Union[float, Dict]]: ...
    # return self._target
