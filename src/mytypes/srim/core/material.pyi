import re
from typing import Dict, List, Mapping, Union

from .utils import (
    check_input,
    is_positive, is_greater_than_zero,
    is_zero_or_one
)
from .element import Element
from .types import elemParamsType


class Material(object):
    """ Material Representation """

    def __init__(
        self, elements: Union[Mapping[Element, elemParamsType], Mapping[str, elemParamsType]],
            density: float, phase: int = 0) -> None:
        """Create Material from elements, density, and phase

        Parameters
        ----------
        elements : :obj:`dict`
             dictionary of elements (:class:`srim.core.elements.Element`, :obj:`str`, or :obj:`int`) with properties
               - ``stoich``  (float, int, required): Stoichiometry of element (fraction)
               - ``E_d``     (float, int, optional): Displacement energy [eV] default 25.0 eV
               - ``lattice`` (float, int, optional): Lattice binding energies [eV] default 0.0 eV
               - ``surface`` (float, int, optional): Surface binding energies [eV] default 3.0 eV
        density : :obj:`float`
             density [g/cm^3] of material
        phase : :obj:`int`
             phase of material (solid = 0, gas = 1). Default solid (0).


        Notes
        -----
        This class is more featureful that `srim.core.layer.Layer`
        would lead you to believe. In general this class will not be
        called by the user.

        Structure of dictionary elements properties:
         - stoich  (required): Stoichiometry of element (fraction)
         - E_d     (optional): Displacement energy [eV] default 25.0 eV
         - lattice (optional): Lattice binding energies [eV] default 0.0 eV
         - surface (optional): Surface binding energies [eV] default 3.0 eV

        dictionary element properties can be:

        float or int: stoich
          all others take default values for now

        dictionary:
          {'stoich', 'E_d', 'lattice', 'surface'}
          stoich is required all others are optional

        elements list structure:
          [stoich, E_d, lattice, surface]
          first element is required all others optional

        For example a single element in elements can be specified as:
          - {'Cu': 1.0}
          - {Element('Cu'): 1.0}
          - {Element('Cu'): [1.0, 25.0]}
          - {'Cu': {'stoich': 1.0}}
          - {Element('Cu'): {'stoich': 1.0, 'E_d': 25.0, 'lattice': 0.0, 'surface': 3.0}

        All stoichiometries will be normalized to 1.0

        Eventually the materials will have better defaults that come
        from databases.
        """
        self.phase = phase
        self.density = density
        self.elements: Dict[Element, Dict[str, float]] = {}

        stoich_sum = 0.0
        for element, values in elements.items():

            if not isinstance(element, Element):
                element = Element(element)

            if isinstance(values, dict):
                stoich = values['stoich']
                e_disp = values.get('E_d', 25.0)
                lattice = values.get('lattice', 0.0)
                surface = values.get('surface', 3.0)
            elif isinstance(values, list):
                default_values = [0.0, 25.0, 0.0, 3.0]
                if len(values) == 0 or len(values) > 4:
                    raise ValueError('list must be 0 < length < 5')
                values = values + default_values[len(values):]
                stoich, e_disp, lattice, surface = values
            elif isinstance(values, (int, float)):
                stoich = values
                e_disp = 25.0
                lattice = 0.0
                surface = 3.0
            else:
                raise ValueError('elements must be of type int, float, list, or dict')

            # Check input
            stoich = check_input(float, is_greater_than_zero, stoich)
            e_disp = check_input(float, is_positive, e_disp)
            lattice = check_input(float, is_positive, lattice)
            surface = check_input(float, is_positive, surface)

            stoich_sum += stoich

            self.elements.update({element: {
                'stoich': stoich, 'E_d': e_disp,
                'lattice': lattice, 'surface': surface
            }})

        # Normalize the Chemical Composisiton to 1.0
        for validated_element in self.elements.values():
            validated_element['stoich'] /= stoich_sum

    @ classmethod
    def from_formula(cls, chemical_formula: str, density: float, phase: int = 0) -> 'Material':
        """ Creation Material from chemical formula string and density

        Parameters
        ----------
        chemical_formula : :obj:`str`
            chemical formula string in specific format
        density : :obj:`float`
            density [g/cm^3] of material
        phase : :obj:`int`, optional
            phase of material (solid = 0, gas = 1). Default solid (0).

        Notes
        -----
        Examples of chemical_formula that can be used:
         - SiC
         - CO2
         - AuFe1.5
         - Al10.0Fe90.0

        Chemical Formula will be normalized to 1.0
        """
        elements = cls._formula_to_elements(chemical_formula)
        return Material(elements, density, phase)

    @ staticmethod
    def _formula_to_elements(chemical_formula: str) -> Dict[Element, Dict[str, float]]:
        """ Convert chemical formula to elements """
        single_element = r'([A-Z][a-z]?)([0-9]*(?:\.[0-9]*)?)?'
        elements: Dict[Element, Dict[str, float]] = {}

        if re.match('^(?:{})+$'.format(single_element), chemical_formula):
            matches = re.findall(single_element, chemical_formula)
        else:
            error_str = 'chemical formula string {} does not match regex'
            raise ValueError(error_str.format(chemical_formula))

        # Check for errors in stoichiometry
        for symbol, fraction in matches:
            element = Element(symbol)

            if element in elements:
                error_str = 'cannot have duplicate elements {} in stoichiometry'
                raise ValueError(error_str.format(element.symbol))

            if fraction == '':
                fraction = 1.0

            elements.update({element: {'stoich': float(fraction)}})
        return elements

    @ property
    def density(self) -> float:
        """Material's density"""
        return self._density

    @ density.setter
    def density(self, value: float) -> None:
        self._density: float = check_input(float, is_positive, value)

    @ property
    def phase(self) -> int:
        """Material's phase"""
        return self._phase

    @ phase.setter
    def phase(self, value: int) -> None:
        self._phase: int = check_input(int, is_zero_or_one, value)

    @ property
    def chemical_formula(self) -> str:
        """Material's chemical formula"""
        return " ".join(f"{elem.symbol} {params['stoich']:1.2f}" for elem, params in self.elements.items())

    def __repr__(self) -> str:
        return f"<Material formula:{self.chemical_formula} density:{self.density:2.3f}>"

    def __eq__(self, material: object) -> bool:

        if not isinstance(material, Material):
            return False

        elif abs(self.density - material.density) > 1e-6:
            return False

        elif len(self.elements) != len(material.elements):
            return False

        for element in self.elements:
            if element not in material.elements:
                return False
            for prop in self.elements[element]:
                if abs(self.elements[element][prop] - material.elements[element][prop]) > 1e-6:
                    return False
        return True
