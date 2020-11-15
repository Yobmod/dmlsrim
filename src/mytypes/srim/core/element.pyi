from typing import Union
from .elementdb import ElementDB


class Element(object):
    """ Element from periodic table

    Parameters
    ----------
    identifier : :obj:`str`, :obj:`int`
        Symbol, Name, or Atomic Number of element
    mass : :obj:`float`, optional
        Mass [amu] of element. Default is most common isotope atomic
        weight

    Examples
    --------
    Constructing a Helium Atom.

    >>> Element('He')
    <Element symbol:He name:Helium mass:4.00>

    >>> Element('Helium')
    <Element symbol:He name:Helium mass:4.00>

    >>> Element(2)
    <Element symbol:He name:Helium mass:4.00>

    >>> Element('He', 4.3)
    <Element symbol:He name:Helium mass:4.30>
    """

    def __init__(self, identifier: Union[str, int], mass: float = 0.0):
        """Initializes element from identifier and mass
        Identifier: str | int (symbol, name or atomic_number)"""
        element_dict = ElementDB.lookup(identifier)

        self._symbol = element_dict['symbol']
        self._name = element_dict['name']
        self._atomic_number = element_dict['z']
        self._mass = mass or element_dict['mass']

    def __eq__(self, other_element: object) -> bool:
        if not isinstance(other_element, Element):
            return False
        elif (self.symbol == other_element.symbol and
              self.name == other_element.name and
              self.atomic_number == other_element.atomic_number and
              self.mass == other_element.mass):
            return True
        else:
            return False

    def __repr__(self) -> str:
        return "<Element symbol:{} name:{} mass:{:2.2f}>".format(
            self.symbol, self.name, self.mass)

    def __hash__(self) -> int:
        return sum(hash(item) for item in [
            self._mass, self._symbol, self._name, self.atomic_number
        ])

    @property
    def symbol(self) -> str:
        """Element's atomic symbol"""
        assert isinstance(self._symbol, str)
        return self._symbol

    @property
    def name(self) -> str:
        """Element's formal name"""
        assert isinstance(self._name, str)
        return self._name

    @property
    def atomic_number(self) -> int:
        """Element's atomic number"""
        assert isinstance(self._atomic_number, int)
        return self._atomic_number

    @property
    def mass(self) -> float:
        """Element's mass"""
        assert isinstance(self._mass, float)
        return self._mass
