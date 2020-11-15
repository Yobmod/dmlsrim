from __future__ import annotations
import yaml
import os
import re
import srim
from typing import Dict, Any, Union, TYPE_CHECKING

from .types import element_TD

if TYPE_CHECKING:  # to prevent import cycle
    from .element import Element


def create_elementdb() -> Dict[str, Any]:
    dbpath = os.path.join(srim.__path__[0], 'data', 'elements.yaml')
    db_dict: Dict[str, element_TD] = yaml.load(open(dbpath, "r"), Loader=yaml.FullLoader)
    return db_dict


class ElementDB(object):
    """Element database at ``srim.data.elements.yaml``"""
    _db = create_elementdb()

    @classmethod
    def lookup(cls, identifier: Union[str, int]) -> element_TD:
        """ Looks up element from symbol, name, or atomic number

        Parameters
        ----------
        identifier : :obj:`str`, :obj:`int`
            Unique symbol, name, or atomic number of element

        Notes
        -----
            This class is used for creation of elements, ions,
            etc. but generally will not be needed by the user.
        """
        if isinstance(identifier, (bytes, str)):
            if re.match("^[A-Z][a-z]?$", identifier):   # Symbol
                return cls._lookup_symbol(identifier)
            elif re.match("^[A-Z][a-z]*$", identifier):  # Name
                return cls._lookup_name(identifier)
        elif isinstance(identifier, int):               # Atomic Number
            return cls._lookup_atomic_number(identifier)
        raise ValueError('identifier of type:{} value:{} not value see doc'.format(
            type(identifier), identifier))

    @classmethod
    def _lookup_symbol(cls, symbol: str) -> element_TD:
        """ Looks up symbol in element database

        :param str symbol: Symbol of atomic element
        """
        db_symbol: element_TD = cls._db[symbol]
        return db_symbol

    @classmethod
    def _lookup_name(cls, name: str) -> element_TD:
        """ Looks element in database by name

        :param str name: (Full) Name of atomic element (British spelling)
        """
        for symbol in cls._db:
            if cls._db[symbol]['name'] == name:
                db_symbol: element_TD = cls._db[symbol]
                return db_symbol
        raise KeyError('name:{} does not exist'.format(name))

    @classmethod
    def _lookup_atomic_number(cls, atomic_number: int) -> element_TD:
        """ Look up element in database by atomic number (Z)

        :param int atomic_number: Atomic number of atomic element
        """
        for symbol in cls._db:
            if cls._db[symbol]['z'] == atomic_number:
                db_symbol: element_TD = cls._db[symbol]
                return db_symbol
        raise IndexError(f'atomic number:{atomic_number} does not exist')
