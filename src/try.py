from typing import NamedTuple
# from typing import Mapping, MutableMapping, Sequence, MutableSequence
from typing_extensions import TypedDict


class TD(TypedDict, total=True):
    x: int
    y: int
    label: str


"""
TD = TypedDict('TD', {
    'x': int,
    'y': int,
    'label': str,
})
"""


def test_TD(inp: TD) -> bool:
    if inp['label'] is True:
        return True
    else:
        return False


my_TD: TD = {'x': 1, 'y': 2, 'label': "str"}
res = test_TD(my_TD)


class NT(NamedTuple):
    x: int
    y: int
    label: str


"""
NT = NamedTuple('NT', [
    ('x', int),
    ('y', int),
    ('label', str),
])
"""


def test_NT(inp: NT) -> bool:
    if inp.label is True:
        return True
    else:
        return False


my_NT: NT = NT(1, 2, 'working')

res = test_NT(my_NT)
