from typing import Tuple, cast      # NOQA
from typing import Union
from typing import Dict
from typing import NamedTuple
# from typing import Mapping, MutableMapping, Sequence, MutableSequence
from typing_extensions import TypedDict, Final   # NOQA
from dataclasses import dataclass, asdict


class tDC(TypedDict):
    y: int
    label: int


f = (5, "str")
e: tDC = {'y': 5, 'label': 5}


@dataclass
class DC:
    y: int
    label: str

    def as_dict(self) -> Dict[str, Union[str, int]]:
        # narrow str, Any to declated dtypes
        return asdict(self)

    def as_typdict(self) -> tDC:
        ad: tDC = tDC(self.as_dict())
        return ad


dc = DC(*f)
my_dc = DC(2, "red")
dc_todict = my_dc.as_typdict()


class TD(TypedDict, total=False):
    x: int
    y: int
    label: str
    a: float


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

d = {'a': 3.2}

Tdd = d.update(**asdict(dc))


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


@dataclass
class Point:
    x: float
    y: float


@dataclass
class Line:
    a: Point
    b: Point


line = Line(Point(1, 2), Point(3, 4))


x = {1, 2, 3, 4, 5}

for a in x:
    print("a")
