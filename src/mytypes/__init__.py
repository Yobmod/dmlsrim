from typing import TYPE_CHECKING, Iterable, Tuple, Union, List, NewType  # , Any, NewType, TypeVar
from typing_extensions import Literal
from PIL import Image
import numpy as np

precisionLitType = Literal['um', 'nm', 'A', 'a', 'micro', 'nano', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms']


if TYPE_CHECKING:  # for mypy
    imageType = NewType('imageType', np._ArrayLike[int])  # requres numpy stubs to be discoverable by mypy
    intArray = np.ndarray[int]
    floatArray = np.ndarray[float]
else:  # python doesnt allow ndarray subscriptable yet
    imageType = np.ndarray
    intArray = np.ndarray
    floatArray = np.ndarray


class PILImage(Image.Image):
    ...


colorType = Union[Tuple[int, int, int, int], Tuple[int, int, int]]
pointType = Tuple[int, int]
