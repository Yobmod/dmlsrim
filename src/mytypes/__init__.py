from typing import Type, TYPE_CHECKING, Tuple, Union, List, NewType  # , Any, NewType, TypeVar
from PIL import Image
import numpy as np


if TYPE_CHECKING:  # for mypy
    class ndarray(np.ndarray): ...
    imageType = NewType('imageType', np._ArrayLike[np.ndarray[np.ndarray[int]]])
    intArray = np.ndarray[int]
    floatArray = np.ndarray[float]
else:  # python doesnt allow ndarray subscriptable yet
    imageType = np.ndarray
    intArray = np.ndarray
    floatArray = np.ndarray


class PILImage(Image.Image): ...


colorType = Union[Tuple[int, int, int, int], Tuple[int, int, int]]
pointType = Tuple[int, int]
