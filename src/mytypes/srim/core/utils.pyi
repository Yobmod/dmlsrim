"""Utility functions that are used to construct Target and Ion

"""


from typing import Any, Callable, Type, TypeVar, Union, cast, overload

inp = Union[str, float, int]
_inp = TypeVar('_inp', bound=inp)


@overload
def check_input(
    input_converter: Callable[[_inp], str],
    condition: Callable[[str], bool],
    value: _inp,
) -> str: ...


@overload
def check_input(
    input_converter: Callable[[_inp], float],
    condition: Callable[[float], bool],
    value: _inp,
) -> float: ...


@overload
def check_input(
    input_converter: Callable[[_inp], int],
    condition: Callable[[int], bool],
    value: _inp,
) -> int: ...


def check_input(
    input_converter: Callable[[_inp], inp],
    condition: Callable[[inp], bool],
    value: _inp,
) -> _inp:
    conv_value = input_converter(value)
    if not condition(conv_value):
        raise ValueError('type of argument does not satisfy condition')
    else:
        return value


def is_zero(value: int) -> bool:
    return value == 0


def is_zero_or_one(value: int) -> bool:
    return value in range(2)


def is_zero_to_two(value: int) -> bool:
    return value in range(3)


def is_zero_to_five(value: int) -> bool:
    return value in range(6)


def is_one_to_seven(value: int) -> bool:
    return value in range(1, 8)


def is_one_to_eight(value: int) -> bool:
    return value in range(1, 9)


def is_srim_degrees(value: float) -> bool:
    return 0.0 <= value < 90.0


def is_positive(value: float) -> bool:
    return value >= 0.0


def is_greater_than_zero(value: float) -> bool:
    return value > 0.0


def is_quoteless(value: str) -> bool:
    return '"' not in value
