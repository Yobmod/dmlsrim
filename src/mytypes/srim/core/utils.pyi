"""Utility functions that are used to construct Target and Ion

"""


def check_input(input_type, condition, value):
    value = input_type(value)
    if not condition(value):
        raise ValueError('type of argument does not satisfy condition')
    return value


def is_zero(value): return True if value == 0 else False


def is_zero_or_one(value): return True if value in range(2) else False


def is_zero_to_two(value): return True if value in range(3) else False


def is_zero_to_five(value): return True if value in range(6) else False


def is_one_to_seven(value): return True if value in range(1, 8) else False
def is_one_to_eight(value): return True if value in range(1, 9) else False
def is_srim_degrees(value): return True if 0.0 <= value < 90.0 else False
def is_positive(value): return True if value >= 0.0 else False
def is_greater_than_zero(value): return True if value > 0.0 else False
def is_quoteless(value): return True if '"' not in value else False
