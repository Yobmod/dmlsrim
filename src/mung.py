
"""."""
from __future__ import annotations
import json
import pickle
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import os
from srim import Ion, Layer, Target  # , output
from srim.srim import TRIM
from srim.output import Results

from concurrent.futures import as_completed, ProcessPoolExecutor
import multiprocessing as mp
from time import sleep
from dataclasses import asdict  # , dataclass as dc
from pydantic.dataclasses import dataclass

from typing import cast, Iterable, Sequence, Set, Union, List, Tuple, Dict, NamedTuple
from typing_extensions import Literal, TypedDict

from .mytypes import floatArray, precisionLitType
from .dcf import SrimData
