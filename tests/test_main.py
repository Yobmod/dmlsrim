
import pytest
import os
import sys
from pathlib import Path
from colorama import init as color_init, Fore, Back  # , Style

my_path = Path(__file__).parent.parent  # print(my_path)
sys.path.insert(0, str(my_path))
sys.path.insert(0, str(my_path) + R"\src")

from ..src.main import ion_Ni, layer_Ni, make_image_path, make_data_path, srim_executable_directory  # NOQA

color_init(autoreset=True)
AssertionErrorColor = Fore.GREEN + Back.RED


@pytest.mark.paths
def test_srim_directory_path() -> None:
    """Test if srim directory used is found / exists"""

    assert srim_executable_directory.exists() is True, AssertionErrorColor + "Srim.exe path not found. Is it installed?"


@pytest.mark.paths
def test_make_data_path() -> None:
    """Test creation of data paths using default args and given paths"""

    out_data_path = make_data_path(layer_Ni, ion_Ni)
    assert os.path.exists(out_data_path) is True, "image path with default path not found/created"

    data_path = Path(R'.\data')
    out_data_path = make_data_path(layer_Ni, ion_Ni, data_path)
    assert os.path.exists(out_data_path) is True, R"data path with '.\data' arg not found/created"


@pytest.mark.paths
def test_make_image_path() -> None:
    """Test creation of image paths using default args and given paths"""

    out_image_path = make_image_path(layer_Ni, ion_Ni)
    assert os.path.exists(out_image_path) is True, "image path with default path not found/created"

    image_path = Path(R'.\images')
    out_image_path = make_image_path(layer_Ni, ion_Ni, image_path)
    assert out_image_path.exists() is True, R"image path with '.\image' arg not found/created"
