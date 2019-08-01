import pytest
import os
from pathlib import Path

# from ../src/main import
from ..src.main import ion_Ni, layer_Ni, make_image_path


def test_make_image_path() -> None:
    out_image_path = make_image_path(layer, ion)
    assert os.path.exists(out_image_path) is True

    image_path = Path(r'.\images')
    out_image_path = make_image_path(layer, ion, image_path)
    assert os.path.exists(out_image_path) is True

    image_path = r'.\images'
    out_image_path = make_image_path(layer, ion, image_path)
    assert os.path.exists(out_image_path) is True
