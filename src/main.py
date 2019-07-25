import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from srim import Ion, Layer, Target, TRIM
from srim.output import Results


from typing import Union
from typing_extensions import Literal


# Construct a 3MeV Nickel ion
ion = Ion('Ni', energy=3.0e6)
# print(ion.symbol, ion.name)
# print(ion.energy)


# Construct a layer of nick 20um thick with a displacement energy of 30 eV
layer = Layer(
    {'Ni': {'stoich': 1.0, 'E_d': 30.0, 'lattice': 0.0, 'surface': 3.0}},
    density=8.9, width=20_000.0
)
"""
print(layer.elements.keys())
print(layer.width)  # angstroms
for element, prop in layer.elements.items():
    print(element.symbol, element.name, element.mass)
    print(prop['stoich'], prop['E_d'], prop['lattice'], prop['surface'])
"""


def make_element_subfolder_name(layer: Layer, ion: Ion,
                                precision: Literal['um', 'nm', 'A', 'a', 'micro', 'nano', 'angstrom'] = 'um') -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy
    data_path default = '.\\data'. precision is units of the layer width, default = 'um' """

    for (element, prop) in layer.elements.items():
        stoich = prop['stoich']
        if stoich == 1.0:
            element_str = element.symbol
        elif stoich.is_integer():
            element_str = f'element.symbol{stoich:.0f}'
        else:
            element_str = f'element.symbol{stoich:.2f}'
        element_list_str = "".join(element_str)

    layer_width_nm = f'{layer.width / 10:.0f}nm'
    layer_width_um = f'{layer.width / 10000:.0f}um'
    ion_energy_kev = f'{ion.energy / 1000:.0f}keV'

    if precision == 'um':
        layer_width = layer_width_um
    elif precision == 'nm':
        layer_width = layer_width_nm
    else:
        layer_width = layer.width

    data_subfolder_name = Path(f"{element_list_str}_{layer_width}_{ion.symbol}@{ion_energy_kev}")
    return data_subfolder_name


def make_data_path(layer: Layer, ion: Ion,
                   data_path: Union[Path, str] = Path(R'.\data'),
                   precision: Literal['um', 'nm', 'A', 'a', 'micro', 'nano', 'angstrom'] = 'um') -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy
    data_path default = '.\\data'. precision is units of the layer width, default = 'um' """

    data_subfolder_name = make_element_subfolder_name(layer, ion, precision)
    output_directory: Path = Path(data_path) / data_subfolder_name
    output_directory.mkdir(parents=True, exist_ok=True)
    return output_directory


def make_image_path(layer: Layer, ion: Ion,
                    image_path: Union[Path, str] = Path(R'.\images'),
                    precision: Literal['um', 'nm', 'A', 'a', 'micro', 'nano', 'angstrom'] = 'um') -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy
    data_path default = '.\\data'. precision is units of the layer width, default = 'um' """

    data_subfolder_name = make_element_subfolder_name(layer, ion, precision)
    outimage_directory: Path = Path(image_path) / data_subfolder_name
    outimage_directory.mkdir(parents=True, exist_ok=True)
    return outimage_directory


# Construct a target of a single layer of Nickel
target = Target([layer])


# Initialize a TRIM calculation with given target and ion for 25 ions, quick˓→calculation
trim = TRIM(target, ion, number_ions=25, calculation=1)

# Specify the directory of SRIM.exe# For windows users the path will include  C://...
srim_executable_directory = 'C://srim'

# takes about 10 seconds on my laptop
results = trim.run(srim_executable_directory)  # If all went successfull you should have seen a TRIM window popup and run 25 ions!


# os.makedirs(output_directory, exist_ok=True)
# top_level_csv_files = Path.cwd().glob('*.csv')
# all_csv_files = Path.cwd().rglob('*.csv')
# from shutil import copyfile
# copyfile(source, destination)
data_path = Path(R'.\data')
image_path = Path(R'.\images')
data_out_dir = make_data_path(layer, ion, data_path)

TRIM.copy_output_files(srim_executable_directory, data_out_dir)
# results = Results(output_directory)


# to units of Angstromsenergy_damage = (phon.ions + phon.recoils)*dx
def plot_damage_energy(folder: Path, ax: plt.axis) -> float:
    results = Results(folder)
    phon = results.phonons
    dx = max(phon.depth) / 100.0  # to units of Angstroms
    energy_damage = (phon.ions + phon.recoils) * dx
    ax.plot(phon.depth, energy_damage / phon.num_ions, label='{}'.format(folder))
    return sum(energy_damage)


folders = [data_out_dir]
image_out_dir = make_image_path(layer, ion)
os.makedirs(image_out_dir, exist_ok=True)

fig, axes = plt.subplots(1, len(folders), sharex=True, sharey=True)

for ax, folder in zip(np.ravel(axes), folders):
    energy_damage = plot_damage_energy(folder, ax)
    energy_damage_kev = energy_damage / 1000
    print("Damage energy: {:.1f} keV".format(energy_damage_kev))
    ax.set_xlabel('Depth [Angstroms]')
    ax.set_ylabel('eV')
    ax.legend()

fig.suptitle('Damage Energy vs. Depth', fontsize=15)
fig.set_size_inches((20, 6))
fig.savefig(os.path.join(image_out_dir, 'damagevsdepth.png'), transparent=True)
