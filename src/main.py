import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from srim import Ion, Layer, Target, TRIM
from srim.output import Results

from concurrent.futures import ThreadPoolExecutor

from typing import Union, List
from typing_extensions import Literal, Final


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

precisionLitType = Literal['um', 'nm', 'A', 'a', 'micro', 'nano', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms']


def make_element_subfolder_name(layer: Layer, ion: Ion,
                                precision: precisionLitType = 'um') -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy.
    precision is units of the layer width, default = 'um' """

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

    if precision == 'um' or precision == 'micro':
        layer_width = layer_width_um
    elif precision == 'nm' or precision == 'nano':
        layer_width = layer_width_nm
    else:
        layer_width = layer.width

    data_subfolder_name = Path(f"{element_list_str}_{layer_width}_{ion.symbol}@{ion_energy_kev}")
    return data_subfolder_name


def make_data_path(layer: Layer,
                   ion: Ion,
                   data_path: Union[Path, str] = Path(R'.\data'),
                   precision: precisionLitType = 'um') -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy
    data_path default = '.\\data'. precision is units of the layer width, default = 'um' """

    data_subfolder_name = make_element_subfolder_name(layer, ion, precision)
    output_directory: Path = Path(data_path) / data_subfolder_name
    output_directory.mkdir(parents=True, exist_ok=True)
    return output_directory


def make_image_path(layer: Layer, ion: Ion,
                    image_path: Union[Path, str] = Path(R'.\images'),
                    precision: precisionLitType = 'um') -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy
    data_path default = '.\\images'. precision is units of the layer width, default = 'um' """

    data_subfolder_name = make_element_subfolder_name(layer, ion, precision)
    outimage_directory: Path = Path(image_path) / data_subfolder_name
    outimage_directory.mkdir(parents=True, exist_ok=True)
    return outimage_directory


# Construct a target of a single layer of Nickel
target = Target([layer])


# Initialize a TRIM calculation with given target and ion for 25 ions, quick˓→calculation
trim = TRIM(target, ion, number_ions=25, calculation=1)

# Specify the directory of SRIM.exe# For windows users the path will include  C://...
srim_executable_directory: Final = R'C:\srim'

# takes about 10 seconds on my laptop
results = trim.run(srim_executable_directory)  # If all went successfull you should have seen a TRIM window popup and run 25 ions!


# os.makedirs(output_directory, exist_ok=True)
# top_level_csv_files = Path.cwd().glob('*.csv')
# all_csv_files = Path.cwd().rglob('*.csv')
# from shutil import copyfile
# copyfile(source, destination)
data_path: Final = Path(R'.\data')
image_path: Final = Path(R'.\images')
data_out_dir = make_data_path(layer, ion, data_path)

TRIM.copy_output_files(srim_executable_directory, data_out_dir)
# results = Results(output_directory)


def get_energy_damage_array(results: Results) -> np.ndarray:
    phon = results.phonons
    dx = max(phon.depth) / 1000.0  # units from pm to nm
    energy_damage = (phon.ions + phon.recoils) * dx
    damage_array_nm = np.array(phon.depth / 1000, energy_damage / phon.num_ions)
    return damage_array_nm


def calc_energy_damage(results: Results, units: precisionLitType = 'nm') -> List[float]:
    phon = results.phonons
    if units in ('nm', 'nano'):
        dx = max(phon.depth) / 1000.0  # units from pm to nm
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        dx = max(phon.depth) / 100.0  # units from pm to Angstroms
    energy_damage: List[float] = (phon.ions + phon.recoils) * dx
    return energy_damage


def plot_damage_energy(results: Results, ax: plt.axis, units: precisionLitType = 'nm') -> None:
    phon = results.phonons
    if units in ('nm', 'nano'):
        units_str = 'nm'
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        units_str = 'Angstroms'
    energy_damage = calc_energy_damage(results, units)
    ax.plot(phon.depth, energy_damage / phon.num_ions, label='{}'.format(folder))
    ax.set_xlabel(f'Depth [{units_str}]')
    ax.set_ylabel('eV')
    ax.legend()





if __name__ == "__main__":
    folders: List[Union[str, Path]] = [data_out_dir]
    image_out_dir = make_image_path(layer, ion)
    os.makedirs(image_out_dir, exist_ok=True)

    fig, axes = plt.subplots(1, len(folders), sharex=True, sharey=True)
    for ax, folder in zip(np.ravel(axes), folders):
        results = Results(folder)
        energy_damage_sum: float = sum(calc_energy_damage(results))
        energy_damage_kev = energy_damage_sum / 1000
        print("Damage energy: {:.1f} keV".format(energy_damage_kev))
        plot_damage_energy(results, ax, units='nm')

    with ThreadPoolExecutor() as exc:
        pass  # add folders to list, get results list?, then map jobs to folders?

    fig.suptitle('Damage Energy vs. Depth', fontsize=15)
    fig.set_size_inches((20, 6))
    fig.savefig(os.path.join(image_out_dir, 'damagevsdepth.png'), transparent=True)
