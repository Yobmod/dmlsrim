import os
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import matplotlib.pyplot as plt
from colorama import init as color_init, Fore

from srim import Ion, Layer, Target  # , TRIM
from srim.srim import TRIM
from srim.output import Results
import srim
import yaml
from srim.core import elementdb, element

import typing as t
from typing import Union, List, cast
from mytypes import floatArray, precisionLitType

if sys.version_info < (3, 8):
    from typing_extensions import Final
else:
    from typing import Final

color_init(autoreset=True)

dbpath = os.path.join(srim.__path__[0], 'data', 'elements.yaml')  # type: ignore
yaml.load(open(dbpath, "r"), Loader=yaml.FullLoader)

# this via instance?
element_database = elementdb.ElementDB()
sulphur: element.Element = element_database.lookup('S')
print(sulphur)
# orthis bacause classmethod?
sulphur = elementdb.ElementDB.lookup('S')
print(sulphur)

srim_executable_directory: Final = Path(R'C:\srim')


# TODO get a .csv of element data, import and extract line as dict?
# TODO element picker GUI
elem_ce_dict = {'E_d': 25.0, 'lattice': 3.0, 'surface': 4.23, 'atomic_num': 58, 'atomic_mass': 140.1}
elem_u_dict = {'E_d': 25.0, 'lattice': 3.0, 'surface': 5.42, 'atomic_num': 92, 'atomic_mass': 238.0}
elem_th_dict = {'E_d': 25.0, 'lattice': 3.0, 'surface': 5.93, 'atomic_num': 90, 'atomic_mass': 232.0}
elem_o_dict = {'E_d': 28.0, 'lattice': 3.0, 'surface': 2.00, 'atomic_num': 8, 'atomic_mass': 15.99}

cmpd_ceo2_dict = {
    'Ce': {'stoich': 1.0, **elem_ce_dict},
    'O': {'stoich': 2.0, **elem_o_dict},
}

energy__kev_list = [100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000]


ion_Ni = Ion('Ni', energy=3.0e6)
ions_He_list = [Ion('He', energy=x * 1000) for x in energy__kev_list]
# print(ion.symbol, ion.name)
# print(ion.energy)

# Construct layers
layer_Ni = Layer(
    {'Ni': {'stoich': 1.0, 'E_d': 30.0, 'lattice': 0.0, 'surface': 3.0}},
    density=8.9, width=20_000.0, name='Ni'
)


"""
print(layer.elements.keys())
print(layer.width)  # angstroms
for element, prop in layer.elements.items():
    print(element.symbol, element.name, element.mass)
    print(prop['stoich'], prop['E_d'], prop['lattice'], prop['surface'])
"""


# TODO - use parser or chempy for this?
def make_element_subfolder_name(layer: Layer, ion: Ion,
                                precision: 'precisionLitType' = 'um') -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy.
    precision is units of the layer width, default = 'um' """

    if layer.name:
        element_list_str = layer.name
    else:
        element_list = []
        for (elem, prop) in layer.elements.items():
            stoich = prop['stoich']
            # print(element.symbol, stoich)
            if stoich == 1.0:
                element_str = elem.symbol
            elif stoich.is_integer():
                element_str = f'{elem.symbol}{stoich:.0f}'
            else:
                element_str = f'{elem.symbol}{stoich:.2f}'
            # print(element_str)
            element_list.append(element_str)
        element_list_str = "-".join(element_list)
        # print(element_list_str)

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
    # print(data_subfolder_name)
    return data_subfolder_name


def make_data_path(layer: Layer,
                   ion: Ion,
                   data_path: Union[Path, str] = Path(R'..\data'),
                   precision: precisionLitType = 'um') -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy
    data_path default = '..\\data'. precision is units of the layer width, default = 'um' """

    data_subfolder_name = make_element_subfolder_name(layer, ion, precision)
    output_directory: Path = Path(data_path) / data_subfolder_name
    output_directory.mkdir(parents=True, exist_ok=True)
    return output_directory


def make_image_path(layer: Layer, ion: Ion,
                    image_path: Union[Path, str] = Path(R'..\images'),
                    precision: precisionLitType = 'um') -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy
    data_path default = '..\\images'. precision is units of the layer width, default = 'um' """

    data_subfolder_name = make_element_subfolder_name(layer, ion, precision)
    outimage_directory: Path = Path(image_path) / data_subfolder_name
    outimage_directory.mkdir(parents=True, exist_ok=True)
    return outimage_directory


def get_energy_damage_array(results: Results) -> 'floatArray':
    phon = results.phonons
    dx = max(phon.depth) / 1000.0  # units from pm to nm
    energy_damage = (phon.ions + phon.recoils) * dx
    damage_array_nm = np.array(phon.depth / 1000, energy_damage / phon.num_ions)
    return damage_array_nm


def calc_energy_damage(results: Results, units: precisionLitType = 'nm') -> 'floatArray':
    phon = results.phonons
    if units in ('nm', 'nano'):
        dx = max(phon.depth) / 100.0  # units from pm to nm
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        dx = max(phon.depth) / 100.0  # units from pm to Angstroms
    energy_damage: floatArray = (phon.ions + phon.recoils) * dx  # add the arrays and multiply
    cast(floatArray, energy_damage)
    # reveal_type(energy_damage)
    print(energy_damage)
    return energy_damage


def save_damage_energy(results: Results, units: precisionLitType = 'nm') -> None:
    phon = results.phonons
    if units in ('nm', 'nano'):
        units_str: precisionLitType = 'nm'
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        units_str = 'A'
    energy_damage = calc_energy_damage(results, units_str)
    damage_table: floatArray = np.array([[phon.depth], [energy_damage / phon.num_ions]])
    print(damage_table)
    print(f"{damage_table.shape}")
    # TODO save to csv, use units_str for header


def plot_damage_energy(results: Results, ax: plt.Axes, units: precisionLitType = 'nm', plot_label: str = 'Collision damage depth') -> None:
    phon = results.phonons
    if units in ('nm', 'nano'):
        units_str = 'nm'
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        units_str = 'Angstroms'
    energy_damage = calc_energy_damage(results, units)
    ax.plot(phon.depth, energy_damage / phon.num_ions, label=f'{plot_label}')
    ax.set_xlabel(f'Depth [{units_str}]')
    ax.set_ylabel('eV')
    ax.legend()


if __name__ == "__main__":

    # Construct a target of a single layer of Nickel
    # Initialize a TRIM calculation with given target and ion for 25 ions, quick˓→calculation
    # Specify the directory of SRIM.exe# For windows users the path will include  C://...
    # takes about 10 seconds on my laptop
    # If all went successfull you should have seen a TRIM window popup and run 25 ions!
    target_Ni = Target([layer_Ni])
    trim_Ni = TRIM(target_Ni, ion_Ni, number_ions=25, calculation=1)
    results = trim_Ni.run(srim_executable_directory)

    # os.makedirs(output_directory, exist_ok=True)
    # top_level_csv_files = Path.cwd().glob('*.csv')
    # all_csv_files = Path.cwd().rglob('*.csv')
    # from shutil import copyfile
    # copyfile(source, destination)

    # better to get cwd.parent then add data or images subfolder. Then src and tests are both using same folders. Also have to fix for tests?
    data_path: Final = Path(R'..\data')
    image_path: Final = Path(R'..\images')
    data_out_dir = make_data_path(layer_Ni, ion_Ni, data_path)
    TRIM.copy_output_files(srim_executable_directory, data_out_dir)
    # results = Results(output_directory)

    width_list = [10_000, 20_000, 50_000]

    layer_ceria_list = [Layer(cmpd_ceo2_dict, density=7.22, width=width, name='ceria') for width in width_list]

    layer_list = layer_ceria_list + [layer_Ni]

    target_list = [Target([layer]) for layer in layer_list]

    results_list = []
    for ion in ions_He_list:
        for i, target in enumerate(target_list):
            data_out_dir = make_data_path(layer_list[i], ion, data_path)
            trim = TRIM(target, ion, number_ions=25, calculation=1)
            results = trim.run(srim_executable_directory)
            results_list.append(results)
            TRIM.copy_output_files(srim_executable_directory, data_out_dir)
            print(f'{ion.symbol}-{ion.energy/1000}kev done')

    """to use threading, need to generate different srim data dir for each thread? Worth it?
    or just threading for result analysis?"""

    folders: List[Union[str, Path]] = [data_out_dir]
    image_out_dir = make_image_path(layer_Ni, ion_Ni)
    os.makedirs(image_out_dir, exist_ok=True)

    fig: plt.Figure
    axes: plt.Axes
    fig, axes = plt.subplots(1, len(folders), sharex=True, sharey=True)

    for ax, folder in zip(np.ravel(axes), folders):
        srim_res = Results(folder)
        energy_damage_sum: float = sum(cast(t.Iterable[float], calc_energy_damage(srim_res)))
        energy_damage_kev = energy_damage_sum / 1000
        print(Fore.GREEN + f"Damage energy: {energy_damage_kev:.1f} keV")
        plot_damage_energy(srim_res, ax, units='nm')

    res_folder_list: List[Union[str, Path]] = []
    ax_list: List[plt.Axes] = []
    with ThreadPoolExecutor() as exc:
        srim_results_list = exc.map(Results, res_folder_list)
        exc.map(plot_damage_energy, srim_results_list, ax_list)  # , units='nm')
        # add folders to list, get results list?, then map jobs to folders?

    fig.suptitle('Damage Energy vs. Depth', fontsize=15)
    fig.set_size_inches((20, 6))
    fig.savefig(os.path.join(image_out_dir, 'damagevsdepth.png'), transparent=True)
