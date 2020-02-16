
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from srim import Ion, Layer, Target, TRIM
from srim.output import Results

# from concurrent.futures import ThreadPoolExecutor

from typing import Union, List, Tuple
from typing_extensions import Final, TypedDict
from mytypes import floatArray, precisionLitType


srim_executable_directory: Final = R'C:\srim'


class Elem_Dict(TypedDict):
    E_d: float
    lattice: float
    surface: float
    atomic_num: int
    atomic_mass: float


elem_ce_dict: Elem_Dict = {'E_d': 25.0, 'lattice': 3.0, 'surface': 4.23, 'atomic_num': 58, 'atomic_mass': 140.1}
elem_u_dict: Elem_Dict = {'E_d': 25.0, 'lattice': 3.0, 'surface': 5.42, 'atomic_num': 92, 'atomic_mass': 238.0}
elem_th_dict: Elem_Dict = {'E_d': 25.0, 'lattice': 3.0, 'surface': 5.93, 'atomic_num': 90, 'atomic_mass': 232.0}
elem_o_dict: Elem_Dict = {'E_d': 28.0, 'lattice': 3.0, 'surface': 2.00, 'atomic_num': 8, 'atomic_mass': 15.99}
elem_si_dict: Elem_Dict = {'E_d': 28.0, 'lattice': 3.0, 'surface': 2.00, 'atomic_num': 7, 'atomic_mass': 15.99}
elem_ti_dict: Elem_Dict = {'E_d': 28.0, 'lattice': 3.0, 'surface': 2.00, 'atomic_num': 22, 'atomic_mass': 15.99}

{'HHI_P': 700.0, 'HHI_R': 1000.0, 'atomic_radii': 100.0, 'boil': 717.799987793,
 'color': '#FFFF30', 'covalent_radii': 102.0, 'd_elec': 0.0, 'density': 2.06699991226,
 'electronegativity': 2.57999992371, 'f_elec': 0.0, 'first_ionization_energy': 10.3599996567,
 'group': 16.0, 'mass': 32.0649986267, 'melt': 388.510009766, 'name': 'Sulfur',
 'p_elec': 4.0, 'period': 3.0, 'production': 350.0, 's_elec': 2.0,
 'scattering_factors':
 {'a1': 6.9053, 'a2': 5.2034, 'a3': 1.4379, 'a4': 1.5863, 'b1': 1.4679, 'b2': 22.215099, 'b3': 0.2536, 'b4': 56.172001, 'c': 0.8669},
 'specific_heat': 0.709999978542, 'symbol': 'S', 'van_der_waals_radii': 180.0,
 'volume': 28.0034999847, 'z': 16}


def make_element_subfolder_name(layer: Layer, ion: Ion,
                                precision: 'precisionLitType' = 'um') -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy.
    precision is units of the layer width, default = 'um' """

    if layer.name:
        element_list_str = layer.name
    else:
        element_list = []
        for (element, prop) in layer.elements.items():
            stoich = prop['stoich']
            # print(element.symbol, stoich)
            if stoich == 1.0:
                element_str = element.symbol
            elif stoich.is_integer():
                element_str = f'{element.symbol}{stoich:.0f}'
            else:
                element_str = f'{element.symbol}{stoich:.2f}'
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
    print(data_subfolder_name)
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


def get_energy_damage_array(results: Results) -> np.ndarray:
    phon = results.phonons
    dx = max(phon.depth) / 1000.0  # units from pm to nm
    energy_damage = (phon.ions + phon.recoils) * dx
    damage_array_nm = np.array(phon.depth / 1000, energy_damage / phon.num_ions)
    return damage_array_nm


def calc_energy_damage(results: Results, depth: int = 0, units: precisionLitType = 'nm') -> floatArray:
    phon = results.phonons

    if units in ('nm', 'nano'):
        dx = max(phon.depth) / 100.0  # units from pm to nm
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        dx = max(phon.depth) / 100.0  # units from pm to Angstroms

    # if depth != 0:
        # dx = depth
    energy_damage: floatArray = (phon.ions + phon.recoils) * dx  # add the arrays and multiply
    return energy_damage


def plot_damage_energy(results: Results, ax: plt.axis, folder: Path, units: precisionLitType = 'nm') -> None:
    phon = results.phonons
    if units in ('nm', 'nano'):
        units_str = 'nm'
        depth = phon.depth / 1000.0
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        units_str = 'Angstroms'
        depth = phon.depth / 100.0
    energy_damage = calc_energy_damage(results, units)
    ax.plot(depth, energy_damage / phon.num_ions, label='{}'.format(folder))
    ax.set_xlabel(f'Depth [{units_str}]')
    ax.set_ylabel('eV')
    ax.legend()


def plot_damage_energy_per_ion(results: Results, folder: Path, units: precisionLitType = 'nm') -> Tuple[np.ndarray, np.ndarray]:
    phon = results.phonons
    if units in ('nm', 'nano'):
        units_str = 'nm'
        depth = phon.depth / 10
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        units_str = 'Angstroms'
        depth = phon.depth
    fig, ax = plt.subplots()
    energy_damage = calc_energy_damage(results, 300, units)
    energy_damage_sum: float = sum(energy_damage)
    # energy_damage_kev = energy_damage_sum / 1000

    # ax.plot(depth, energy_damage / phon.num_ions, label='{}'.format(folder))
    legend = f'{folder.name}, {energy_damage_sum} eV'
    ax.plot(depth, energy_damage / phon.num_ions, label='{}'.format(legend))
    ax.set_xlabel(f'Depth [{units_str}]')
    ax.set_ylabel('Collision damage [eV / ion]')
    ax.legend()
    fig.suptitle('Damage Energy vs. Depth', fontsize=15)
    fig.set_size_inches((10, 6))
    fig.savefig(os.path.join(image_out_dir, 'damagevsdepth_per_ion.png'), transparent=True)

    damage_depth_data = (phon.depth, energy_damage / phon.num_ions)
    return damage_depth_data


def plot_damage_energy_total(results: Results, folder: Path, units: precisionLitType = 'nm') -> Tuple[np.ndarray, np.ndarray]:
    phon = results.phonons
    if units in ('nm', 'nano'):
        units_str = 'nm'
        depth = phon.depth / 10
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        units_str = 'Angstroms'
        depth = phon.depth
    fig, ax = plt.subplots()
    energy_damage = calc_energy_damage(results, 300, units)
    energy_damage_sum: float = sum(energy_damage)
    # energy_damage_kev = energy_damage_sum / 1000

    # ax.plot(depth, energy_damage / phon.num_ions, label='{}'.format(folder))
    legend = f'{folder.name}, {energy_damage_sum} eV'
    ax.plot(depth, energy_damage, label='{}'.format(legend))
    ax.set_xlabel(f'Depth [{units_str}]')
    ax.set_ylabel(f'Collision damage [eV] (total from {phon.num_ions} ions')
    ax.legend()
    fig.suptitle('Damage Energy vs. Depth', fontsize=15)
    fig.set_size_inches((10, 6))
    fig.savefig(os.path.join(image_out_dir, 'damagevsdepth_total.png'), transparent=True)

    damage_depth_data = (phon.depth, energy_damage)
    return damage_depth_data


def get_peak_damage_depth(results: Results) -> float:
    phon = results.phonons
    energy_damage = calc_energy_damage(results)

    peak_damage = max(energy_damage)
    peak_damage_index = np.argmax(energy_damage)
    peak_damage_depth = phon.depth[peak_damage_index] / 10  # convert to nm
    return peak_damage_depth


"""
def plot_multi_damage_energy(data: List[Tuple[np.ndarray, np.ndarray]], units: str) -> None:
    depth = data[0]
    damage_per_ion = data[1]
    if units in ('nm', 'nano'):
        units_str = 'nm'
        depth = data[0] / 10
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        units_str = 'Angstroms'
        depth = data[0]

    fig, ax = plt.subplots()
    ax.plot(depth, damage_per_ion, label='{}'.format())
    ax.set_xlabel(f'Depth [{units_str}]')
    ax.set_ylabel('Collision damage [eV]')
    ax.legend()
    fig.suptitle('Damage Energy vs. Depth', fontsize=15)
    fig.set_size_inches((10, 6))
    fig.savefig(os.path.join(image_out_dir, 'damagevsdepth.png'), transparent=True)
"""

if __name__ == "__main__":
    data_path: Final = Path(R'.\data')
    image_path: Final = Path(R'.\images')

    fig, ax = plt.subplots()
    ax.set_xlabel(f'Depth [nm]')
    ax.set_ylabel(f'Collision damage [eV] (total from 1000000 ions)')
    ax.legend()

    for folder in data_path.iterdir():
        results = Results(folder)
        phon = results.phonons
        energy_damage = calc_energy_damage(results)
        energy_damage_sum: float = sum(energy_damage)
        energy_damage_kev = energy_damage_sum / 1000
        print("{} Damage energy: {:.1f} eV".format(folder.name, energy_damage_sum))

        peak_damage = max(energy_damage)
        peak_damage_depth = get_peak_damage_depth(results)
        print(f"maximum damage is {peak_damage} at {peak_damage_depth}")
        legend = f'{folder.name}, {energy_damage_sum} eV'
        ax.plot(phon.depth / 10, energy_damage, label='{}'.format(legend))

    fig.suptitle('Damage Energy vs. Depth', fontsize=15)
    fig.set_size_inches((10, 6))
    axes = plt.gca()
    axes.set_xlim([0, 2000])
    fig.savefig(os.path.join(image_path, 'damagevsdepth_total.png'), transparent=True)
