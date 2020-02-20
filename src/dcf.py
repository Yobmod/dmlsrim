from __future__ import annotations

import os
from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

from srim import Ion, Layer, Target
from srim.srim import TRIM
from srim.output import Results

from concurrent.futures import ProcessPoolExecutor
from concurrent import futures
from dataclasses import dataclass, asdict
from typing import Iterable, Sequence, Set, Union, List, Tuple, cast, Dict
from typing_extensions import Final, Literal
from mytypes import floatArray, precisionLitType


@dataclass
class SrimData:
    folder: Path  # folder the results is saved to
    ion: Ion
    target: Target
    num_ions: int
    damage_total: float
    damage_array: Tuple[floatArray, floatArray]


@dataclass(frozen=True)
class ElemClass:
    atomic_num: int
    atomic_mass: float
    E_d: float = 25.0
    lattice: float = 0.0
    surface: float = 3.0

    def __post_init__(self) -> None:
        if self.E_d <= 0:
            raise ValueError('Invalid E_d (negative)')
        assert self.lattice >= 0


# TODO see main.py for getting element classes. Need to convert to ElemClass or not? Use Dacite for convert via dict?
# or inherit from it?
elem_ce_dict = ElemClass(E_d=25.0, lattice=3.0, surface=4.23, atomic_num=58, atomic_mass=140.1)
elem_u_dict: Dict = {'E_d': 25.0, 'lattice': 3.0, 'surface': 5.42, 'atomic_num': 92, 'atomic_mass': 238.0}
elem_th_dict = {'E_d': 25.0, 'lattice': 3.0, 'surface': 5.93, 'atomic_num': 90, 'atomic_mass': 232.0}
elem_o_dict = {'E_d': 28.0, 'lattice': 3.0, 'surface': 2.00, 'atomic_num': 8, 'atomic_mass': 15.99}
elem_si_dict = {'E_d': 15.0, 'lattice': 2.0, 'surface': 4.70, 'atomic_num': 14, 'atomic_mass': 28.08}
elem_ti_dict = {'E_d': 28.0, 'lattice': 3.0, 'surface': 2.00, 'atomic_num': 22, 'atomic_mass': 15.99}


def make_element_subfolder_name(layer: Layer, ion: Ion,
                                precision: precisionLitType = 'um') -> Path:
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


def calc_energy_damage(results: Results, units: precisionLitType = 'nm', depth: int = 0) -> floatArray:
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


def plot_damage_energy_per_ion(results: Results, folder: Path, units: precisionLitType = 'nm') -> Tuple[floatArray, floatArray]:
    phon = results.phonons
    if units in ('nm', 'nano'):
        units_str = 'nm'
        depth = phon.depth / 10
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        units_str = 'Angstroms'
        depth = phon.depth
    fig, ax = plt.subplots()
    energy_damage: floatArray = calc_energy_damage(results, units, 300)
    energy_damage_sum = sum(energy_damage)
    # energy_damage_kev = energy_damage_sum / 1000

    # ax.plot(depth, energy_damage / phon.num_ions, label='{}'.format(folder))
    legend = f'{folder.name}, {energy_damage_sum} eV'
    ax.plot(depth, energy_damage / phon.num_ions, label='{}'.format(legend))
    ax.set_xlabel(f'Depth [{units_str}]')
    ax.set_ylabel('Collision damage [eV / ion]')
    ax.legend()
    fig.suptitle('Damage Energy vs. Depth', fontsize=15)
    fig.set_size_inches((10, 6))
    fig.savefig(os.path.join(folder, 'damagevsdepth_per_ion.png'), transparent=True)

    damage_depth_data = (phon.depth, cast(floatArray, energy_damage / phon.num_ions))
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
    energy_damage: floatArray = calc_energy_damage(results, units, 300)
    energy_damage_sum: float = sum(cast(Iterable[float], energy_damage))
    # energy_damage_kev = energy_damage_sum / 1000

    # ax.plot(depth, energy_damage / phon.num_ions, label='{}'.format(folder))
    legend = f'{folder.name}, {energy_damage_sum} eV'
    ax.plot(depth, energy_damage, label='{}'.format(legend))
    ax.set_xlabel(f'Depth [{units_str}]')
    ax.set_ylabel(f'Collision damage [eV] (total from {phon.num_ions} ions')
    ax.legend()
    fig.suptitle('Damage Energy vs. Depth', fontsize=15)
    fig.set_size_inches((10, 6))
    fig.savefig(os.path.join(folder, 'damagevsdepth_total.png'), transparent=True)

    damage_depth_data = (phon.depth, energy_damage)
    return damage_depth_data


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


def run_srim(ion: Ion,
             target: Target,
             data_out_dir: Path,
             num_ions: int, srim_dir: Path) -> Results:
    # use layer, data_path and iob to create out_dir
    # run trim, return out_dir and result
    # copy result to out_dir from srim_dir
    trim = TRIM(target, ion, number_ions=num_ions, calculation=1)  # 1 million -> about 5 hours
    results = trim.run(srim_dir)
    TRIM.copy_output_files(srim_dir, data_out_dir)
    print(f'{ion.symbol}-{ion.energy/1000}kev done')
    return results


def mung_srim(results: Results) -> float:
    energy_damage_sum: float = sum(cast(Iterable[float], calc_energy_damage(results)))
    # energy_damage_kev = energy_damage_sum / 1000
    print("Damage energy: {:.1f} eV".format(energy_damage_sum))
    # print("Damage energy: {:.1f} keV".format(energy_damage_kev))
    return energy_damage_sum


def plot_srim(results: Results, image_out_dir: Path) -> Tuple[floatArray, floatArray]:
    # damage_per_ion_vs_depth = plot_damage_energy_per_ion(results, data_out_dir, units='nm')
    damage_per_ion_vs_depth = plot_damage_energy_total(results, image_out_dir, units='nm')
    return damage_per_ion_vs_depth


def combined_srim(ion: Ion,
                  target: Target,
                  data_path: Path,
                  num_ions: int,
                  srim_dir: Path) -> SrimData:
    # run ions in list against layer and datapath
    # get out_dir and result
    # create list of folders and list of results
    start = datetime.now()
    data_out_dir = make_data_path(target[0], ion, data_path)
    image_out_dir = data_out_dir / 'images'
    result = run_srim(ion, target, data_out_dir, num_ions, srim_dir)
    damage_total = mung_srim(result)
    damage_array = plot_srim(result, image_out_dir)
    datum = SrimData(data_out_dir, ion, target, num_ions, damage_total, damage_array)
    end = datetime.now()
    duration = end - start
    print(duration)

    return datum


if __name__ == "__main__":

    data_path: Final[Path] = Path(R'.\data\x')
    srim_executable_directory: Final[Path] = R'C:\srim'
    energy_kev_list = [2000, 2500, 3000, 4000, 5000]

    def create_ion_list(ion_name: Literal['H', 'He'],
                        energy_list: Union[Sequence[int], Set[int]]) -> List[Ion]:
        ion_list = [Ion(f'ion_name', energy=x * 1000) for x in energy_list]
        return ion_list

    ions_He_list = create_ion_list('He', energy_kev_list)

    layer_CeO2_2um = Layer(
        {'Ce': {'stoich': 1.0, **asdict(elem_ce_dict)},
         'O': {'stoich': 2.0, **elem_o_dict},
         },
        density=7.22, width=20_000.0, name='ceria'
    )

    layer_SiO2_10um = Layer(
        {'Si': {'stoich': 1.0, **elem_si_dict},
         'O': {'stoich': 2.0, **elem_o_dict},
         },
        density=2.65, width=100_000.0, name='silica'
    )

    target = Target([layer_CeO2_2um, layer_SiO2_10um])

    def pool_srim(ion_list: Union[Sequence[Ion], Set[Ion]],
                  target: Target, data_path: Path, num_ions: int, srim_dir: Path) -> None:  # List[SrimData]
        ...

    with ProcessPoolExecutor(max_workers=6) as ppexc:
        SrimData_futures = [ppexc.submit(combined_srim,
                                         ion,
                                         target,
                                         data_path,
                                         num_ions=1_000_000,  # 1 million -> about 5 hours
                                         srim_dir=srim_executable_directory)
                            for ion in ions_He_list]

        # Srimdata_futuress = ppexc.map(combined_srim,
        # ions_He_list,
        # repeat(target),
        # repeat(data_path),
        # num_ions=1_000_000)  # returns results in order done

    for f in futures.as_completed(SrimData_futures):
        print('main: result: {}'.format(f.result()))

    Srim_data_list: List[SrimData] = [f.result() for f in futures.as_completed(SrimData_futures)]
