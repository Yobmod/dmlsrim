from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import os
from srim import Ion, Layer, Target
from srim.srim import TRIM
from srim.output import Results
from concurrent.futures import ThreadPoolExecutor, as_completed  # , ProcessPoolExecutor
import multiprocessing as mp
from time import sleep
from dataclasses import dataclass, asdict
from typing import Iterable, Sequence, Set, Union, List, Tuple, cast, Dict
from typing_extensions import Final, Literal
from mytypes import floatArray, precisionLitType
from matplotlib import use
use('Agg')  # NoQa


@dataclass
class SrimData:
    folder: Path  # folder the results is saved to
    ion: Ion
    target: Target
    num_ions: int
    damage_total: float
    damage_array: Tuple[floatArray, floatArray]

    def __post_init__(self) -> None:
        self.layers: List[Layer] = self.target.layers


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

    def as_dict(self) -> Dict[str, float]:
        # narrow str, Any to declared dtypes
        # get types from self.__annotation__ and make union in another function?
        return asdict(self)

    # def as_typdict(self) -> 'ElemTD':
        # return asdict(self)


# TODO see main.py for getting element classes. Need to convert to ElemClass or not? Use Dacite for convert via dict?
# or inherit from it?
elem_ce_dict = ElemClass(E_d=25.0, lattice=3.0, surface=4.23, atomic_num=58, atomic_mass=140.1)
elem_u_dict = {'E_d': 25.0, 'lattice': 3.0, 'surface': 5.42, 'atomic_num': 92, 'atomic_mass': 238.0}
elem_th_dict = {'E_d': 25.0, 'lattice': 3.0, 'surface': 5.93, 'atomic_num': 90, 'atomic_mass': 232.0}
elem_o_dict = ElemClass(E_d=28.0, lattice=3.0, surface=2.00, atomic_num=8, atomic_mass=15.99)
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
    # print(data_subfolder_name)
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
    """<depth> given in <units>"""
    phon = results.phonons

    # TODO is this needed? use depth if given, and if not, use phonon.depth.max. Units matter?
    if units in ('nm', 'nano'):
        dx = int(max(phon.depth) / 100.0)  # units from pm to nm
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        dx = int(max(phon.depth) / 100.0)  # units from pm to Angstroms
    else:
        raise ValueError

    print(f"phonon depth is int or array? {dx}")
    if depth > 0:
        dx = depth

    energy_damage: floatArray = (phon.ions + phon.recoils) * dx  # add the arrays and multiply
    return energy_damage  # up to depth if given otherwise all


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


def mung_srim(results: Results, depth: int = 0) -> float:
    energy_damage_sum: float = sum(cast(Iterable[float], calc_energy_damage(results)))
    # energy_damage_kev = energy_damage_sum / 1000
    print("Damage energy: {:.1f} eV".format(energy_damage_sum))
    # print("Damage energy: {:.1f} keV".format(energy_damage_kev))
    # TODO tuple total damage	max damage	depth of max / nm damage up to depth

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
    # pid = os.getpid() if using processpool
    data_out_dir = make_data_path(target.layers[0], ion, data_path)
    image_out_dir = data_out_dir  # make_image_path(target.layers[0], ion, data_path)
    print(f"{data_out_dir.name} started")  # using PID {pid}")

    result = run_srim(ion, target, data_out_dir, num_ions, srim_dir)
    damage_total = mung_srim(result)
    damage_array = plot_srim(result, image_out_dir)
    datum = SrimData(data_out_dir, ion, target, num_ions, damage_total, damage_array)

    end = datetime.now()
    duration = end - start
    print(f"{data_out_dir.name} done in {str(duration).split('.', 2)[0]}")  # " using PID {pid}")

    return datum


def create_ion_list(ion_name: Literal['H', 'He', 'Li'],
                    energy_list: Union[Sequence[int], Set[int]],
                    units: Literal['ev', 'kev', 'mev']
                    ) -> List[Ion]:
    ion_list = [Ion(f'{ion_name}', energy=x * 1000) for x in energy_list]
    return ion_list


def pool_srim(ions: Union[Sequence[Ion], Set[Ion]],
              target: Target, data_path: Path, num_ions: int, srim_dir: Path) -> List[SrimData]:  # List[SrimData]

    # with ProcessPoolExecutor(max_workers=mp.cpu_count() - 1) as ppexc:
    with ThreadPoolExecutor(max_workers=mp.cpu_count() * 5) as ppexc:

        """# using submit() and list comprehension
        SrimData_futures = [ppexc.submit(combined_srim,
                                         ion,
                                         target,
                                         data_path,
                                         num_ions=1_000_000,  # 1 million -> about 5 hours
                                         srim_dir=srim_executable_directory)
                            for ion in ions_He_list]
        """

        SrimData_futures = []
        for ion in ions:
            res = ppexc.submit(combined_srim,
                               ion,
                               target,
                               data_path,
                               num_ions,  # 1 million -> about 5 hours
                               srim_dir)
            sleep(1)
            SrimData_futures.append(res)

        """
        # alternate using map() and repeat().  # returns results in order done
        SrimData_futures = ppexc.map(combined_srim,
                                     [Ion('He', energy=1000000), Ion('He', energy=2000000)],
                                     repeat(target),
                                     repeat(data_path),
                                     repeat(1_000_000),  # 1 million -> about 5 hours
                                     repeat(srim_executable_directory))
        """

    Srim_data_list: List[SrimData] = [f.result() for f in as_completed(SrimData_futures)]

    print(f"{len(Srim_data_list)} jobs done")
    return Srim_data_list


if __name__ == "__main__":

    # set path to save data
    data_path: Final[Path] = Path(R'.\data\uo2pure')

    # poin to srim exec
    srim_executable_directory: Final[Path] = Path(R'C:\srim')

    # create list of ions from list of energies in keV
    energy_kev_list = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                       1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000]
    ions_He_list = create_ion_list('He', energy_kev_list, units='kev')

    # create target from list of layers
    layer_UO2_2um = Layer(
        {'U': {'stoich': 1.0, **elem_u_dict},
         'O': {'stoich': 2.0, **asdict(elem_o_dict)},
         },
        density=7.22, width=100_000.0, name='urania'
    )

    layer_SiO2_10um = Layer(
        {'Si': {'stoich': 1.0, **elem_si_dict},
         'O': {'stoich': 2.0, **asdict(elem_o_dict)},
         },
        density=2.65, width=100_000.0, name='silica'
    )

    target = Target([layer_UO2_2um])

    data_list = pool_srim(ions_He_list, target, data_path,
                          num_ions=1000_000, srim_dir=srim_executable_directory)
