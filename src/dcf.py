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


from typing import Iterable, Sequence, Set, Union, List, Tuple, cast, Dict, NamedTuple
from typing_extensions import Literal, TypedDict
from mytypes import floatArray, precisionLitType
from matplotlib import use
use('Agg')  # NoQa


class PydanticConfig:
    arbitrary_types_allowed = True


@dataclass(config=PydanticConfig)
class SrimData:
    folder: Path  # folder the results is saved to
    ion: Ion
    num_ion: int
    target: Target
    damage_total: float
    damage_array: Tuple[floatArray, floatArray]

    def __post_init__(self) -> None:

        ...

    def __post_init_post_parse__(self) -> None:
        self.results = Results(self.folder)
        import re

        if not self.ion:
            self.ion = self.results.ioniz.ion
        self.num_ions: int = self.results.ioniz.num_ions

        if not self.target:
            with open(R".\data\ceria_on_silica\ceria_2um_He@400keV\tdata.txt", 'r') as f:
                f.read()
                """===============Target material =======================
                Layer 1 """
                match_target = re.search(r'(?<=====\r\n)Layer\s+\d+\s+:.*?(?=====)', f.read(), re.DOTALL)
                #match_target = re.search(r'(?<=====\r\n)Layer\s+\d+\s+:.*?(?=====)', f.read(), re.DOTALL)

                if match_target:
                    print(match_target.group(0))
                else:
                    print("target not found")
                # out = output.SRIM_Output()
                # output.SRIM_Output._read_target(out, f.read())

            # self.target = Target.
        # self.layers: List[Layer] = self.target.layers


class ElemTD(TypedDict):
    atomic_num: int
    atomic_mass: float
    E_d: float
    lattice: float
    surface: float


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

    def as_typdict(self) -> ElemTD:
        #dic = {str(k): float(v) for k, v in asdict(self).items()}
        #ret: ElemTD = dic
        #ret = cast(ElemTD, asdict(self))
        return ElemTD(atomic_num=self.atomic_num,
                      atomic_mass=self.atomic_mass,
                      E_d=self.E_d,
                      lattice=self.lattice,
                      surface=self.surface)


class DamageStats(NamedTuple):
    total: float
    max_damage: float
    max_index: int
    max_depth: float


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


def get_depth_damage_array(results: Results, units: str = 'nm') -> np.ndarray[float]:
    """get array of [0] depths in nm and [damage] for whole target"""
    if units in ('nm', 'nano'):
        ratio_A_to_units = 10
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        ratio_A_to_units = 1
    else:
        raise ValueError

    phon = results.phonons
    dx = max(phon.depth) / 100
    energy_damage = np.array((phon.ions + phon.recoils) * dx)
    depth_array = np.array(phon.depth / ratio_A_to_units)
    damage_array_nm: np.ndarray[float] = np.stack((depth_array, energy_damage))
    return damage_array_nm


def trunc_depth_damage_array(results: Results, units: precisionLitType = 'nm', depth: int = 0) -> floatArray:
    """Get list of damage up to given depth. <depth> given in <units>"""
    depth_damage_array = get_depth_damage_array(results, units=units)

    if depth > 0:
        # print(depth_damage_array[0][depth_damage_array[0][:] <= depth])
        depth_damage = depth_damage_array[:, depth_damage_array[0][:] <= depth]
    else:
        depth_damage = depth_damage_array[:]
    return depth_damage  # up to depth if given otherwise all


def get_damage_array(results: Results, units: precisionLitType = 'nm', depth: int = 0) -> floatArray:
    depth_damage = trunc_depth_damage_array(results, units=units, depth=depth)
    damage_array = depth_damage[1]
    return damage_array


def get_damage_stats(results: Results, units: precisionLitType = 'nm', depth: int = 0) -> DamageStats:
    array = trunc_depth_damage_array(results, units=units, depth=depth)
    total_damage: int = int(sum(cast(Iterable[float], array[1])))
    max_damage: int = int(max(array[1]))
    max_ind: int = np.argmin(array[1])
    depth_of_max: float = array[0][max_ind]
    return DamageStats(total_damage, max_damage, max_ind, depth_of_max)


def plot_damage_multi(results: List[Results],
                      save_dir: Path,
                      units: precisionLitType = 'nm',
                      depth: int = 0
                      ) -> None:
    # phon = results.phonons
    if units in ('nm', 'nano'):
        units_str = 'nm'
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        units_str = 'Angstroms'

    if depth > 0:
        pass
        # add doted line at depth
    if isinstance(results, Results):
        results = [results]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(f'Depth [{units_str}]')
    ax.set_ylabel('Collision damage [eV]')

    for res in results:
        depth_damage_array = trunc_depth_damage_array(res, units=units, depth=depth)
        damage_stats = get_damage_stats(res, units=units, depth=depth)
        ion_name = res.ioniz.ion.symbol
        ion_energy = int(res.ioniz.ion.energy / 1000)
        legend = f'{ion_name} @ {ion_energy} keV, damage {damage_stats.total} eV'
        ax.plot(depth_damage_array[0], depth_damage_array[1], label='{}'.format(legend))

    ax.legend()
    fig.suptitle('Damage Energy vs. Depth', fontsize=15)
    fig.set_size_inches((10, 6))
    fig.savefig(os.path.join(save_dir, 'damagevsdepth_multi.png'), transparent=True)
    # return fig


def plot_damage_multi_from_path(data_parent: Path,
                                units: precisionLitType = 'nm',
                                depth: int = 0,
                                ) -> None:
    loaded_data = [Results(dp) for dp in data_parent.iterdir() if dp.is_dir()]
    plot_damage_multi(loaded_data, data_parent, units=units, depth=depth)


def plot_damage_energy_per_ion(results: Results, folder: Path, units: precisionLitType = 'nm') -> Tuple[floatArray, floatArray]:
    phon = results.phonons
    if units in ('nm', 'nano'):
        units_str = 'nm'
        depth = phon.depth / 10
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        units_str = 'Angstroms'
        depth = phon.depth
    fig, ax = plt.subplots()
    energy_damage: floatArray = get_damage_array(results, units, 0)
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


def plot_damage_energy_total(results: Results, folder: Path, units: precisionLitType = 'nm') -> Tuple[np.ndarray[float], np.ndarray[float]]:
    phon = results.phonons
    if units in ('nm', 'nano'):
        units_str = 'nm'
        depth = phon.depth / 10
    elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
        units_str = 'Angstroms'
        depth = phon.depth
    fig, ax = plt.subplots()
    energy_damage: floatArray = get_damage_array(results, units, 0)
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
    energy_damage_sum: float = sum(cast(Iterable[float], get_damage_array(results)))
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
    pid = os.getpid()  # if using processpool
    data_out_dir = make_data_path(target.layers[0], ion, data_path)
    image_out_dir = data_out_dir  # make_image_path(target.layers[0], ion, data_path)
    print(f"{data_out_dir.name} started) using PID {pid}")

    result = run_srim(ion, target, data_out_dir, num_ions, srim_dir)
    damage_total = mung_srim(result)
    damage_array = plot_srim(result, image_out_dir)
    datum = SrimData(data_out_dir, ion, num_ions, target, damage_total, damage_array)

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
    with ProcessPoolExecutor(max_workers=mp.cpu_count() * 5) as ppexc:

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


def pickle_srim(srimdata: Union[SrimData, Sequence[SrimData]]) -> None:
    # sequence or iterable?
    if isinstance(srimdata, SrimData):
        srimdata = [srimdata]

    for srim_x in srimdata:
        datapath = srim_x.folder / "result.pkl"

        with open(datapath, "w+b") as pkl_f:
            pickle.dump(srim_x, pkl_f)
            print(f"Data pickled to {datapath}")


def json_srim(srimdata: List[SrimData]) -> None:
    data_ref = srimdata if isinstance(srimdata, SrimData) else srimdata[0]
    datapath = data_ref.folder.parent / "result.json"

    with open(datapath, "w+") as json_f:
        json.dump(srimdata, json_f)
        print(f"Data save as json to {datapath}")
