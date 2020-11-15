from __future__ import annotations
import json
import pickle
import matplotlib as mpl
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
from pydantic import BaseModel
from tabulate import tabulate

import typing as typ
from typing import cast, Iterable, Sequence, Set, Union, List, Tuple, Dict, NamedTuple, Optional
from typing_extensions import Literal, TypedDict
from mytypes import floatArray, precisionLitType
mpl.use('tkAgg')  # NoQa


class PydanticConfig:
    arbitrary_types_allowed = True


@dataclass(config=PydanticConfig)
class SrimData:
    folder: Path  # folder the results is saved to
    ion: Ion
    num_ion: int
    target: Target
    damage_total: float
    damage_array: floatArray

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
                # match_target = re.search(r'(?<=====\r\n)Layer\s+\d+\s+:.*?(?=====)', f.read(), re.DOTALL)

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
        # dic = {str(k): float(v) for k, v in asdict(self).items()}
        # ret: ElemTD = dic
        # ret = cast(ElemTD, asdict(self))
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
elem_u_dict: ElemTD = {'E_d': 25.0, 'lattice': 3.0, 'surface': 5.42, 'atomic_num': 92, 'atomic_mass': 238.0}
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

    if precision in ['um', 'micro']:
        layer_width = layer_width_um
    elif precision in ['nm', 'nano']:
        layer_width = layer_width_nm
    else:
        layer_width = layer.width

    # print(data_subfolder_name)
    return Path(f"{element_list_str}_{layer_width}_{ion.symbol}@{ion_energy_kev}")


cwd_path = Path(R'.')


def make_data_path(layer: Layer,
                   ion: Ion,
                   data_path: Union[Path, str] = cwd_path / 'data',
                   precision: precisionLitType = 'um',
                   ) -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy
    data_path default = '.\\data'. precision is units of the layer width, default = 'um' """

    data_subfolder_name = make_element_subfolder_name(layer, ion, precision)
    output_directory: Path = Path(data_path) / data_subfolder_name
    output_directory.mkdir(parents=True, exist_ok=True)
    return output_directory


def make_image_path(layer: Layer, ion: Ion,
                    image_path: Union[Path, str] = cwd_path / 'images',
                    precision: precisionLitType = 'um') -> Path:
    """create a folder from layer elements and stoichiometries and ion type and energy
    data_path default = '.\\images'. precision is units of the layer width, default = 'um' """

    data_subfolder_name = make_element_subfolder_name(layer, ion, precision)
    outimage_directory: Path = Path(image_path) / data_subfolder_name
    outimage_directory.mkdir(parents=True, exist_ok=True)
    return outimage_directory


class SrimResults():
    """Class with methods:
    *__str__(self) -> str
    *_valid_units(self, units: precisionLitType) -> Tuple[Literal['nm', 'A'], int]
    *get_depth_damage_array(self) -> floatArray
    *trunc_depth_damage_array(self) -> floatArray
    *get_damage_array(self) -> floatArray
    *get_damage_stats(self) -> DamageStats:
    *plot_damage_energy_per_ion(self, folder: Path) -> None
    *plot_damage_energy_total(self, folder: Path) -> None

    """

    def __init__(self,
                 inp: Union[Path, Results],
                 target: Optional[Target] = None,
                 units: precisionLitType = 'nm',
                 depth: int = 0,
                 savepath: Union[Path, str] = R".",
                 ):

        try:
            if isinstance(inp, Path):
                self.results = Results(inp)
            elif isinstance(inp, Results):
                self.results = inp
        except FileNotFoundError:
            print("Data files not found. Check paths for TRIM output files")

        self._target = target
        self.depth = depth

        if isinstance(savepath, str):
            savepath = Path(savepath)
        self.savepath = savepath  # TODO make valid function that creates path

        self.units = self._valid_units(units)[0]
        self.ratio_A_to_units = self._valid_units(units)[1]

        self.ion = self.results.ioniz.ion
        self.ion_energy = self.ion.energy / 1_000  # in keV
        self.ion_num = self.results.phonons.num_ions

    def __str__(self) -> str:
        if self._target is not None and self.layers is not None:
            target_str = " ".join([str(lay) for lay in self.layers])
        else:
            target_str = "[Target not defined]"
        return f"Results class for target={target_str}, ion={self.ion.name} @ {self.ion_energy} kev. Analysis to {self.depth}{self.units}"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SrimResults):
            return NotImplemented
        else:
            return bool(self.ion.symbol == other.ion.symbol and self.ion_energy == other.ion_energy)

    def __lt__(self, other: SrimResults) -> bool:  # less than
        return bool((self.ion.symbol, self.ion_energy) < (other.ion.symbol, other.ion_energy))

    def _valid_units(self, units: precisionLitType) -> Tuple[Literal['nm', 'A'], int]:
        if units in ('nm', 'nano'):
            units = 'nm'
            ratio_A_to_units = 10
        elif units in ('a', 'A', 'angstrom', 'angstroms', 'Angstrom', 'Angstroms'):
            units = 'A'
            ratio_A_to_units = 1
        else:
            raise ValueError
        return units, ratio_A_to_units

    @property
    def target(self) -> Target:
        return self._target

    @target.setter
    def target(self, target: Target) -> None:
        self._target = target

    @property
    def layers(self) -> Optional[Target.layers]:
        if self._target is not None:
            return self._target.layers
        else:
            return None

    def _get_layer_depths(self) -> List[int]:
        if self.layers is not None:
            layer_depths = [layer.width / self.ratio_A_to_units for layer in self.layers]  # sourcery skip
            return layer_depths
        else:
            raise AttributeError("Layer depth cannot be calculated unless Target is given")

    def _get_total_thickness(self) -> int:
        return sum(self._get_layer_depths())

    def get_depth_damage_array(self) -> floatArray:
        """get array of [0] depths in nm and [damage] for whole target"""

        phon = self.results.phonons
        depths: np.ndarray[float] = phon.depth
        dx = depths.max() / 100  # ratio for eV/A to eV per measurement
        energy_damage = np.array((phon.ions + phon.recoils) * dx)
        depth_array = np.array(phon.depth / self.ratio_A_to_units)
        damage_array_nm: np.ndarray[float] = np.stack((depth_array, energy_damage))
        return damage_array_nm

    def trunc_depth_damage_array(self) -> floatArray:
        """Get list of damage up to given depth. <depth> and <units> defined at initialisation
        If no depth given, returns all depth / damage data"""

        depth_damage_array = self.get_depth_damage_array()

        if self.depth > 0:
            # print(depth_damage_array[0][depth_damage_array[0][:] <= depth])
            depth_damage = depth_damage_array[:, depth_damage_array[0][:] <= self.depth]
        else:
            depth_damage = depth_damage_array[:]
        return cast(floatArray, depth_damage)  # up to depth if given otherwise all

    def make_depth_damage_table(self) -> str:
        # sourcery skip: inline-immediately-returned-variable
        """convert data to table from printing"""
        headers = ['depth', 'damage']
        table = tabulate(self.trunc_depth_damage_array().transpose(),
                         headers, tablefmt="fancy_grid", colalign=("right",))
        return table

    def make_depth_damage_csv(self, filename: Union[Path, str] = "depth_damage.csv") -> str:
        filename = Path(f"{self.ion.symbol}@{self.ion_energy}kev.csv")
        filepath = self.savepath / filename
        np.savetxt(filepath,
                   self.trunc_depth_damage_array().transpose(),
                   header='depth\t damage',
                   footer=str(self),
                   delimiter="\t", comments="")
        return self.make_depth_damage_table()

    def get_damage_array(self) -> floatArray:
        """1D array of damage at given self.<depth>"""
        depth_damage = self.trunc_depth_damage_array()
        damage_array = depth_damage[1]
        return cast(floatArray, damage_array)

    def get_depth_array(self) -> floatArray:
        """1D array of damage at given self.<depth>"""
        depth_damage = self.trunc_depth_damage_array()
        damage_array = depth_damage[0]
        return cast(floatArray, damage_array)

    def get_index_from_depth(self, chosen_depth: int) -> int:
        depth_array = self.get_depth_array()
        # subtract chosen_depth from each, 0, then must be argmin. np.abs in case float precision gives -0.x
        idx: int = np.abs(depth_array - chosen_depth).argmin()
        return idx

    def get_damage_stats(self) -> DamageStats:
        """Get stats to given <depth> as namedTuple. total_damage / max_damage / max_ind / depth_of_max"""
        array = self.trunc_depth_damage_array()
        total_damage = int(array[1].sum())
        max_damage = int(array[1].max())
        max_ind = np.argmin(array[1])
        depth_of_max = float(array[0][max_ind])
        return DamageStats(total_damage, max_damage, max_ind, depth_of_max)

    def _create_damage_depth_fig(self,
                                 depth_marker: int = 0,
                                 plot_type: Literal['total', 'per_ion'] = 'total'
                                 ) -> Tuple[mpl.figure.Figure, mpl.axes.Axes]:

        if plot_type == 'total':
            y_units = f"[eV] (total from {self.ion_num} ions)"
        elif plot_type == 'per_ion':
            y_units = '[eV / ion]'
        else:
            y_units = '[eV]'

        fig, ax = plt.subplots()
        fig.suptitle('Damage Energy vs. Depth', fontsize=15)
        fig.set_size_inches((10, 6))

        ax.set_xlabel(f'Depth [{self.units}]')
        ax.set_ylabel(f'Collision damage {y_units}')

        if depth_marker != 0:
            plt.axvline(x=depth_marker, ymin=0.0, ymax=200.0, linestyle='--', linewidth=1,
                        color='g', label=f'Layer boundary, d = {depth_marker} nm')
        return fig, ax

    def _add_damage_depth_line(self,
                               ax: Optional[mpl.axes.Axes] = None,
                               plot_type: Literal['total', 'per_ion'] = 'total'
                               ) -> mpl.axes.Axes:
        if ax is None:
            ax = plt.gca()

        energy_damage: floatArray = self.get_damage_array()
        energy_damage_sum: float = sum(cast(Iterable[float], energy_damage))
        limit = len(energy_damage)  # number of data point up to self.depth

        phon = self.results.phonons
        depth_array = phon.depth[:limit] / self.ratio_A_to_units

        if plot_type == 'per_ion':
            damage = np.divide(energy_damage, phon.num_ions)
        elif plot_type == 'total':
            damage = energy_damage
        else:
            raise NameError('plot_type must be "total" or "per_ion"')

        legend = f'{self.ion.symbol} @ {self.ion_energy: .0f} keV, d = {round(energy_damage_sum, 0)} eV'
        ax.plot(depth_array, damage, label=f'{legend}')
        ax.legend()
        return ax

    def plot_damage_depth(self,
                          filename: Union[Path, str, None] = None,
                          depth_marker: int = 0,
                          plot_type: Literal['total', 'per_ion'] = 'total'
                          ) -> None:
        """plot damage energy against depth, up to given <depth>. Add vertical line at <depth_marker>"""

        fig, ax = self._create_damage_depth_fig(plot_type=plot_type)
        self._add_damage_depth_line(ax=ax, plot_type=plot_type)

        depth_txt = f"_to_{self.depth}{self.units}" if self.depth != 0 else ""

        if filename is None:
            filename = Path(f'{self.ion.symbol}@{self.ion_energy}_damagevsdepth_{plot_type}{depth_txt}.png')
        elif isinstance(filename, (str, Path)):
            filename = Path(filename)

        print(self.savepath / filename)
        fig.savefig(self.savepath / filename, transparent=True)
        # plt.show()
        plt.close(fig)

    def plot_srim(self,
                  filename: Union[Path, str, None] = None,
                  depth_marker: int = 0,
                  total: bool = True,
                  per_ion: bool = True,
                  ) -> None:

        if total:
            self.plot_damage_depth(filename, depth_marker=depth_marker, plot_type='total')
        if per_ion:
            self.plot_damage_depth(filename, depth_marker=depth_marker, plot_type='per_ion')


class MultiSrimResults():

    def __init__(self,
                 inp: Union[str, Path, typ.Sequence[Path], typ.Sequence[Results], typ.Sequence[SrimResults]],
                 target: Optional[Target] = None,
                 units: precisionLitType = 'nm',
                 depth: int = 0,
                 savepath: Union[Path, str] = R".\output",
                 ):

        self._target = target
        self.depth = depth
        self.units = units
        self.savepath = Path(savepath)
        self.result_list: Sequence[SrimResults]
        try:
            if isinstance(inp, (Path, str)):
                inp = Path(inp)
                self.result_list = [SrimResults(x, depth=self.depth) for x in inp.iterdir() if inp.is_dir()]
            elif isinstance(inp, typ.Sequence) and all(isinstance(x, (Path, Results)) for x in inp):
                self.result_list = [SrimResults(cast(Union[Results, Path], x), depth=self.depth) for x in inp]
            elif isinstance(inp, typ.Sequence) and all(isinstance(x, SrimResults) for x in inp):
                self.result_list = cast(Sequence[SrimResults], inp)

            else:
                raise TypeError("MultiSrimResults() takes Path or Sequence")
        except FileNotFoundError:
            print("Data files not found. Check paths for TRIM output files")

        self.proxy = self.result_list[0]

    def plot_damage_multi(self,
                          depth_marker: int = 0,
                          plot_type: Literal['total', 'per_ion'] = 'total',
                          ) -> plt.figure.Figure:
        """
        plot_damage_multi([res.results], res.path, res.units, res.depth)
        """

        proxy_res = self.result_list[0]

        fig, ax = proxy_res._create_damage_depth_fig(depth_marker, plot_type)

        for result in self.result_list:
            result._add_damage_depth_line(ax, plot_type=plot_type)
        fig.savefig(self.savepath / 'damagevsdepth_multi.png', transparent=True)
        plt.close(fig)
        return fig

    def get_depth_damage_arrays(self) -> List[floatArray]:
        return [res.get_depth_damage_array() for res in self.result_list]

    def get_damage_total_energy_array(self, depth_marker: int = 0) -> floatArray:
        depth_damage_arrays = self.get_depth_damage_arrays()
        if depth_marker > 0:
            # print(depth_damage_array[0][depth_damage_array[0][:] <= depth])
            depth_damages = [dd[:, dd[0][:] <= depth_marker] for dd in depth_damage_arrays]
        else:
            depth_damages = [dd[:] for dd in depth_damage_arrays]

        dmges = [depth_damage[1] for depth_damage in depth_damages]
        total_dmg = [sum(d) for d in dmges]
        print(total_dmg)
        ion_energies = [res.ion_energy for res in self.result_list][:len(total_dmg)]

        ion_energy_total_dmg: floatArray = np.stack((ion_energies, total_dmg))
        return ion_energy_total_dmg

    def plot_max_and_total_dmg(self,
                               depth_marker: int = 0,
                               x_max: int = 0,
                               depth_max: int = 0,
                               damage_max: int = 0,
                               title: bool = False,
                               ) -> mpl.figure.Figure:

        fig, ax = plt.subplots()
        data = self.get_damage_total_energy_array(depth_marker)

        ax.set_xlabel('Ion Energy [ev]')
        ax.set_ylabel(f'Collision damage [eV] (total from {self.proxy.ion_num} ions)')
        legend = f'Total Damage in {depth_marker} nm'
        line = ax.plot(data[0], data[1], label='{}'.format(legend))
        # ax.legend()

        if x_max > 0:
            plt.xlim((0, x_max))

        if damage_max > 0:
            plt.ylim((0, damage_max))

        if title:
            fig.suptitle('Total Damage Energy and depth of maximum damage vs. Ion Energy', fontsize=15)

        proxy = self.result_list[0]
        total_thickness = proxy._get_total_thickness()

        sex_ax = ax.twinx()
        sex_ax.set_ylabel(f'Depth of max damage [{self.units}]', rotation=-90, labelpad=5)

        legend2 = 'Depth of max damage'
        ion_energies = [res.ion_energy for res in self.result_list]
        max_damage_depth = [res.get_damage_stats().max_depth for res in self.result_list]
        line2 = sex_ax.plot(
            ion_energies[: len(max_damage_depth)],
            max_damage_depth, color='g', label='{}'.format(legend2))

        if depth_marker != 0:
            max_max_depth = max(max_damage_depth)
            depth_marker_txt = f"_marked_{depth_marker}nm"
            d_line = plt.axhline(y=depth_marker, xmin=0.0, xmax=max_max_depth*1.2, linestyle='--', linewidth=1,
                                 color='g', label=f'Layer boundary, d = {depth_marker} nm')
        else:
            depth_marker_txt = ""
            d_line = None

        if depth_max > 0:
            plt.ylim((0, depth_max))
        else:
            plt.ylim((0, total_thickness))

        lines = line + line2 + [d_line]
        labels = [lin.get_label() for lin in lines]
        sex_ax.legend(lines, labels, loc='upper right')
        # 'upper left', 'upper right', 'lower left', 'lower right'
        # 'upper center', 'lower center', 'center left', 'center right'
        # 'center', 'best'
        fig.set_size_inches((10, 6))
        fig.savefig(self.savepath / f'total_dmg_and_depth_of_max_vs_ion_energy{depth_marker_txt}.png', transparent=True)

        return fig

    def plot_damage_total_energy(self, x_max: int = 0, y_max: int = 0, title: bool = False) -> mpl.figure.Figure:

        fig, ax = plt.subplots()
        data = self.get_damage_total_energy_array()

        ax.set_xlabel('Ion Energy [ev]')
        ax.set_ylabel(f'Collision damage [eV] (total from {self.proxy.ion_num} ions)')
        legend = f'{self.proxy.ion.symbol}'
        ax.plot(data[0], data[1], label='{}'.format(legend))
        ax.legend()

        if x_max > 0:
            plt.xlim((0, x_max))

        if y_max > 0:
            plt.ylim((0, y_max))

        if title:
            fig.suptitle('Total Damage Energy vs. Ion Energy', fontsize=15)
        fig.set_size_inches((10, 6))
        fig.savefig(self.savepath / f'total_damage_vs_ion_energy_to {self.depth}nm.png', transparent=True)
        return fig, ax

    def plot_max_dmgdepth_energy(self,
                                 depth_marker: int = 0,
                                 x_max: int = 0,
                                 y_max: int = 0,
                                 ) -> Tuple[mpl.figure.Figure, mpl.axes.Axes]:
        proxy = self.result_list[0]
        total_thickness = proxy._get_total_thickness()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Ion Energy [ev]')
        ax.set_ylabel(f'Depth of max damage [{self.units}]')

        legend = f'{self.proxy.ion.symbol}'
        ion_energies = [res.ion_energy for res in self.result_list]
        max_damage_depth = [res.get_damage_stats().max_depth for res in self.result_list]
        ax.plot(ion_energies[:len(max_damage_depth)], max_damage_depth, label='{}'.format(legend))

        if depth_marker != 0:
            max_max_depth = max(max_damage_depth)
            depth_marker_txt = f"_marked_{depth_marker}nm"
            plt.axhline(y=depth_marker, xmin=0.0, xmax=max_max_depth*1.2, linestyle='--', linewidth=1,
                        color='g', label=f'Layer boundary, d = {depth_marker} nm')
        else:
            depth_marker_txt = ""

        if x_max > 0:
            plt.xlim((0, x_max))

        if y_max > 0:
            plt.ylim((0, y_max))
        else:
            plt.ylim((0, total_thickness))

        ax.legend()
        fig.suptitle('Depth of maximum damage vs. Ion energy', fontsize=15)
        fig.set_size_inches((10, 6))
        fig.savefig(self.savepath / f'depth_of_max_vs_ion_energy{depth_marker_txt}.png', transparent=True)
        return fig, ax

    def plot_surface_damage(self, chosen_depth: int) -> None:
        # calc and plot damage at surface

        idx = self.proxy.get_index_from_depth(chosen_depth)
        dmg_nm = [res.get_damage_array()[idx] for res in self.result_list]
        ion_energies = [res.ion_energy for res in self.result_list]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Ion Energy [ev]')
        ax.set_ylabel(f'Collision damage [eV] at {chosen_depth} nm (total from {self.proxy.ion_num} ions)')
        legend = f'{self.proxy.ion.symbol}'

        ax.plot(ion_energies, dmg_nm, label='{}'.format(legend))
        ax.legend()
        fig.suptitle(f'Total Damage Energy at {chosen_depth} nm vs. Ion Energy', fontsize=15)
        fig.set_size_inches((10, 6))
        fig.savefig(self.savepath / f'damage_{chosen_depth: .0f}nm_vs_ion_energy.png', transparent=True)
        # plt.close(fig)


############################################################################################################
# TODO Make wrapper class that takes parent folder


def plot_damage_multi(results: Sequence[Results],
                      save_dir: Path,
                      units: precisionLitType = 'nm',
                      depth: int = 0
                      ) -> None:
    """
    plot_damage_multi([res.results], res.path, res.units, res.depth)
    """
    if depth > 0:
        pass
        # add doted line at depth
    if isinstance(results, Results):
        results = [results]

    res = SrimResults(results[0])
    units_str = res._valid_units(units)[0]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(f'Depth [{units_str}]')
    ax.set_ylabel('Collision damage [eV]')

    for result in results:
        res = SrimResults(result)
        depth_damage_array = res.trunc_depth_damage_array()
        damage_stats = res.get_damage_stats()
        ion_symbol = res.ion.symbol
        ion_energy = int(res.ion.energy / 1000)
        legend = f'{ion_symbol} @ {ion_energy} keV, damage {damage_stats.total} eV'
        ax.plot(depth_damage_array[0], depth_damage_array[1], label='{}'.format(legend))

    ax.legend()
    fig.suptitle('Damage Energy vs. Depth', fontsize=15)
    fig.set_size_inches((10, 6))
    fig.savefig(os.path.join(save_dir, 'damagevsdepth_multi.png'), transparent=True)
    # return fig


def get_ion_energy_damage_total(results: Sequence[SrimResults]) -> floatArray:
    ion_energies = [data.ion_energy for data in results]
    total_dmg = [data.get_damage_stats().total for data in results]
    ion_energy_total_dmg: floatArray = np.stack((ion_energies, total_dmg))
    return ion_energy_total_dmg


def plot_ion_energy_damage(results: Sequence[SrimResults],
                           save_dir: Path,
                           ) -> None:
    data = get_ion_energy_damage_total(results)

    fig = plt.figure()
    legend = ""

    ax = fig.add_subplot(111)
    ax.plot(data[0], data[1], label='{}'.format(legend))
    ax.set_xlabel('Ion Energy [keV]')
    ax.set_ylabel('Collision damage [eV]')
    ax.legend()

    fig.suptitle('Damage Energy vs. Ion Energy', fontsize=15)
    fig.set_size_inches((10, 6))
    fig.savefig(save_dir / 'damagevsion.png', transparent=True)


def plot_ion_energy_maxdmg_depth(results: Sequence[SrimResults],
                                 save_dir: Path,
                                 ) -> None:
    ion_energies = [data.ion_energy for data in results]
    dmg_depth = [data.get_damage_stats().max_depth for data in results]
    ion_energy_dmg_depth: floatArray = np.stack((ion_energies, dmg_depth))
    units = results[0].units
    fig = plt.figure()
    legend = ""

    ax = fig.add_subplot(111)
    ax.plot(ion_energy_dmg_depth[0], ion_energy_dmg_depth[1], label='{}'.format(legend))
    ax.set_xlabel('Ion Energy [keV]')
    ax.set_ylabel(f'Depth of max damage [{units}]')
    ax.legend()

    fig.suptitle('Depth of Maximum Damage Energy vs. Ion Energy', fontsize=15)
    fig.set_size_inches((10, 6))
    fig.savefig(save_dir / 'depthvsion.png', transparent=True)


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
    srim_res = SrimResults(result)
    damage_stats = srim_res.get_damage_stats()
    damage_total = damage_stats.total
    damage_array = srim_res.get_depth_damage_array()
    srim_res.plot_srim(image_out_dir)
    datum = SrimData(data_out_dir, ion, num_ions, target, damage_total, damage_array)

    end = datetime.now()
    duration = end - start
    print(f"{data_out_dir.name} done in {str(duration).split('.', 2)[0]}")  # " using PID {pid}")

    return datum


def create_ion_list(ion_name: Literal['H', 'He', 'Li'],
                    energy_list: Union[Sequence[int], Set[int]],
                    units: Literal['ev', 'kev', 'mev'],
                    ) -> List[Ion]:
    ion_list = [Ion(f'{ion_name}', energy=x * 1000) for x in energy_list]  # sourcery skip
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
