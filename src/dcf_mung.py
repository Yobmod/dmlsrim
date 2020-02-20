
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from srim.output import Results

from typing import cast, Iterable
from typing_extensions import Final
from mytypes import floatArray, precisionLitType


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


def get_peak_damage_depth(results: Results) -> float:
    phon = results.phonons
    energy_damage = calc_energy_damage(results)
    # peak_damage = max(energy_damage)
    peak_damage_index: int = cast(int, np.argmax(energy_damage))
    peak_damage_depth: float = cast(float, phon.depth[peak_damage_index] / 10)  # convert to nm
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
        energy_damage_sum: float = sum(cast(Iterable[float], energy_damage))
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
