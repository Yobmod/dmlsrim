from pathlib import Path
from typing_extensions import Final
from srim import Layer, Target
from .dcf import create_ion_list, pool_srim, pickle_srim, elem_u_dict, elem_o_dict, elem_si_dict

if __name__ == "__main__":

    # set path to save data
    data_path: Final[Path] = Path(r".\data\uo2pure")

    # poin to srim exec
    srim_executable_directory: Final[Path] = Path(r"C:\srim")

    # create list of ions from list of energies in keV
    energy_kev_list = [
        100,
        200,
        300,
        400,
        500,
        600,
        700,
        800,
        900,
        1000,
        1250,
        1500,
        1750,
        2000,
        2500,
        3000,
        4000,
        5000,
    ]

    ions_He_list = create_ion_list("He", energy_kev_list, units="kev")

    # create target from list of layers
    layer_UO2_10um = Layer(
        {
            "U": {"stoich": 1.0, **elem_u_dict},
            "O": {"stoich": 2.0, **elem_o_dict.as_dict()},
        },
        density=10.97,
        width=100_000.0,
        name="urania",
    )

    layer_SiO2_10um = Layer(
        {
            "Si": {"stoich": 1.0, **elem_si_dict},
            "O": {"stoich": 2.0, **elem_o_dict.as_dict()},
        },
        density=2.65,
        width=100_000.0,
        name="silica",
    )

    target = Target([layer_UO2_10um])

    data_list = pool_srim(
        ions_He_list,
        target,
        data_path,
        num_ions=1_000_000,
        srim_dir=srim_executable_directory,
    )

    pickle_srim(data_list)
