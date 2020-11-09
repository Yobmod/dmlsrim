from dcfclass import (
    pool_srim,
    pickle_srim,
    SrimResults,)

from dcf import (create_ion_list,
                 elem_ce_dict, elem_o_dict, elem_si_dict,
                 plot_damage_multi_from_path, )
from srim import Layer, Target  # , output
# from srim.output import Results
from pathlib import Path
# import matplotlib.pyplot as plt

energy_kev_list = [
    150, 250, 350, 450
]

ions_He_list = create_ion_list("He", energy_kev_list, units="kev")


# create target from list of layers
layer_CO2_2um = Layer(
    {
        "U": {"stoich": 1.0, **elem_ce_dict.as_dict()},
        "O": {"stoich": 2.0, **elem_o_dict.as_dict()},
    },
    density=7.22,
    width=20_000.0,
    name="ceria",
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

ceria_2um_on_silica = Target([layer_CO2_2um, layer_SiO2_10um])


if __name__ == "__main__":

    # poin to srim exec
    srim_exe_dir = Path(r"C:\srim")
    data_parent = Path(R".\data\ceria_pure__2um")

    data_list = pool_srim(
        ions_He_list,
        ceria_2um_on_silica,
        data_parent,
        num_ions=100,
        srim_dir=srim_exe_dir,
    )

    pickle_srim(data_list)

    res = SrimResults(data_list[0].results)
    res.plot_srim(Path('.'))
    plot_damage_multi_from_path(data_parent)
