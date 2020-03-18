from dcfclass import (
    pool_srim,
    plot_damage_multi_from_path,
    pickle_srim,
    SrimResults,)

from dcf import create_ion_list, elem_ce_dict, elem_o_dict, elem_si_dict
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


# TODO
"""
plot graphs individually
plot graphs from parent path
plot graphs with limit on x axis
"""

#res = SrimResults(data_list[0].results)
# res.plot_srim(Path('.'))
# plot_damage_multi_from_path(data_parent)
# damage_depth_array = get_depth_damage_array(loaded_data)

# damage_array = get_damage_array(loaded_data, depth=2000)

#fig = plt.figure()
#plot_damage_multi(loaded_data, fig, data_parent)
# data_array = plot_damage_energy_per_ion(loaded_data, data_path, 'nm')
# plot_damage_energy_total(results: Results, folder: Path, units: precisionLitType='nm') -> Tuple[np.ndarray[float], np.ndarray[float]]:
# sum_of_damage = mung_srim(results: Results, depth: int=0)
# plot_srim(results: Results, image_out_dir: Path) -> Tuple[floatArray, floatArray]:

# damage_depth_array = trunc_depth_damage_array(loaded_data, depth=2000)
# stats = get_damage_stats(loaded_data, depth=2000)
# print(stats)
##
