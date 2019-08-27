""" Module for automating srim calculations"""


from .output import Results, SRResults
from .input import AutoTRIM, TRIMInput
from .config import DEFAULT_SRIM_DIRECTORY

from pathlib import Path
from typing import TYPE_CHECKING, Dict, Union

if TYPE_CHECKING:
    from .core.ion import Ion
    from .core.target import Target
    from .core.layer import Layer


class TRIMSettings(object):
    """ TRIM Settings

    This object can construct all options available when running a TRIM calculation.

    Parameters
    ----------
    description : :obj:`str`, optional
       A name to give calculation. Has no effect on the actual
       calculation.
    reminders : :obj:`str`, optional
       TODO: could not find description. default 0
    autosave : :obj:`int`, optional
       save calculations after every `autosave` steps. default 0 will
       not autosave except at end
    plot_mode : :obj:`int`, optional
       Default 5.
       (0) ion distribution with recoils projected on y-plane
       (1) ion distribution with recoils projected on z-plane
       (2) ion distribution without recoils projected on y-plane
       (3) transverse plot of ions + recoil cascades, yz-plane
       (4) all four (0-3) on one screen
       (5) no graphics (default and at least 5X faster than others)
    plot_xmin : :obj:`float`, optional
       minimum x depth to plot only really matters if ``plot_mode``
       between 0-4. Default 0.0.
    plot_xmax : :obj:`float`, optional
       maximum x depth to plot only really matters if ``plot_mode``
       between 0-4. Default 0.0.
    ranges : :obj:`bool`, optional
       whether include ``RANGES.txt``, ``RANGE_3D.txt`` to output
       files. Default (0) False
    backscattered : :obj:`bool`, optional
       whether include ``BACKSCAT.txt`` to output files. Default (0)
       False
    transmit : :obj:`bool`, optional
       whether include ``TRANSMIT.txt`` to output files. Default (0)
       False
    sputtered : :obj:`bool`, optional
       whether include ``SPUTTER.txt`` to output files. Default (0)
       False
    collisions : :obj:`bool`, optional
       whether include ``COLLISON.txt`` to output files. Yes they did
       mispell collisions. Default (0) False
    exyz : int
       increment in eV to use for ``EXYZ.txt`` file. Default (0)
    angle_ions : :obj:`float`, optional
       angle of incidence of the ion with respect to the target
       surface. Default (0) perpendicular to the target surface along
       x-axis. Values 0 - 89.9.
    bragg_correction : :obj:`float`, optional
       bragg correction to stopping. Default (0) no correction
    random_seed : :obj:`int`, optional
       a random seed to start calculation with. Default random integer
       between 0 and 100,000. Thus all calculations by default are random.
    version : :obj:`int`, optional
       SRIM-2008 or SRIM-2008 so not really much choice. Default (0)

    Notes
    -----
        This class should never explicitely created. Instead set as
        kwargs in :class:`srim.srim.TRIM`
    """

    def __init__(self, **kwargs: Union[str, int, float]) -> None:
        """Initialize settings for a TRIM running"""
        self._settings: Dict[str, Union[str, int, float]] = ...

    def __getattr__(self, attr: str) -> Union[str, int, float]: ...
    # return self._settings[attr]


class TRIM(object):

    def __init__(self, target: Target, ion: Ion, calculation: int = 1, number_ions: int = 1000, **kwargs: Union[str, int, float]) -> None:
        self.settings: TRIMSettings = ...
        self.calculation: int = ...
        self.number_ions: int = ...
        self.target: Target = ...
        self.ion: Ion = ...

    def _write_input_files(self) -> None:
        """ Write necissary TRIM input files for calculation """
        AutoTRIM().write()
        TRIMInput(self).write()

    @staticmethod
    def copy_output_files(src_directory: Union[str, Path], dest_directory: Union[str, Path], check_srim_output: bool = True) -> None: ...

    def run(self, srim_directory: Union[str, Path] = DEFAULT_SRIM_DIRECTORY) -> Results: ...


class SRSettings(object):
    """ SR Settings

    Parameters
    ----------
    energy_min : :obj:`float`, optional
       lowest energy in [eV] to calculation range
    output_type : :obj:`int`, optional
       specify units for output table
       (1) eV/Angstrom
       (2) keV/micron
       (3) MeV/mm
       (4) keV / (ug/cm2)
       (5) MeV / (mg/cm2)
       (6) keV / (mg/cm2)
       (7) eV / (1E15 atoms/cm2)
       (8) L.S.S reduced units
    output_filename : :obj:`str`, optional
       filename to give for SR output from calcualtion
    correction : :obj:`float`, optional
       Bragg rule correction. Usually no correction needed for heavy
       elements. Default 1.0 implies 100% of value (no change). 1.1
       will increase by 10%.

    Notes
    -----
        This class should never explicitely created. Instead set as
        kwargs in :class:`srim.srim.SR`
    """

    def __init__(self, **args: Union[float, int, str]) -> None:
        self._settings: Dict[str, Union[float, int, str]] = ...

    def __getattr__(self, attr: str) -> Union[str, int, float]: ...
    # return self._settings[attr]


class SR(object):
    """ Automate SR Calculations

    Parameters
    ----------
    leyer : :class:`srim.core.layer.Layer`
        constructed layer for SR calculation
    ion : :class:`srim.core.ion.Ion`
        constructed ion for SR calculation
    kwargs :
        See :class:`srim.srim.SRSettings` for available SR
        options. There are a few and none are required. Defaults are
        appropriate for most cases.
    """

    def __init__(self, layer: Layer, ion: Ion, **kwargs: Union[str, int, float]) -> None:
        self.settings: Dict = ...
        self.layer: Layer = ...
        self.ion: Ion = ...

    def _write_input_file(self) -> None: ...

    def run(self, srim_directory: Union[str, Path] = DEFAULT_SRIM_DIRECTORY) -> SRResults: ...
