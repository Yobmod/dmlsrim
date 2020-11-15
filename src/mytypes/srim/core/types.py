from typing import Union, List, Dict

try:
    from typing import TypedDict  # type: ignore
except ImportError:
    from typing_extensions import TypedDict


elemParamsType = Union[Dict[str, float], List[float], int, float]


class element_TD(TypedDict):
    mass: float
    name: str
    symbol: str
    z: int          # atomic_number


"""
Ac:
  atomic_radii: 195.0
  boil: 3471.0
  color: '#70ABFA'
  covalent_radii: -1.0
  d_elec: 1.0
  density: 10.0699996948
  electronegativity: 1.10000002384
  f_elec: 0.0
  first_ionization_energy: 5.17000007629
  group: 0.0
  mass: 227.0
  melt: 1323.15002441
  name: Actinium
  p_elec: 0.0
  period: 7.0
  production: 0.0010000000475
  s_elec: 2.0
  scattering_factors: {a1: 35.659698, a2: 23.103201, a3: 12.5977, a4: 4.08655, b1: 0.589092,
                       b2: 3.65155, b3: 18.599001, b4: 117.019997, c: 13.5266}
  specific_heat: 0.119999997318
  symbol: Ac
  van_der_waals_radii: -1.0
  volume: 44.8699989319
  z: 89
  """
