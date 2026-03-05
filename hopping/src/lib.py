"""
@author: Karol Kulinowski
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cst
import json
import pandas as pd
import sqlite3

KB = cst.k/cst.eV  # Boltzmann constant [eV]
E0 = cst.e  # electron charge [C]


def dataframe_to_json(df: 'pandas.DataFrame', file: str) -> None:
    """Convert dataframe to JSON
    :param df: DataFrame
    :param file: File name to save JSON
    """
    # Generate JSON string
    json_string = df.to_json(orient='index')
    with open(file, 'w') as f:
        f.write(json_string)


def json_to_dataframe(file: str) -> 'pandas.DataFrame':
    """Convert JSON to dataframe
    :param file: File name with JSON data
    """
    with open(file, 'r') as f:
        json_string = f.read()
    data_df = pd.DataFrame.from_dict(json.loads(json_string), orient='index')
    return data_df


def sciNotation(num: float) -> 'tuple[float, int]':
    """Formats number to SCI notation m*10^n
    :return: number in SCI notation as tuple
    """
    num_sci: list[str] = ('{:E}'.format(num)).split('E')
    return float(num_sci[0]), int(num_sci[1])


def calcDerivative(data_x: 'numpy.ndarray[float]', data_y: 'numpy.ndarray[float]', k: int)\
        -> tuple[np.ndarray[float], np.ndarray[float]]:
    """Calculates derivative
    Ref: A. Möbius, Crit. Rev. Solid State Mater. Sci., 44, 1 (2019).
    :param data_x: X data
    :param data_y: Y data
    :param k: window size
    :return: derivative X and Y as tuple
    """
    assert k % 2 != 0, "Window size must be an odd number."
    n = data_y.size
    arr_xa, der = [], []
    for i in range(k, n):
        x, y = data_x[i - k:i], data_y[i - k:i]
        # print(i, x, y)
        b0 = (k * sum(x * y) - sum(x) * sum(y)) / (k * sum(x ** 2) - sum(x) ** 2)
        b1 = (-k * sum(x ** 3) + sum(x) * sum(x ** 2)) / (k * sum(x ** 2) - sum(x) ** 2)
        xa = -b1 / 2.
        arr_xa.append(xa), der.append(b0)
        # print('{0:5d}{1:10.3f}{2:10.3f}{3:10.3f}'.format(i, b0, np.cos(xa), xa))
    return np.array(arr_xa), np.array(der)


def calcLogDerivative(data_x: 'numpy.ndarray[float]', data_y: 'numpy.ndarray[float]', k: int)\
        -> tuple[np.ndarray[float], np.ndarray[float]]:
    """Calculates logarithmic derivative using Möbius scheme
    Ref: A. Möbius, Crit. Rev. Solid State Mater. Sci., 44, 1 (2019).
    :param data_x: X data
    :param data_y: Y data
    :param k: window size
    :return: logarithmic derivative X and Y as tuple
    """
    x, y = calcDerivative(data_x, np.log(data_y), k)
    y *= x
    return x, y


def calcLogDerivativeFD(x: 'numpy.ndarray[float]', y: 'numpy.ndarray[float]')\
        -> tuple[np.ndarray[float], np.ndarray[float]]:
    """Calculates logarithmic derivative using finite difference differentiation algorithm
    :param x: X data
    :param y: Y data
    :return: derivative X and Y as tuple
    """
    y1 = np.array([(x[i - 1] + x[i + 1]) / 2 * (np.log(y[i + 1]) - np.log(y[i - 1])) / (x[i + 1] - x[i - 1])
                   for i in range(1, x.size - 1)])
    return x[1:-1], y1


def calc_rHopES(tES: float, t: float, xi: float = 1.) -> float:
    """"Calculates average hopping distance in ES regime
    :param tES: Efros-Shklovskii (ES) characteristic temperature [K]
    :param t: temperature [K]
    :param xi: localization length [nm]
    :return: hopping distance in ES regime
    """
    return 1/4. * (tES / t) ** .5 * xi


def calc_dHopES(tES: float, t: float) -> float:
    """"Calculates average hopping energy in ES regime
    :param tES: Efros-Shklovskii (ES) characteristic temperature [K]
    :param t: temperature [K]
    :return: hopping energy in ES regime [eV]
    """
    return 1/2. * cst.k / cst.eV * t * (tES / t) ** .5 / (cst.k / cst.eV * t)


def calc_rHopM(tM: float, t: float, xi: float = 1) -> float:
    """"Calculates average hopping distance in Mott regime
    :param tM: Mott (M) characteristic temperature [K]
    :param t: temperature [K]
    :param xi: localization length [nm]
    :return: hopping distance in Mott regime [nm]
    """
    return 3/8. * (tM / t) ** .25 * xi


def calc_dHopM(tM: float, t: float) -> float:
    """"Calculates average hopping energy in Mott regime
    :param tM: Mott (M) characteristic temperature [K]
    :param t: temperature [K]
    :return: hopping energy in Mott regime [eV]
    """
    return 1/4. * cst.k / cst.eV * t * (tM / t) ** .25


def bohr_radius(m: float, eps: float) -> float: 
    """Calculates Bohr radius
    :param m: effective mass
    :param eps: permittivity
    :return: Bohr radius [m]
    """
    return eps/m * 4*np.pi*cst.epsilon_0*cst.hbar**2/cst.m_e/cst.e**2  # Bohr r. in [m]


def fermi_level_dos(xi: float, tM: float) -> float: 
    """Calculates DOS at Fermi level
    :param xi: localization length (or Bohr radius) [m]
    :param tM: characteristic Mott temperature [K]
    :return: DOS at Fermi level [eV^-1m^-3]
    """
    return 18.1/KB/tM/xi**3


def calc_xi(tES: float, eps:float) -> float:
    """Calculates localization length
    :param eps: permittivity
    :param tES: characteristic ES temperature [K]
    "return localization lenght [m]
    """
    return 2.8*cst.e**2/(epsilon*cst.epsilon_0*tES*KB*cst.eV)


def fitVRH(x_data: float, y_data: float, v: float) -> 'numpy.array[float]':
    """Fit Variable-Range Hopping"""
    x = x_data ** (-v)
    y = np.log(y_data * x_data ** (2. * v))
    fit, cov = np.polyfit(x, y, 1, cov=True)  # linear fit
    t0, sig0 = (-fit[0]) ** (1. / v), np.exp(fit[1])
    sig0_err = np.exp(fit[1]) * np.sqrt(cov[1, 1])
    t0_err = 1 / v * t0 ** (1 - v) * np.sqrt(cov[0, 0])
    return sig0, t0, sig0_err, t0_err


def fitNNH(x_data, y_data) -> 'numpy.array[float]':
    """Fit Nearest-Neighbour Hopping model"""
    x = 1 / x_data / KB
    y = np.log(y_data / th)
    fit, cov = np.polyfit(x, y, 1, cov=True)
    return np.exp(fit[1]), -fit[0], np.exp(fit[1])*np.sqrt(cov[1, 1]), np.sqrt(cov[0, 0])


class ThinFilmDatabase:
    """Class for handling database of TL10 samples"""
    __df = pd.read_sql(
        'SELECT * FROM analysis_dataset',
        sqlite3.connect('src/database/thin_films.db'),
        index_col='sample_id')

    # Plot parameters
    __plot_params = {
        'TL10_5': {'m': '->', 'c': 'C6', 'label': '1'},
        'TL10_7.5': {'m': '-<', 'c': 'C7', 'label': '0.5'},
        'TL10_10': {'m': '-o', 'c': 'C0', 'label': r'$6.6\times10^{-2}$'},
        'TL10_15': {'m': '-s', 'c': 'C1', 'label': r'$7.4\times10^{-5}$'},
        'TL10_17.5': {'m': '-^', 'c': 'C2', 'label': r'$8.2\times10^{-6}$'},
        'TL10_20': {'m': '-v', 'c': 'C3', 'label': r'$7.3\times10^{-7}$'},
        'TL10_25': {'m': '-D', 'c': 'C4', 'label': r'$3.3\times10^{-7}$'},
        'TL10_30': {'m': '-X', 'c': 'C5', 'label': '0'}
    }

    # Temperature ranges
    __ht_range = {
        'TL10_5': [180, 210], 'TL10_7.5': [210, 270],
        'TL10_10': [220, 270], 'TL10_15': [210, 270],
        'TL10_17.5': [170, 270], 'TL10_20': [220, 280],
        'TL10_25': [230, 290], 'TL10_30': [220, 280],
    }
    __lt_range = {
        'TL10_5': [30, 42], 'TL10_7.5': [30, 40],
        'TL10_10': [30, 40], 'TL10_15': [44, 60],
        'TL10_17.5': [100, 120], 'TL10_20': [130, 140],
        'TL10_25': [150, 170], 'TL10_30': [170,180],
    }

    @classmethod
    def load_df(cls) -> 'pandas.DataFrame':
        """Get dataframe"""
        return cls.__df

    @classmethod
    def plt_params(cls, sample) -> dict:
        return cls.__plot_params[sample]

    @classmethod
    def ht_range(cls, sample) -> dict:
        return cls.__ht_range[sample]

    @classmethod
    def lt_range(cls, sample) -> dict:
        return cls.__lt_range[sample]

    @classmethod
    def read_conductivity(cls, sample) -> 'pandas.Series':
        """Read conductivity data in [1/Ohm/cm]"""
        # TODO: filenames from dictionary
        th = cls.__df.at[sample, 'thickness_nm']*1e-7
        df = pd.read_csv(f'src/datafiles/{sample}.csv', 
            comment='#', index_col=0, header=None, names=['temperature', 'conductivity'])/th
        return pd.Series(df['conductivity'], index=df.index)
