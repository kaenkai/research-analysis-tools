"""
@author: Karol Kulinowski
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cst
import json
# import sys
import pandas as pd
# from matplotlib import ticker
# from matplotlib.colors import LogNorm
# from scipy.optimize import curve_fit
# import matplotlib as mpl
# import json

KB = cst.k/cst.eV  # Boltzmann constant [eV]
E0 = cst.e  # electron charge [C]


def dataframe_to_json(df: 'pd.DataFrame', file: str) -> None:
    """Convert dataframe to JSON
    :param df: DataFrame
    :param file: File name to save JSON
    """
    # Generate JSON string
    json_string = df.to_json(orient='index')
    with open(file, 'w') as f:
        f.write(json_string)


def json_to_dict(file: str) -> dict:
    """Convert JSON to dictionary"""
    with open(file, 'r') as f:
        json_string = f.read()
    data_dict = json.loads(json_string)
    return data_dict


def json_to_dataframe(file: str) -> pd.DataFrame:
    """Convert JSON to dataframe
    :param file: File name with JSON data
    """
    data_dict = json_to_dict(file)
    data_df = pd.DataFrame.from_dict(data_dict, orient='index')
    return data_df

def sciNotation(num: float) -> tuple[float, int]:
    """Formats number to SCI notation m*10^n
    :return: number in SCI notation as tuple
    """
    num_sci: list[str] = ('{:E}'.format(num)).split('E')
    return float(num_sci[0]), int(num_sci[1])


def calcDerivative(data_x: 'np.ndarray[float]', data_y: 'np.ndarray[float]', k: int)\
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


def calcLogDerivative(data_x: 'np.ndarray[float]', data_y: 'np.ndarray[float]', k: int)\
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


def calcLogDerivativeFD(x: 'np.ndarray[float]', y: 'np.ndarray[float]')\
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


def effective_mass(x: float) -> float:
    """Calculates effective mass of anatase-rutile phases mix
    :param x: rutile phase weigth content (0:1)
    :return: effective mass of anatase-rutile phases mix
    TODO change values to getters
    """
    return 20*x + 1*(1-x)


def permittivity(x: float, cryst_conn: str) -> float:
    """Calculates permittivity of anatase-rutile phases mix
    :param x: rutile phase weigth content (0:1)
    :param cryst_conn: crystallite connection type (parallel or series)
    :return: permittivity of anatase-rutile phases mix
    TODO change values to getters
    """
    if cryst_conn == 'parallel':
        return x*127 + (1-x)*45
    elif cryst_conn == 'series':
        return 1/(x/127 + (1-x)/45)
    else:
        raise ValueError('cryst_conn has to be either parallel or series')


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
    """Calculates localization lenght
    :param eps: permittivity
    :param tES: characteristic ES temperature [K]
    "return localization lenght [m]
    """
    return 2.8*cst.e**2/(epsilon*cst.epsilon_0*tES*KB*cst.eV)


def fitVRH(x_data, y_data, v=0.25) -> 'np.array[float]':
    """Fit VRH
    :param x_data conductivity:
    :param y_data temperature:
    :param v: VRH conduction model (1/4, 1/3, 1/2 or 1)
    :return: sig0 [Ohm^-1m^-1], t0 [K], sig0_err, t0_err
    """
    x = x_data ** (-v)
    y = np.log(y_data * x_data ** (2. * v))
    fit, cov = np.polyfit(x, y, 1, cov=True)  # linear fit
    t0, sig0 = (-fit[0]) ** (1. / v), np.exp(fit[1])
    sig0_err = np.exp(fit[1]) * np.sqrt(cov[1, 1])
    t0_err = 1 / v * t0 ** (1 - v) * np.sqrt(cov[0, 0])
    return sig0, t0, sig0_err, t0_err


def fitNNH(x_data, y_data) -> 'np.array[float]':
    """Fit NNH
    :param sample: sample id
    :param t_min: min of temperature range
    :param t_max: max of temperature range
    :param cov: calculate covariance matrix?
    :return: sigma_0, E_A, sigma_0 error, E_A error
    """
    x = 1 / x_data / KB
    y = np.log(y_data / th)
    fit, cov = np.polyfit(x, y, 1, cov=True)
    return np.exp(fit[1]), -fit[0], np.exp(fit[1])*np.sqrt(cov[1, 1]), np.sqrt(cov[0, 0])
