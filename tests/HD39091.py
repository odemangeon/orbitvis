import sys

if ".." not in sys.path:
    sys.path.append("..")
import orbitvis.keplerian
import importlib
importlib.reload(orbitvis.keplerian)
from orbitvis.keplerian import plot_keplerian_orbit
import astropy.constants as cst
import astropy.constants as unt
from math import pi

gm = cst.GM_sun.value  # m3 / s2
mjup2sun = unt.M_jup.to(unt.M_sun).value
mearth2jup = unt.M_earth.to(unt.M_jup).value
au = cst.au.value


def geta(P, Ms, Mp):
    """Return semi-major axis in meters (SI) using kepler equation

    :param float/np.ndarray P: Planetary orbital period in days
    :param float/np.ndarray Ms: Stellar mass in solar mass
    :param float/np.ndarray Mp: Planetary mass in jupiter mass
    :return float/np.ndarray a: Planetary orbital semi-major axis in au
    """
    Ps = P * 24. * 3600.0  # change P to seconds for SI
    return ((Ps / (2. * pi))**2. * gm * (Ms + Mp * mjup2sun))**(1. / 3.) / au


per_c = 6.2682
a_c = geta(P=per_c, Ms=1.094, Mp=4.82 * mearth2jup)
inc_c = 87.27
per_b = 2093.07
a_b = geta(P=per_b, Ms=1.094, Mp=10.02)
delta_inc_b = 2.5
inc_b = inc_c + delta_inc_b
print("inc_b = {} deg".format(inc_b))
fig = plot_keplerian_orbit(a=a_b, per=per_b, ecc=0.637, t_p=2445852.0, OMEGA=180., omega=-330.61, inc=inc_b, fig=None, show=False)

plot_keplerian_orbit(a=a_c, per=per_c, ecc=0., t_p=2458325.5034, OMEGA=180., omega=-90, inc=87.27, fig=fig, show=True)
