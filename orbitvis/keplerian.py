import matplotlib.pyplot as pl
import astropy.units as unt

from PyAstronomy.pyasl import KeplerEllipse
from matplotlib.collections import PatchCollection
from numpy import linspace, array


def plot_keplerian_orbit(a, per, ecc=0, t_p=0, OMEGA=180, omega=90, inc=90, R_main=1, R_sec=1, fig=None, show=True):
    """Plot a keplerian orbit.

    Show the orbit in the plan of the three planes: xy (plane of the sky), xz and yzselfself.
    Also show a zoom on the transit and eclipse regions.

    :param float      a: Semi-major axis [astropy Quantity/au]
    :param float    per: Orbital period [astropy Quantity/day]
    :param float    ecc: Orbital eccentricity (0-1)
    :param float    t_p: Time of periapsis passage [day]
    :param float  OMEGA: Longitude of the ascending node [astropy Quantity/deg]
    :param float  omega: Argument of periapsis [astropy Quantity/deg]
        Note that the longitude if periapsis is given by OMEGA + omega.
    :param float    inc: Orbit inclination [astropy Quantity/deg].
    :param float R_main: Radius of the main object (for ex. the star) [astropy Quantity/R_sun]
    :param float  R_sec: Radius of the secondary object (for ex. the planet) [astropy Quantity/R_jup]
    :param Figure  fig: Figure
    :param bool   show: If True show the figure before exiting the function
    """
    # Check parameters units
    d_par = {"a": a, "per": per, "ecc": ecc, "t_p": t_p, "OMEGA": OMEGA, "omega": omega, "inc": inc,
             "R_main": R_main, "R_sec": R_sec}
    for param, unit in zip(["a", "per", "t_p", "OMEGA", "omega", "inc", "R_main", "R_sec"],
                           [unt.au, unt.d, unt.d, unt.deg, unt.deg, unt.deg, unt.R_sun, unt.R_jup]):
        if isinstance(d_par[param], unt.Quantity):
            d_par[param] = d_par[param].to(unit).value
    # Main radius in AU
    R_main_au = d_par["R_main"] * unt.R_sun.to(unt.au)
    # Initialise the axes and figure
    if fig is None:
        # Create figure and axes
        fig, axes = pl.subplots(nrows=2, ncols=2)
        # Put labels and titles
        axes[1, 0].set_title("Plan of Sky - full xy")
        axes[1, 0].set_xlabel("x - East ->")
        axes[1, 0].set_ylabel("y - North ->")
        axes[0, 0].set_title("Plan xz")
        axes[1, 1].set_title("Plan zy")
        axes[0, 1].set_title("Plan of Sky - Eclipses Zoom")
        axes[0, 1].set_xlabel("x - East ->")
        axes[0, 1].set_ylabel("y - North ->")
        # Ensure equal scale axis
        for ax in axes.flat:
            ax.axis("equal")
    else:
        axes = array(fig.get_axes()).reshape((2, 2))
    # Create the Keplerian elliptical orbit
    # ke = KeplerEllipse(1, 10, e=0.0, Omega=180., i=90.0, w=90.0)
    ke = KeplerEllipse(d_par["a"], d_par["per"], e=d_par["ecc"], tau=d_par["t_p"], Omega=d_par["OMEGA"],
                       i=d_par["inc"], w=d_par["omega"])
    # Create a time vector to sample the full period
    time = linspace(t_p, t_p + per, 200)
    # Compute the position xyz for the time vector
    pos = ke.xyzPos(time)
    # Find the nodes of the orbit (Observer at -z)
    ascn, descn = ke.xyzNodes_LOSZ()
    # Plot the Main object
    for ax in axes.flat:
        star = pl.Circle((0, 0), radius=R_main_au, facecolor="C8", edgecolor="black")
        ax.add_artist(star)
    # Plot orbit in xy plans
    axes[1, 0].plot(pos[::, 0], pos[::, 1])
    # Plot orbit in xy plans Eclipses zoom
    axes[0, 1].plot(pos[::, 0], pos[::, 1])
    # Do the Zoom
    axes[0, 1].set_xlim((-R_main_au * 1.5, R_main_au * 1.5))
    axes[0, 1].set_ylim((-R_main_au * 1.5, R_main_au * 1.5))
    # Plot orbit in xz plan
    axes[0, 0].plot(pos[::, 0], pos[::, 2])
    # Plot orbit in zy plan
    axes[1, 1].plot(pos[::, 2], pos[::, 1])
    if show:
        pl.show()
    return fig
