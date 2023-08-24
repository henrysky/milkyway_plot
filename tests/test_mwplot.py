import os
import unittest

import numpy.testing as npt

import numpy as np
import pylab as plt
from astropy import units as u
import astropy.coordinates as apycoords
from mw_plot import MWSkyMap, MWPlot, MWSkyMapBokeh, MWPlotBokeh
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014


class MatplotlibTestCases(unittest.TestCase):
    def test_mw_skymap(self):
        # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
        # grayscale: whether to turn the background image to grayscale
        plot_instance = MWSkyMap(projection="aitoff", grayscale=False)

        # so that the colorbar will has a better contract
        # plot_instance.clim = (5., 15.)

        # alpha value for the milkyway image
        plot_instance.imalpha = 1.0

        # set up plot title
        plot_instance.title = "LMC and SMC in red dots"
        plot_instance.s = 200

        # LMC and SMC coordinates, get coordinates with galpy from_name
        lsmc_ra = [Orbit.from_name("LMC").ra(), Orbit.from_name("SMC").ra()] * u.degree
        lsmc_dec = [
            Orbit.from_name("LMC").dec(),
            Orbit.from_name("SMC").dec(),
        ] * u.degree

        # use mw_scatter instead of scatter
        plot_instance.mw_scatter(lsmc_ra, lsmc_dec, c="r")

        plot_instance.savefig(file="lmc_smc_projection.png")

    def test_mw_plot(self):
        # Orbit Integration using galpy for the Sun
        op = Orbit([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], radec=True, ro=8.0, vo=220.0)
        ts = np.linspace(0, 5, 10000) * u.Gyr
        op.integrate(ts, MWPotential2014)
        x = op.x(ts) * u.kpc
        y = op.y(ts) * u.kpc
        z = op.z(ts)

        # setup a MWPlot instance
        plot_instance = MWPlot(
            radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True
        )
        plot_instance.imalpha = 1.0
        plot_instance.s = 10  # make the scatter points bigger

        # set up plot title
        plot_instance.title = (
            "Orbit of Sun in 5Gyr using galpy colored by kpc above galactic plane"
        )

        # use mw_scatter instead of scatter because we want a colorbar
        plot_instance.mw_scatter(x, y, c=[z, "kpc above galactic plane"])

        plot_instance.savefig(file="gaia.png")

    def test_mw_one_annotation(self):
        mw1 = MWPlot(
            radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True
        )
        mw1.title = "Annotation"
        mw1.scatter(8 * u.kpc, 0 * u.kpc, c="r", s=200)
        mw1.ax.annotate(
            "Earth",
            xy=(8.0, 0.0),
            xytext=(10.0, 1.0),
            arrowprops=dict(arrowstyle="->", shrinkA=0.15),
            fontsize=20,
            bbox=dict(pad=2),
        )

        # Save the figure
        mw1.savefig(file="annotate.jpg")

    def test_mwdkymap_one_scatter_annotation(self):
        # plot
        mw1 = MWSkyMap(projection="equirectangular", grayscale=False, dpi=100)
        mw1.title = "Samples"

        # scatter points
        mw1.scatter(
            [23, 80] * u.degree,
            [23, 80] * u.degree,
            c="r",
            s=20,
            facecolor="none",
            edgecolor="r",
        )

        # annotated scatter
        coords = apycoords.SkyCoord(
            l=[128, 173] * u.degree, b=[-1.0, -1.0] * u.degree, frame="galactic"
        )
        names = ["NGC 1", "NGC 2"]
        mw1.scatter_annotate(names, coords, arrowprops=dict(color="C0"))
        mw1.savefig(file="mwskymap_scatter_annotate_1.jpg")

    def test_mw_one_scatter_annotation(self):
        mw1 = MWPlot(
            radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True
        )

        mw1.title = "Annotation"
        mw1.scatter_annotate(
            ["Earth", "Galactic \n Center"], [[8.0, 0.0], [0.0, 0.0]] * u.kpc
        )
        mw1.savefig(file="mwplot_scatter_annotate_1.jpg")

    def test_mw_skymap_impossible_transform(self):
        mw1 = MWSkyMap(projection="aitoff", grayscale=False)
        fig, ax = plt.subplots(figsize=(10, 5))
        # should raise error
        self.assertRaises(TypeError, mw1.transform, ax)

        mw1 = MWSkyMap(projection="equirectangular", grayscale=False)
        fig, ax = plt.subplots(figsize=(10, 5))
        # should not raise error
        mw1.transform(ax)
        

class BokehTestCases(unittest.TestCase):
    def test_mw_skymap_bokeh(self):
        # setup a MWSkyMapBokeh instance
        # grayscale: whether to turn the background image to grayscale
        plot_instance = MWSkyMapBokeh(grayscale=False)

        # so that the colorbar will has a better contract
        # plot_instance.clim = (5., 15.)

        # alpha value for the milkyway image
        plot_instance.imalpha = 1.0

        # set up plot title
        plot_instance.title = "LMC and SMC in red dots"
        plot_instance.s = 200

        # LMC and SMC coordinates, get coordinates with galpy from_name
        lsmc_ra = [Orbit.from_name("LMC").ra(), Orbit.from_name("SMC").ra()] * u.degree
        lsmc_dec = [
            Orbit.from_name("LMC").dec(),
            Orbit.from_name("SMC").dec(),
        ] * u.degree

        # use mw_scatter instead of scatter
        plot_instance.scatter(lsmc_ra, lsmc_dec)

    def test_mw_plot_bokeh(self):
        # Orbit Integration using galpy for the Sun
        op = Orbit([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], radec=True, ro=8.0, vo=220.0)
        ts = np.linspace(0, 5, 10000) * u.Gyr
        op.integrate(ts, MWPotential2014)
        x = op.x(ts) * u.kpc
        y = op.y(ts) * u.kpc
        z = op.z(ts)

        # setup a MWPlotBokeh instance
        plot_instance = MWPlotBokeh(
            radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True
        )
        plot_instance.imalpha = 1.0
        plot_instance.s = 10  # make the scatter points bigger

        # set up plot title
        plot_instance.title = (
            "Orbit of Sun in 5Gyr using galpy colored by kpc above galactic plane"
        )

        # use mw_scatter instead of scatter because we want a colorbar
        plot_instance.scatter(x, y)


if __name__ == "__main__":
    unittest.main()
