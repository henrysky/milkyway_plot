import astropy.coordinates as apycoords
import matplotlib.pyplot as plt
import numpy as np
import pytest
from astropy import units as u

from mw_plot import MWFaceOn, MWSkyMap


def test_hips():
    ls_result = MWSkyMap.search_sky_background("gaia")
    assert any("Gaia" in i for i in ls_result)

    ls_result = MWSkyMap.search_sky_background("gaia dr3")
    assert any("Gaia DR3" in i for i in ls_result)


@pytest.mark.parametrize(
    "projection,grayscale,grid,background",
    [
        ("equirectangular", False, "equatorial", "gamma"),
        ("equirectangular", True, "equatorial", "CDS/P/DM/flux-color-Rp-G-Bp/I/355/gaiadr3"),
        ("aitoff", False, "galactic", "optical"),
        ("hammer", True, "ecliptic", "far-infrared"),
        ("mollweide", True, None, "infrared"),
        ("mollweide", False, None, "Mellinger color optical survey")
    ],
)
def test_mw_skymap(simbad, projection, grayscale, grid, background):
    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # grayscale: whether to turn the background image to grayscale
    plot_instance = MWSkyMap(
        projection=projection, grayscale=grayscale, grid=grid, background=background
    )

    # so that the colorbar will has a better contract
    # plot_instance.clim = (5., 15.)

    # alpha value for the milkyway image
    plot_instance.imalpha = 1.0

    # set up plot title
    plot_instance.title = "LMC and SMC in red dots"
    plot_instance.s = 200

    # LMC and SMC coordinates, get coordinates from Simbad
    result = simbad.query_objects(["LMC", "SMC"])

    # use mw_scatter instead of scatter
    plot_instance.mw_scatter(u.Quantity(result["RA_d"]), u.Quantity(result["DEC_d"]), c="r")

    plot_instance.savefig(file="lmc_smc_projection.png")


def test_mw_skymap_objname():
    # without parallax but have distance
    MWSkyMap(center="M31", radius=(5000, 5000) * u.arcsec, background="Gaia DR3 density map")
    # with parallax
    MWSkyMap(center="LMC", radius=(5000, 5000) * u.arcsec, background="Gaia DR3 density map")
    # without parallax and without distance
    MWSkyMap(center="Galactic Center", radius=(5000, 5000) * u.arcsec, background="Gaia DR3 density map")

    with pytest.raises(Exception):
        MWSkyMap(center="Alien Spaceship", radius=(5000, 5000) * u.arcsec, background="Gaia DR3 density map")

    with pytest.raises(Exception):
        MWSkyMap(center="M31", radius=(5000, 5000) * u.arcsec, background="Alien Spaceship")


@pytest.mark.parametrize(
    "projection,grid,background",
    [
        ("equirectangular", "supergalactic", "optical"),
        ("lambert", "galactic", "optical"),
        ("sinusoidal", "ecliptic", "optical"),
        ("polar", None, "optical"),
        ("rectilinear", "supergalactic", "optical"),
        ("mollweide", None, "rainbow"),
    ],
)
def test_skymap_bad_config(projection, grid, background):
    """
    Test bad configuration for MWSkyMap will raise an exception
    """
    with pytest.raises(Exception):
        MWSkyMap(projection=projection, grid=grid, background=background)


def test_mw_plot():
    # setup a MWPlot instance
    plot_instance = MWFaceOn(
        radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True
    )
    plot_instance.imalpha = 1.0
    plot_instance.s = 10  # make the scatter points bigger

    # set up plot title
    plot_instance.title = (
        "Orbit of Sun in 5Gyr using galpy colored by kpc above galactic plane"
    )

    # use mw_scatter instead of scatter because we want a colorbar
    plot_instance.mw_scatter(
        [1, 2, 3] * u.kpc, [1, 2, 3] * u.kpc, c=[[1, 2, 3] * u.pc, "kpc above galactic plane"]
    )

    plot_instance.savefig(file="gaia.png")


def test_mw_one_annotation():
    mw1 = MWFaceOn(radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True)
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


def test_mwdkymap_one_scatter_annotation():
    # plot
    mw1 = MWSkyMap(projection="equirectangular", grayscale=False)
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


def test_mw_one_scatter_annotation():
    mw1 = MWFaceOn(radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True)

    mw1.title = "Annotation"
    mw1.scatter_annotate(
        ["Earth", "Galactic \n Center"], [[8.0, 0.0], [0.0, 0.0]] * u.kpc
    )
    mw1.savefig(file="mwplot_scatter_annotate_1.jpg")


def test_faceon_transform():
    mw1 = MWFaceOn(radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True)
    fig, ax = plt.subplots(figsize=(10, 5))
    # should raise error if not axes or figure
    with pytest.raises(TypeError):
        mw1.transform(np.array([1,2,3]))

    mw1 = MWFaceOn(radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True)
    fig, ax = plt.subplots(figsize=(10, 5))
    # should not raise error
    mw1.transform(ax)

    # transform all subplots with fig
    fig, ax = plt.subplots(2, 2, figsize=(10, 5))
    mw1.transform(fig)

    # transform all subplots with axes
    fig, ax = plt.subplots(2, 2, figsize=(10, 5))
    mw1.transform(ax[0])


def test_skymap_transform():
    mw1 = MWSkyMap(projection="aitoff", grayscale=False)
    fig, ax = plt.subplots(figsize=(10, 5))
    # should raise error because projection is different
    with pytest.raises(TypeError):
        mw1.transform(ax)
    # should raise error if not axes or figure
    with pytest.raises(TypeError):
        mw1.transform(np.array([1,2,3]))

    mw1 = MWSkyMap(projection="equirectangular", grayscale=False)
    fig, ax = plt.subplots(figsize=(10, 5))
    # should not raise error
    mw1.transform(ax)

    # transform all subplots with fig
    fig, ax = plt.subplots(2, 2, figsize=(10, 5))
    mw1.transform(fig)

    # transform all subplots with axes
    fig, ax = plt.subplots(2, 2, figsize=(10, 5))
    mw1.transform(ax[0])
