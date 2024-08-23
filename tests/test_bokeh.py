from astropy import units as u

from mw_plot.bokeh_backend import MWFaceOnBokeh, MWSkyMapBokeh


def test_mw_skymap_bokeh(simbad):
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

    # LMC and SMC coordinates, get coordinates from Simbad
    result = simbad.query_objects(["LMC", "SMC"])

    # use mw_scatter instead of scatter
    plot_instance.scatter(u.Quantity(result["RA_d"]), u.Quantity(result["DEC_d"]))


def test_mw_plot_bokeh():
    # setup a MWPlotBokeh instance
    plot_instance = MWFaceOnBokeh(
        radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True
    )
    plot_instance.imalpha = 1.0
    plot_instance.s = 10  # make the scatter points bigger

    # set up plot title
    plot_instance.title = "Testing Testing"

    # use mw_scatter instead of scatter because we want a colorbar
    plot_instance.scatter([1, 2, 3] * u.kpc, [1, 2, 3] * u.kpc)
