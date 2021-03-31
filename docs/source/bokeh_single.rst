.. automodule:: mw_plot.mw_plot_bokeh

Single Bokeh Plot
=====================

Classes API
-------------

.. autoclass:: mw_plot.MWPlotBokeh
    :members:

.. autoclass:: mw_plot.MWSkyMapBokeh
    :members:

Single MWPlotBokeh
----------------------

.. bokeh-plot::
    :source-position: above
    :linenos:

    from mw_plot import MWPlotBokeh
    from astropy import units as  u
    from bokeh.io import output_file, show

    # setup a mw-plot instance of bird's eyes view of the disc
    mw1 = MWPlotBokeh(radius=5 * u.kpc, center=(0, 0)*u.kpc, unit=u.kpc, coord='galactic', grayscale=False, annotation=True)
    # you should use mw1.show(), I do it this way for the docs to compile correctly
    show(mw1.bokeh_fig)

Single MWPlotBokeh
----------------------
