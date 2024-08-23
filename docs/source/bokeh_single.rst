.. automodule:: mw_plot.mw_plot_bokeh

Interactive Single Plot
==========================

Classes API
-------------

.. autoclass:: mw_plot.MWFaceOnBokeh
    :members:

.. autoclass:: mw_plot.MWSkyMapBokeh
    :members:

Interactive MilkywWay Bird's Eye View
--------------------------------------

.. .. bokeh-plot::
..     :source-position: above
.. 

..     from mw_plot import MWFaceOnBokeh
..     from astropy import units as  u
..     from bokeh.io import output_file, show

..     # setup a mw-plot instance of bird's eye view of the disc
..     mw1 = MWPlotBokeh(radius=5 * u.kpc, center=(0, 0)*u.kpc, unit=u.kpc, coord='galactic', grayscale=False, annotation=True)
..     # you should use mw1.show(), I do it this way for the docs to compile correctly
..     show(mw1.bokeh_fig)

Interactive MilkyWay Sky Map
----------------------------------

.. .. bokeh-plot::
..     :source-position: above
.. 

..     from mw_plot import MWSkyMapBokeh
..     from astropy import units as  u
..     from bokeh.io import output_file, show

..     # setup a mw-plot instance of the sky
..     mw1 = MWSkyMapBokeh()
..     # you should use mw1.show(), I do it this way for the docs to compile correctly
..     show(mw1.bokeh_fig)
