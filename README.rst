.. image:: https://raw.githubusercontent.com/henrysky/milkyway_plot/master/mw-plot-logo-b.png
   :height: 132 px
   :width: 500 px

|

|docs||license||ci||coverage|

.. |docs| image:: https://readthedocs.org/projects/milkyway-plot/badge/?version=latest
   :target: https://milkyway-plot.readthedocs.io/en/latest/

.. |license| image:: https://img.shields.io/github/license/henrysky/milkyway_plot.svg
   :alt: GitHub license
   :target: https://github.com/henrysky/milkyway_plot/blob/master/LICENSE

.. |ci| image:: https://github.com/henrysky/astroNN/workflows/CI/badge.svg
   :alt: Build Status
   :target: https://github.com/henrysky/astroNN/actions

.. |coverage| image:: https://codecov.io/gh/henrysky/milkyway_plot/graph/badge.svg?token=NqMxk1j3PQ 
   :target: https://codecov.io/gh/henrysky/milkyway_plot

Getting Started
=================

mw-plot handy python package to do plotting on a face-on/edge-on milkyway/skymap with `matplotlib` (https://matplotlib.org/) for 
static plots or `bokeh` (https://bokeh.org/) for interactive plots (`bokeh` module under develop).

You can set the center and radius of the plot anywhere on a milkyway galaxy image with galactic or galactocentric coordinates.

``mw_plot`` will fill pixels for region outside the pre-compiled images. No acknowledgement to ``mw_plot``
is required if you generate plots for your non-commerical publication, but you **must acknowledgement the origin of
the background images**. The relevant reference text can by obtained by the property ``reference`` of a ``mw_plot`` class instance.

Documentation is available at https://milkyway-plot.readthedocs.io/

Author
---------------

-  | **Henry Leung** - henrysky_
   | Department of Astronomy & Astrophysics, University of Toronto
   | Contact Henry: henrysky.leung [at] utoronto.ca

.. _henrysky: https://github.com/henrysky

License
---------------------------------------------------------

This project is licensed under the MIT License - see the `LICENSE`_ file for details

.. _LICENSE: LICENSE
