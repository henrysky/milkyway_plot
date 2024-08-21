.. image:: https://raw.githubusercontent.com/henrysky/milkyway_plot/master/mw-plot-logo-b.png
   :height: 132 px
   :width: 500 px

|

|pypi| |docs| |license| |ci| |coverage|

Getting Started
=================

``mw-plot`` is a handy Python package for plotting face-on and all-sky maps of the Milky Way using ``matplotlib`` (https://matplotlib.org/) for 
static plots and ``bokeh`` (https://bokeh.org/) for interactive plots.

You can set the center and radius of the plot anywhere on a milkyway galaxy image with galactic or galactocentric coordinates.

Detailed documentation is available at https://milkyway-plot.readthedocs.io/

No acknowledgement to ``mw_plot`` is required if you generate plots for your non-commerical publication, but you **must acknowledgement the origin of
the background images**. The relevant reference text can by obtained by the property ``citation`` of a ``mw_plot`` class instance.

For example, to get the citation text for the background image of a face-on Milky Way plot:

.. code-block:: python

    from mw_plot import MWFaceOn
    mw = MWFaceOn()
    print(mw.citation)

or for the background infrared image of an all-sky Milky Way plot:

.. code-block:: python

    from mw_plot import MWSkyMap
    mw = MWSkyMap(background="infrared")
    print(mw.citation)

.. |docs| image:: https://readthedocs.org/projects/milkyway-plot/badge/?version=latest
   :alt: Documentation Status
   :target: https://milkyway-plot.readthedocs.io/en/latest/

.. |license| image:: https://img.shields.io/github/license/henrysky/milkyway_plot
   :alt: GitHub License
   :target: https://github.com/henrysky/milkyway_plot/blob/master/LICENSE

.. |ci| image:: https://img.shields.io/github/actions/workflow/status/henrysky/milkyway_plot/ci_tests.yml
   :alt: GitHub Actions Workflow Status
   :target: https://github.com/henrysky/astroNN/actions

.. |coverage| image:: https://codecov.io/gh/henrysky/milkyway_plot/graph/badge.svg?token=NqMxk1j3PQ
   :alt: Codecov 
   :target: https://codecov.io/gh/henrysky/milkyway_plot

.. |pypi| image:: https://img.shields.io/pypi/v/mw_plot
   :alt: PyPI - Version
   :target: https://pypi.org/project/mw-plot/

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
