.. milkyway_plot documentation master file, created by
   sphinx-quickstart on Sun Mar 21 14:34:36 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to mw-plot's documentation!
=========================================

``mw-plot`` is a handy Python package for plotting face-on and all-sky maps of the Milky Way using ``matplotlib`` (https://matplotlib.org/) for 
static plots and ``bokeh`` (https://bokeh.org/) for interactive plots.

You can set the center and radius of the plot anywhere on a milkyway galaxy image with galactic or galactocentric coordinates.

Detailed documentation is available at https://milkyway-plot.readthedocs.io/

No acknowledgement to ``mw_plot`` is required if you generate plots for your non-commerical publication, but you **must acknowledgement the origin of
the background images**. The relevant reference text can by obtained by the property ``citation`` of a ``mw_plot`` class instance.

For example, to get the citation text for the background image of a face-on Milky Way plot:

.. code-block:: python

    >>> from mw_plot import MWFaceOn
    >>> mw = MWFaceOn()
    >>> mw.citation
    'NASA/JPL-Caltech/R. Hurt (SSC/Caltech)'

or for the background infrared image of an all-sky Milky Way plot:

.. code-block:: python

    >>> from mw_plot import MWSkyMap
    >>> mw = MWSkyMap()
    >>> mw.citation
    'ESA/Gaia/DPAC'

Install
---------------------

To install via ``pip``

.. prompt:: bash $
   
   python -m pip install mw_plot


If something is not working properly, try to upgrade first and then report it as an issue on Github if the problem persists

.. prompt:: bash $

   python -m pip install mw_plot --upgrade


You can also clone the latest commit of mw_plot from Github by running the following command to install the package in editable mode

.. prompt:: bash $

   git clone --depth=1 https://github.com/henrysky/milkyway_plot.git
   python -m pip install -e .

System Requirement
---------------------

``mw-plot`` requires Python 3.9 or above. The following packages are required which will be installed automatically when you install ``mw_plot``:

.. literalinclude:: ../../requirements.txt

.. toctree::
   :maxdepth: 1
   :caption: Getting Started
   
   basics
   changelog

.. toctree::
   :maxdepth: 2
   :caption: Matplotlib backend
   
   matplotlib_faceon
   matplotlib_skymap
   matplotlib_multi

.. toctree::
   :maxdepth: 2
   :caption: Bokeh backend

   bokeh_single
   

Author
---------------

-  | **Henry Leung** - henrysky_
   | Department of Astronomy & Astrophysics, University of Toronto
   | Contact Henry: henrysky.leung [at] utoronto.ca

.. _henrysky: https://github.com/henrysky

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
