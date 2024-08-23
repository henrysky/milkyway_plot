Basic Usage
=================

This python package consists of two major abstract classes - ``MWPlotBase`` and ``MWSkyMapBase`` with Matplotlib or Bokeh backend as well as  a few useful utilities. 

``MWPlotBase`` is used to plot things on a face-on milkyway with galactic and galactocentric coordinates. 
``MWSkyMapBase`` is used to plot skymap with milkyway background with RA/DEC.

Useful constants
-------------------

.. code-block:: python
    :caption: A few usage constants

    >>> from mw_plot import center_radec, anti_center_radec, northpole_radec, southpole_radec  # constants
    >>> from mw_plot import mw_radec # milkyway plane in RA/DEC

    >>> center_radec  # refers to the [RA, DEC] of galactic center in deg
    [266.4167, -29.0078]
    >>> anti_center_radec  # refers to the [RA, DEC] of galactic anti-center in deg
    [86.4167, 28.0078]
    >>> northpole_radec  # refers to the [RA, DEC] of galactic north pole in deg
    [192.7667, 27.1167]
    >>> southpole_radec  # refers to the [RA, DEC] of galactic south pole in deg
    [12.7667, -27.1167]

    >>> mw_plane_ra, mw_plane_dec = mw_radec(deg=True)  # RA/DEC arrays of milkyway plane

.. image:: mpl_imgs/mw_radec_constants.jpg
    :width: 500
    :align: center