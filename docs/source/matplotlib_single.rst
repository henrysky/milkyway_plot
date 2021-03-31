.. automodule:: mw_plot.mw_plot_matplotlib

Single Plot
=============

Classes API
---------------

.. autoclass:: mw_plot.MWPlot
    :members:

.. autoclass:: mw_plot.MWSkyMap
    :members:


MilkyWay Bird's Eyes 
---------------------



MilkyWay Sky Map
------------------

.. code-block:: python
    :linenos:

    import numpy as np
    from astropy import units as  u
    from mw_plot import MWSkyMap

    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # grayscale: whether to turn the background image to grayscale
    mw1 = MWSkyMap(projection='aitoff', grayscale=False)

    # set up plot title
    mw1.title = 'LMC and SMC in red dots'

    # LMC and SMC coordinates
    lsmc_ra = [78.77, 16.26] * u.degree
    lsmc_dec = [-69.01, -72.42] * u.degree

    mw1.scatter(lsmc_ra, lsmc_dec, c='r', s=200)

.. image:: matplotlib_imgs/single_lmcsmc.jpg
