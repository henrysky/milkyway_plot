.. automodule:: mw_plot.mw_plot_matplotlib

Face-On View of Milky Way
==============================

Classes API
---------------

.. autoclass:: mw_plot.MWPlot
    :members:

MilkyWay Bird's Eye
---------------------

.. code-block:: python
    :linenos:

    import numpy as np
    from astropy import units as u
    from mw_plot import MWFaceOn

    # setup a mw-plot instance of bird's eye view of the disc
    mw1 = MWFaceOn(
        radius=20 * u.kpc,
        unit=u.kpc,
        coord="galactocentric",
        annotation=True,
        figsize=(10, 8),
    )

    # set up plot title
    mw1.title = "Bird's Eyes View"

    mw1.scatter(8 * u.kpc, 0 * u.kpc, c="r", s=200)

Annotation
^^^^^^^^^^^

.. code-block:: python
    :linenos:

    import numpy as np
    from astropy import units as u
    from mw_plot import MWFaceOn

    mw1 = MWFaceOn(radius=20 * u.kpc, unit=u.kpc, coord="galactocentric", annotation=True, figsize=(10, 8),)

    # set up plot title
    mw1.title = "Annotation"

    mw1.scatter_annotate(["Earth", "Galactic \n Center"], [[8.0, 0.0], [0.0, 0.0]] * u.kpc)
