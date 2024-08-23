Matplotlib Gallery
=====================

Orbit of Sun 
--------------

.. code-block:: python


    from mw_plot import MWFaceOn

    from galpy.potential import MWPotential2014
    from galpy.orbit import Orbit
    import numpy as np
    from astropy import units as u

    # Orbit Integration using galpy for the Sun
    op = Orbit([0., 0., 0., 0., 0., 0.], radec=True, ro=8., vo=220.)
    ts = np.linspace(0, 5, 10000) * u.Gyr
    op.integrate(ts, MWPotential2014)
    x = op.x(ts) * u.kpc
    y = op.y(ts) * u.kpc
    z = op.z(ts)

    # setup a mw-plot instance of bird's eye view of the disc
    mw1 = MWFaceOn(radius=20 * u.kpc, unit=u.kpc, coord='galactocentric', annotation=True, figsize=(15, 12), r0=8)

    # set up plot title
    mw1.title = 'Orbit of Sun in 5Gyr'

    # use mw_scatter instead of scatter because we want a colorbar
    mw1.mw_scatter(x, y, c=[z, 'kpc above galactic plane'], s=1)

Orbit of Sun 2 
---------------

.. code-block:: python


    import matplotlib.pyplot as plt
    from mw_plot import MWFaceOn
    from astropy import units as  u
    from galpy.potential import MWPotential2014
    from galpy.orbit import Orbit
    import numpy as np

    # Orbit Integration using galpy for the Sun
    op = Orbit([0., 0., 0., 0., 0., 0.], radec=True, ro=8., vo=220.)
    ts = np.linspace(0, 5, 10000) * u.Gyr
    op.integrate(ts, MWPotential2014)
    x = op.x(ts) * u.kpc
    y = op.y(ts) * u.kpc
    z = op.z(ts)

    # setup a mw-plot instance of bird's eye view of the disc
    mw1 = MWFaceOn(radius=20 * u.kpc, center=(0, 0)*u.kpc, unit=u.kpc, coord='galactocentric', annotation=False, grayscale=True)
    mw2 = MWFaceOn(radius=10 * u.kpc, mode="edge-on", center=(0, 0)*u.kpc, unit=u.kpc, coord='galactocentric', annotation=False, grayscale=True)

    # setup subplots with matplotlib
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5))

    # transform the whole figure with mw-plot
    # mw1.transform([ax1, ax2]) will have the same effect
    mw1.transform(ax1)
    mw2.transform(ax2)

    # you can plot something on top of the transformed subplot
    ax1.scatter(x, y, c='r', s=0.1)
    ax2.scatter(x, z, c='r', s=0.1)
    ax1.set_title("Orbit of the Sun in XY plane", fontsize=20)
    ax2.set_title("Orbit of the Sun in XZ plane", fontsize=20)
