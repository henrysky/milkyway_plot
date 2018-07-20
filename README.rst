mw_plot
========

A handy python package to do plotting on a face-on/edge-on milkyway with matplotlib.
You can set the center and radius of the plot anywhere on a milkyway galaxy image with galactic or galactocentric coordinates.

Both ``MW_bg_annotate.jpg`` and ``MW_bg_unannotate.jpg`` are modified from an images by **NASA/JPL-Caltech/R. Hurt (SSC/Caltech)**
Both images are 5600x5600px with resolution of 24.2 light years per pixel.

``MW_edgeon_unannotate.jpg`` is modified from an images by **ESA/Gaia/DPAC**.
The image is 6500x6500px with resolution of 15.38 light years per pixel taken by ESA Gaia DR2.

mw_plot will fill black pixel for region outside the pre-compiled images.

Author
---------------

-  | **Henry Leung** - *Initial work and developer* - henrysky_
   | Contact Henry: henrysky.leung [at] mail.utoronto.ca

.. _henrysky: https://github.com/henrysky

System Requirement
---------------------

-  | **Python** 3.6 or above
-  | **astropy** 2.0 or above
-  | **Numpy** 1.12.0 or above
-  | **Matplotlib** 2.1.0 above

Install
---------------------

To install via ``pip``

.. code-block:: bash

   $ pip install mw_plot

If something is not working properly, try to upgrade first and then report it as an issue

.. code-block:: bash

   $ pip install mw_plot --upgrade

OR clone the latest commit of mw_plot from github and install

.. code-block:: bash

   $ git clone --depth=1 git://github.com/henrysky/milkyway_plot
   $ python setup.py install

Basic Usage
---------------------

This python package consists of 2 classes - `MWPlot` and `MWSkyMap`. `MWPlot` is used to plot things on a face-on/edge-on milkyway
with galactic/galactocentric coordinates. `MWSkyMap` is used to plot skymap with milkyway background with RA/DEC.

For `MWPlot`:

.. code:: python

   from astropy import units as u
   from mw_plot import MWPlot

   # setup MWPlot instance, you have to specify center, radius, unit with astropy unit and choice of coord
   # or not specifying any to use default value shown below
   # center: Coordinates of the center of the plot
   # radius: Radius of the plot
   # coord: can be 'galactocentric' or 'galactic'
   # annotation: whether use a milkyway background with annotation
   # mode: can be 'face-on' or 'edge-on'

   plot_instance = MWPlot(mode='face-on', center=(0, 0)*u.kpc, radius=90750*u.lyr,
                          unit=u.kpc, coord='galactic', annotation=True, rot180=False)

   # Here are some setting you can set after setting up a MWPlot instance
   plot_instance.title = 'you title here'  # plot title, or it can be None to show no title
   plot_instance.fontsize = 35  # fontsize for matplotlib plotting
   plot_instance.figsize = (20, 20)  # figsize for matplotlib plotting
   plot_instance.dpi = 200  # dpi for matplotlib plotting
   plot_instance.cmap = 'viridis'  # matplotlib cmap: https://matplotlib.org/examples/color/colormaps_reference.html
   plot_instance.clim = (vmin, vmax) # colorbar range
   plot_instance.imalpha = 0.85  # alpha value for the milkyway image
   plot_instance.s = 50.0  # make the scatter points bigger
   plot_instance.tight_layout = True # whether plt.tight_layout() will be run

   # Here is the mw_scatter if you have an array to color the point
   # x and y must both carry astropy unit
   plot_instance.mw_scatter(x, y, [z, 'colorbar_title'])

   # To show
   plot_instance.show()

   # To save
   plot_instance.savefig('name.png')

For `MWSkyMap`:

.. code:: python

   from astropy import units as u
   from mw_plot import MWSkyMap

   # setup MWSkyMap instance, you have to specify grid

   plot_instance = MWSkyMap(grid='galactic')

   # Here are some setting you can set after setting up a MWPlot instance
   plot_instance.title = 'you title here'  # plot title, or it can be None to show no title
   plot_instance.fontsize = 35  # fontsize for matplotlib plotting
   plot_instance.figsize = (20, 20)  # figsize for matplotlib plotting
   plot_instance.dpi = 200  # dpi for matplotlib plotting
   plot_instance.cmap = 'viridis'  # matplotlib cmap: https://matplotlib.org/examples/color/colormaps_reference.html
   plot_instance.clim = (vmin, vmax) # colorbar range
   plot_instance.imalpha = 0.85  # alpha value for the milkyway image
   plot_instance.s = 50.0  # make the scatter points bigger
   plot_instance.tight_layout = True # whether plt.tight_layout() will be run

   # Here is the mw_scatter if you have an array to color the point
   # ra and dec must both carry astropy unit
   plot_instance.mw_scatter(ra, dec, [z, 'colorbar_title'])

   # To show
   plot_instance.show()

   # To save
   plot_instance.savefig('name.png')

There are also some handy constants you can import

.. code:: python

   from mw_plot import center_coord, anti_center_coord
   # center_coord refers to the [RA, DEC] of galactic center in deg
   # anti_center_coord refers to the [RA, DEC] of galactic anti-center in deg

Example 1: Plot Gaia DR1 and DR2 Observation with astroNN in Galactic coordinates
------------------------------------------------------------------------------------

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/example_plot_gaia.png?raw=true

You can set the coord to ``galactic`` to plot observation from Gaia. Please notice if you are using astropy's
coordinates transformation, they will transform under left handed frame, you have to set x = -x to flip it to
right handed which is also the expectation of ``mw_plot``

.. code:: python

    from mw_plot import MWPlot

    from astropy import units as  u
    import astropy.coordinates as apycoords
    import numpy as np

    from astroNN.gaia import gaiadr2_parallax
    from astroNN.gaia import tgas_load

    # To load Gaia DR2 - APOGEE DR14 matches, indices corresponds to APOGEE allstar DR14 file
    ra, dec, parallax, parallax_error = gaiadr2_parallax(cuts=True, keepdims=False)
    distance = 1 / parallax * u.kpc
    ra = ra * u.deg
    dec = dec * u.deg
    distance_err = parallax_error / parallax
    c = apycoords.SkyCoord(ra=ra, dec=dec, distance=distance, frame='icrs')

    # Gaia DR1
    # To load the tgas DR1 files and return a dictionary of ra(J2015), dec(J2015), pmra, pmdec,
    # parallax, parallax error, g-band mag
    # cuts=True to cut bad data (negative parallax and percentage error more than 20%)
    output = tgas_load(cuts=True)
    ra1 = output['ra'] * u.deg  # ra(J2015)
    dec1 = output['dec'] * u.deg  # dec(J2015)
    distance1 = 1 / output['parallax'] * u.kpc
    distance_err1 = output['parallax_err'] / output['parallax']
    c_dr1 = apycoords.SkyCoord(ra=ra1, dec=dec1, distance=distance1, frame='icrs')

    # setup a MWPlot instance
    plot_instance = MWPlot(radius=12 * u.kpc, unit=u.kpc, coord='galactic')

    # so that the colorbar will has a better contract
    plot_instance.clim = (5., 15.)

    # alpha value for the milkyway image
    plot_instance.imalpha = 0.5

    # set up plot title
    plot_instance.title = 'Gaia DR2-APOGEE DR14 matches Distance with 20% error cuts'

    # use mw_scatter instead of scatter because we want a colorbar
    # need to flip the sign of x because astropy is left-handed but mw_plot is right-handed
    plot_instance.mw_scatter(-c.galactic.cartesian.x, c.galactic.cartesian.y,
                             [distance_err * 100, 'Gaia DR2 Distance Precentage Error'])

    # On top of the main plot for DR2, plot DR1 too, need to flip the sign of x because astropy is l
    # eft-handed but mw_plot is right-handed
    plot_instance.scatter(-c_dr1.galactic.cartesian.x, c_dr1.galactic.cartesian.y, c='r',
                          label='Gaia DR1 with 20% distances error cut (Red)')

    # Save the figure
    plot_instance.savefig(file='gaia.png')

Or plotting with Gaia Source with RV catalog (No code is provided but you should be able to modify the code above to do that)

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/gaiadr2_rv_combined.png?raw=true

Example 2: Plot Dynamical Modeling of Tidal Stream using galpy
-----------------------------------------------------------------

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/tidal_streams_plot.png?raw=true

You can plot the orbit which are some scatter points on a edge-on milkyway

.. code:: python

    from mw_plot import MWPlot

    from galpy.df import streamdf
    from galpy.orbit import Orbit
    from galpy.potential import LogarithmicHaloPotential
    from galpy.actionAngle import actionAngleIsochroneApprox
    from galpy.util import bovy_conversion  # for unit conversions
    from astropy import units as u

    # setup potential
    lp = LogarithmicHaloPotential(normalize=1., q=0.9)

    # galpy tidal streams modeling
    aAI = actionAngleIsochroneApprox(pot=lp, b=0.8)
    obs = Orbit([0.16148083, 0.35081535, -0.15481504, 0.48719443, -0.27713334, 0.12019596])
    sigv = 0.365  # km/s
    sdf = streamdf(sigv / 220., progenitor=obs, pot=lp, aA=aAI, leading=True, nTrackChunks=11,
                   tdisrupt=40. / bovy_conversion.time_in_Gyr(220., 8.))

    x = sdf._parse_track_dim('x', interp=True, phys=True)
    y = sdf._parse_track_dim('y', interp=True, phys=True) * u.kpc
    z = sdf._parse_track_dim('z', interp=True, phys=True) * u.kpc

    # setup a MWPlot instance
    plot_instance = MWPlot(mode='edge-on', radius=8. * u.kpc, unit=u.kpc, coord='galactocentric')
    plot_instance.s = 10.  # make the scatter points bigger
    plot_instance.imalpha = 1.0

    # set up plot title
    plot_instance.title = 'Orbit of Sun in 20Gyr using galpy colored by kpc above galactic plane'

    # plot line of the orbit with red color and thicker line
    plot_instance.plot(y, z, c='r', linewidth=4.0)

    # Save the figure
    plot_instance.savefig(file='tidal_streams_plot.png')

Example 3: Plot Orbit of Sun Integrated by galpy
-------------------------------------------------------

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/example_plot_1.png?raw=true

You can plot the orbit which are some scatter points on a face-on milkyway

.. code:: python

    from mw_plot import MWPlot

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

    # setup a MWPlot instance
    plot_instance = MWPlot(radius=20 * u.kpc, unit=u.kpc, coord='galactocentric', annotation=True)
    plot_instance.imalpha = 1.0
    plot_instance.s = 10  # make the scatter points bigger

    # set up plot title
    plot_instance.title = 'Orbit of Sun in 5Gyr using galpy colored by kpc above galactic plane'

    # use mw_scatter instead of scatter because we want a colorbar
    plot_instance.mw_scatter(x, y, [z, 'kpc above galactic plane'])

    # Save the figure
    plot_instance.savefig(file='mw_plot.png')

    # Show the figure
    plot_instance.show()

You can turn off the annotation by putting ``annotation=False`` when creating an instance

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/example_plot_1_unannotation.png?raw=true

Example 4: Change the Center and Radius of the Plot
---------------------------------------------------------

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/example_plot_2.png?raw=true

You can set the center point and radius of the plot. In this case, we set (-16, -2.5) in galactic coordinates
such that the plot centered at the Sun at the end of 10Gyr orbit, and set the radius as 6 kpc to close up. We will
just set the color to red without color bar title because there is no color bar needed. Please notice the plot assumed
the milkyway is not moving.

.. code:: python

    from mw_plot import MWPlot

    from galpy.potential import MWPotential2014
    from galpy.orbit import Orbit
    import numpy as np
    from astropy import units as u

    # Orbit Integration using galpy for the Sun
    op = Orbit([0., 0., 0., 0., 0., 0.], radec=True, ro=8., vo=220.)
    ts = np.linspace(0, 0.5, 10000) * u.Gyr
    op.integrate(ts, MWPotential2014)
    x = op.x(ts) * u.kpc
    y = op.y(ts) * u.kpc
    z = op.z(ts)

    # setup a MWPlot instance with a certain center and radius
    plot_instance = MWPlot(center=(-16, -2.5) * u.kpc, radius=5 * u.kpc)

    # set up plot title
    plot_instance.title = 'Orbit of Sun in 0.5 Gyr using galpy'

    # plot, need to subtract 8kpc to shift to galactic coordinates in right hands frame
    plot_instance.plot(x - 8. * u.kpc, y, c='r', linewidth=8.0)

    # Save the figure
    plot_instance.savefig(file='mw_plot_zoomed.png')

    # Show the figure
    plot_instance.show()

Example 5: Plot all sky map
---------------------------------------------------------

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/adr14_gdr2_skymap.png?raw=true

You can also plot all sky map with mw_plot's MWSkyMap class

.. code:: python

    from mw_plot import MWSkyMap

    import numpy as np
    from astropy import units as  u
    import astropy.coordinates as apycoords
    from astroNN.gaia import gaiadr2_parallax

    ra, dec, parallax, parallax_error = gaiadr2_parallax(cuts=.20, keepdims=False, offset=0.00)

    # setup a MWSkyMap instance
    plot_instance = MWSkyMap(grid='galactic')

    parallax[parallax>1] = 1.

    # so that the colorbar will has a better contract
    # plot_instance.clim = (5., 15.)

    # alpha value for the milkyway image
    plot_instance.imalpha = 1.

    # setup colormap
    plot_instance.cmap='jet'

    # set up plot title
    plot_instance.title = 'APOGEE DR14 coloured by 20% error cuts Gaia Parallax'

    # use mw_scatter instead of scatter because we want a colorbar
    plot_instance.mw_scatter(ra * u.degree, dec * u.degree, [parallax, 'Gaia DR2 Parallax'])

    plot_instance.savefig(file='adr14_gdr2_skymap.png')

    # Show the figure
    plot_instance.show()

Or plotting with Gaia Source with RV catalog for APOGEE DR15 (No code is provided but you should be able to modify the code above to do that)

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/adr15_gdr2_skymap.png?raw=true

License
---------------------------------------------------------

This project is licensed under the MIT License - see the `LICENSE`_ file for details

.. _LICENSE: LICENSE
