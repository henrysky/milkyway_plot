mw_plot
========

A handy python script to plot scatters (so far, will support more later) on a face-on milkyway using pylab.
You can set the center and radius of the plot anywhere on a milkyway galaxy image with galactic or galactocentric coordinates.

Both ``MW_bg_annotate.jpg`` and ``MW_bg_unannotate.jpg`` are modified from an images by **NASA/JPL-Caltech/R. Hurt (SSC/Caltech)**
Both images are 5600x5600px with resolution of 24.2 light years per pixel, mw_plot will fill black pixel for region
outside the pre-compiled images.

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

OR clone the latest commit of mw_plot from github and install

.. code-block:: bash

   $ git clone --depth=1 git://github.com/henrysky/milkyway_plot
   $ python setup.py install

Basic Usage
---------------------

.. code:: python

   from astropy import units as u
   from mw_plot import MWPlot

   # setup a MWPlot instance
   plot_instance = MWPlot()

   # Here are some setting you can set
   plot_instance.fontsize = 25  # fontsize for pylab plotting
   plot_instance.unit = u.kpc  # units of the plot (astropy.units)
   plot_instance.coord = 'galactocentric'  # can be 'galactocentric' or 'galactic'
   plot_instance.center = (0, 0) * u.kpc  # Coordinates of the center of the plot
   plot_instance.radius = 90000 * u.lyr  # Radius of the plot
   plot_instance.figsize = (20, 20)
   plot_instance.dpi = 200
   plot_instance.cmap = 'viridis'  # matplotlib cmap: https://matplotlib.org/examples/color/colormaps_reference.html
   plot_instance.imalpha = 0.85  # alpha value for the milkyway image
   plot_instance.s = 50.0  # make the scatter points bigger
   plot_instance.tight_layout = True # whether plt.tight_layout() will be run
   plot_instance.mw_annotation = True # whether use a milkyway background with annotation

   # Here is the mw_plot if you have an array to color the point
   # x and y must both carry astropy unit
   plot_instance.mw_plot(x, y, [z, 'colorbar_title'], 'Title of the plot here')

   # Here is the mw_plot if you do not have array to color the point
   # x and y must both carry astropy unit
   plot_instance.mw_plot(x, y, 'scatter_point_color_here', 'Title of the plot here')

   # To show
   plot_instance.show()

   # To save
   plot_instance.savefig('name.png')

Example 1: plotting orbit of Sun integrated by galpy
---------------------------------------------------------

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/example_plot_1.png?raw=true

You can plot the orbit which are some scatter points on a face-on milkyway

.. code:: python

   from galpy.potential import MWPotential2014
   from galpy.orbit import Orbit
   import numpy as np
   from astropy import units as u
   from mw_plot import MWPlot

   # Orbit Integration using galpy for the Sun
   op= Orbit([0.,0.,0.,0.,0.,0.], radec=True, ro=8., vo=220.)
   ts = np.linspace(0, 20, 10000)
   op.integrate(ts, MWPotential2014)
   x = - op.x(ts) * u.kpc
   y = op.y(ts) * u.kpc
   z = op.z(ts)

   # setup a MWPlot instance
   plot_instance = MWPlot()
   plot_instance.unit = u.kpc
   plot_instance.coord = 'galactocentric'
   plot_instance.radius = 20 * u.kpc
   plot_instance.imalpha = 1.0
   plot_instance.s = 10

   # plot
   plot_instance.mw_plot(x, y, [z, 'kpc above galactic plane'],
   'Orbit of Sun in 20Gyr using galpy colored by kpc above galactic plane')

   # Save the figure
   plot_instance.savefig(file='mw_plot.png')

   # Show the figure
   plot_instance.show()

You can turn off the annotation by ``plot_instance.mw_annotation = False``

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/example_plot_1_unannotation.png?raw=true

Example 2: Change the center and radius of the plot
---------------------------------------------------------

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/example_plot_2.png?raw=true

You can set the center point and radius of the plot. In this case, we set (16, -3) in galactic coordinates
such that the plot centered at the Sun at the end of 10Gyr orbit, and set the radius as 6 kpc to close up. We will
just set the color to red without color bar title because there is no color bar needed. Please notice the plot assumed
the milkyway is not moving.

.. code:: python

   from galpy.potential import MWPotential2014
   from galpy.orbit import Orbit
   import numpy as np
   from astropy import units as u
   from mw_plot import MWPlot

   # Orbit Integration using galpy for the Sun
   op= Orbit([0.,0.,0.,0.,0.,0.], radec=True, ro=8., vo=220.)
   ts = np.linspace(0, 10, 10000)
   op.integrate(ts, MWPotential2014)
   x = - op.x(ts) * u.kpc
   y = op.y(ts) * u.kpc
   z = op.z(ts)

   # setup a MWPlot instance
   plot_instance = MWPlot()
   plot_instance.unit = u.kpc
   plot_instance.coord = 'galactic'

   # Set the center and radius of the plot
   plot_instance.center = (16, -2.5) * u.kpc
   plot_instance.radius = 6 * u.kpc
   plot_instance.s = 50.0  # make the scatter points bigger

   # plot
   plot_instance.mw_plot(x + 8.*u.kpc, y, 'r', 'Orbit of Sun in 10Gyr using galpy')

   # Save the figure
   plot_instance.savefig(file='mw_plot_zoomed.png')

   # Show the figure
   plot_instance.show()

Example 3: plotting Gaia DR1 observation with astroNN in Galactic coordinates
------------------------------------------------------------------------------------

.. image:: https://github.com/henrysky/milkyway_plot/blob/master/readme_images/example_plot_gaia.png?raw=true

You can set the coord to ``galactic`` to plot observation from Gaia

.. code:: python

   from mw_plot import MWPlot
   from astroNN.gaia import tgas_load
   from astropy import units as  u
   import astropy.coordinates as apycoords

   # Use astroNN to load Gaia TGAS DR1 data files
   # cuts=True to cut bad data (negative parallax and percentage error more than 20%)
   output = tgas_load(dr=1, cuts=True)

   # outout dictionary
   ra = output['ra'] * u.deg  # ra(J2015)
   dec = output['dec'] * u.deg  # dec(J2015)
   parallax = output['parallax']  # parallax
   distance = 1 / parallax * u.kpc

   # error propagation to parsec
   distance_err = (1 / parallax) * output['parallax_err'] / output['parallax'] * 1000

   # use astropy coordinates tranformation
   c = apycoords.SkyCoord(ra=ra, dec=dec, distance=distance, frame='icrs')

   # setup a MWPlot instance
   plot_instance = MWPlot()
   plot_instance.unit = u.kpc
   plot_instance.s = 0.0001
   plot_instance.coord = 'galactic'  # use galactic coordinates because Gaia observations are from Earth

   # Set the center and radius of the plot
   plot_instance.radius = 5 * u.kpc

   plot_instance.s = 50.0  # make the scatter points bigger

   # plot
   plot_instance.mw_plot(c.galactic.cartesian.x, c.galactic.cartesian.y, [distance_err, 'Gaia Distance Error [parsec]'],
                      'Gaia TGAS Distance with 20% error cuts')

   # Save the figure
   plot_instance.savefig(file='gaia.png')

License
---------------------------------------------------------

This project is licensed under the MIT License - see the `LICENSE`_ file for details

.. _LICENSE: LICENSE
