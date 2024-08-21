.. automodule:: mw_plot.mw_plot_matplotlib

Sky Map of Milky Way
==============================

API
---------------

.. autoclass:: mw_plot.MWSkyMap
    :members:

MilkyWay Sky Map
------------------

By default, ``mw_plot`` has a few background images included within the package which represent ``optical``, ``gamma``, ``far-infrared`` and ``infrared``. 
You can also use other background images from Hierarchical Progressive Surveys (HiPS). This part of the code has made use of the hips2fits, a tool developed at CDS, Strasbourg, France aiming at extracting
  FITS images from HiPS sky maps with respect to a WCS.

You can search for HiPS images with keywords.

.. code-block:: python
    :linenos:

    from mw_plot import MWSkyMap

    # search for HiPS images with keywords
    MWSkyMap.search_sky_background(keywords=None)

which will return a list of all available HiPS images with ``keywords=None``.

You can search for HiPS images with specific keywords. For example

.. code-block:: python
    :linenos:

    from mw_plot import MWSkyMap

    # search for HiPS images with keywords
    MWSkyMap.search_sky_background(keywords="Gaia DR3")

will return a list of all available HiPS images came from Gaia DR3.

.. code-block:: python
    :linenos:

    import numpy as np
    from astropy import units as u
    from mw_plot import MWSkyMap

    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # grayscale: whether to turn the background image to grayscale
    mw1 = MWSkyMap(projection="aitoff", grayscale=False)

    # set up plot title
    mw1.title = "LMC and SMC in red dots"

    # LMC and SMC coordinates
    lsmc_ra = [78.77, 16.26] * u.degree
    lsmc_dec = [-69.01, -72.42] * u.degree

    mw1.scatter(lsmc_ra, lsmc_dec, c="r", s=200)

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from astropy import units as u
    from mw_plot import MWSkyMap

    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # grayscale: whether to turn the background image to grayscale
    mw1 = MWSkyMap(projection="aitoff", grayscale=False)

    # set up plot title
    mw1.title = "LMC and SMC in red dots"

    # LMC and SMC coordinates
    lsmc_ra = [78.77, 16.26] * u.degree
    lsmc_dec = [-69.01, -72.42] * u.degree

    mw1.scatter(lsmc_ra, lsmc_dec, c="r", s=200)
    plt.tight_layout()

You can also plot in different background
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
    :linenos:

    import numpy as np
    from astropy import units as u
    from mw_plot import MWSkyMap

    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # grid: whether to show the Galactic grid
    mw1 = MWSkyMap(projection="aitoff", background="gamma")

    # set up plot title
    mw1.title = "LMC and SMC in Green dots"

    # LMC and SMC coordinates
    lsmc_ra = [78.77, 16.26] * u.degree
    lsmc_dec = [-69.01, -72.42] * u.degree

    mw1.scatter(lsmc_ra, lsmc_dec, c="g", s=200)

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from astropy import units as u
    from mw_plot import MWSkyMap

    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # grid: whether to show the Galactic grid
    mw1 = MWSkyMap(projection="aitoff", background="gamma")

    # set up plot title
    mw1.title = "LMC and SMC in red dots with Galactic Grid"

    # LMC and SMC coordinates
    lsmc_ra = [78.77, 16.26] * u.degree
    lsmc_dec = [-69.01, -72.42] * u.degree

    mw1.scatter(lsmc_ra, lsmc_dec, c="g", s=200)
    plt.tight_layout()

You can also plot with grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
    :linenos:

    import numpy as np
    from astropy import units as u
    from mw_plot import MWSkyMap

    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # grid: whether to show the Galactic grid
    mw1 = MWSkyMap(projection="aitoff", grid="galactic")

    # set up plot title
    mw1.title = "LMC and SMC in red dots with Galactic Grid"

    # LMC and SMC coordinates
    lsmc_ra = [78.77, 16.26] * u.degree
    lsmc_dec = [-69.01, -72.42] * u.degree

    mw1.scatter(lsmc_ra, lsmc_dec, c="r", s=200)

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from astropy import units as u
    from mw_plot import MWSkyMap

    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # grid: whether to show the Galactic grid
    mw1 = MWSkyMap(projection="aitoff", grid="galactic")

    # set up plot title
    mw1.title = "LMC and SMC in red dots with Galactic Grid"

    # LMC and SMC coordinates
    lsmc_ra = [78.77, 16.26] * u.degree
    lsmc_dec = [-69.01, -72.42] * u.degree

    mw1.scatter(lsmc_ra, lsmc_dec, c="r", s=200)
    plt.tight_layout()

.. code-block:: python
    :linenos:

    import numpy as np
    from astropy import units as u
    from mw_plot import MWSkyMap

    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # radecgrid: whether to show the RA/DEC grid
    mw1 = MWSkyMap(projection="aitoff", grid="equatorial")

    # set up plot title
    mw1.title = "LMC and SMC in red dots with RA/DEC Grid"

    # LMC and SMC coordinates
    lsmc_ra = [78.77, 16.26] * u.degree
    lsmc_dec = [-69.01, -72.42] * u.degree

    mw1.scatter(lsmc_ra, lsmc_dec, c="r", s=200)

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from astropy import units as u
    from mw_plot import MWSkyMap

    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # radecgrid: whether to show the RA/DEC grid
    mw1 = MWSkyMap(projection="aitoff", grid="equatorial")

    # set up plot title
    mw1.title = "LMC and SMC in red dots with RA/DEC Grid"

    # LMC and SMC coordinates
    lsmc_ra = [78.77, 16.26] * u.degree
    lsmc_dec = [-69.01, -72.42] * u.degree

    mw1.scatter(lsmc_ra, lsmc_dec, c="r", s=200)
    plt.tight_layout()


.. code-block:: python
    :linenos:

    import numpy as np
    from astropy import units as u
    from mw_plot import MWSkyMap

    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # radecgrid: whether to show the RA/DEC grid
    mw1 = MWSkyMap(projection="aitoff", grid="ecliptic")

    # set up plot title
    mw1.title = "LMC and SMC in red dots with Ecliptic Grid"

    # LMC and SMC coordinates
    lsmc_ra = [78.77, 16.26] * u.degree
    lsmc_dec = [-69.01, -72.42] * u.degree

    mw1.scatter(lsmc_ra, lsmc_dec, c="r", s=200)

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from astropy import units as u
    from mw_plot import MWSkyMap

    # setup a MWSkyMap instance with projection, other projection can be 'hammer', 'mollweide' etc
    # radecgrid: whether to show the RA/DEC grid
    mw1 = MWSkyMap(projection="aitoff", grid="ecliptic")

    # set up plot title
    mw1.title = "LMC and SMC in red dots with Ecliptic Grid"

    # LMC and SMC coordinates
    lsmc_ra = [78.77, 16.26] * u.degree
    lsmc_dec = [-69.01, -72.42] * u.degree

    mw1.scatter(lsmc_ra, lsmc_dec, c="r", s=200)
    plt.tight_layout()