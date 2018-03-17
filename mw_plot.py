import pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import units as u


class MWPlot():
    """
    NAME: MWPlot
    PURPOSE:
    INPUT:
    OUTPUT:
    HISTORY:
        2018-Mar-17 - Written - Henry Leung (University of Toronto)
    """
    def __init__(self):
        self.fontsize = 25
        self.unit = u.lyr
        self.coord = None
        self.s = 1.0
        self.figsize = (20, 20)
        self.dpi = 200
        self.cmap = "viridis"
        self.center = (0, 0) * u.lyr
        self.radius = 90000 * u.lyr

        # Fixed value
        self.__resolution = 24.2 * u.lyr
        self.__fig = None

    def plot(self, *args, **kwargs):
        plt.plot(*args, **kwargs)

    def show(self, *args, **kwargs):
        plt.show(*args, **kwargs)

    def savefig(self, file='MWPlot.png'):
        plt.tight_layout()
        self.__fig.savefig(file)

    def images_read(self):
        if self.coord == 'galactic':
            img = plt.imread("MW.png")
            coord_english = 'Galactic Coordinates'
        elif self.coord == 'galactocentric':
            img = plt.imread("MW_galactocentric.png")
            coord_english = 'Galactocentric Coordinates'
        else:
            raise ValueError("Unknown coordinates, can only be 'galactic' or 'galactocentric'")

        if not type(self.center) == u.quantity.Quantity and not type(self.radius) == u.quantity.Quantity:
            print(f"You did not specify units for center and radius, assuming the unit is {self.unit.long_names[0]}")
            self.center = self.center * self.unit
            self.radius = self.radius * self.unit
            self.__resolution = self.__resolution.to(self.unit)

        self.center = self.center.to(self.unit)
        self.radius = self.radius.to(self.unit)

        pixel_radius = self.radius / self.__resolution
        pixel_center = (3750 + self.center / self.__resolution)

        img = img[int((pixel_center[0]-pixel_radius).value):int((pixel_center[0]+pixel_radius).value),
              int((pixel_center[1]-pixel_radius).value):int((pixel_center[1]+pixel_radius).value)]

        ext = [(self.center[1] - self.radius).value, (self.center[1] + self.radius).value,
               (self.center[0] - self.radius).value, (self.center[0] + self.radius).value]

        return img, coord_english, ext

    def mw_plot(self, x, y, c, title=None):
        """
        NAME: mw_plot
        PURPOSE:
        INPUT:
        OUTPUT:
        HISTORY:
            2018-Mar-17 - Written - Henry Leung (University of Toronto)
        """
        cbar_flag = False

        unit_english = self.unit.long_names[0]

        if not type(x) == u.quantity.Quantity or not type(y) == u.quantity.Quantity:
            raise TypeError("Both x and y must carry astropy's unit")
        else:
            x = x.to(self.unit)
            y = y.to(self.unit)

        img, coord_english, ext = self.images_read()

        # decide whether we need colorbar or not
        if isinstance(c, list):
            color = c[0]
            cbar_label = c[1]
            cbar_flag = True
        else:
            color = c

        self.__fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.title(title, fontsize=self.fontsize)
        plt.scatter(x, y, zorder=1, s=self.s, c=color, cmap=plt.get_cmap(self.cmap))
        plt.xlabel(f'{coord_english} ({unit_english})', fontsize=self.fontsize)
        plt.ylabel(f'{coord_english} ({unit_english})', fontsize=self.fontsize)
        aspect = img.shape[0] / float(img.shape[1]) * ((ext[1] - ext[0]) / (ext[3] - ext[2]))
        ax = plt.gca()
        ax.set_aspect(aspect)
        ax.imshow(img, zorder=0, extent=ext)

        if cbar_flag is True:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(cax=cax)
            cbar.ax.tick_params(labelsize=self.fontsize)
            cbar.set_label(f"{cbar_label}", size=self.fontsize)

        ax.tick_params(labelsize=self.fontsize)
