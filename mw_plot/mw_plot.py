import os

import numpy as np
import pylab as plt
from astropy import units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable


class MWPlot:
    """
    NAME: MWPlot
    PURPOSE:
    INPUT:
    OUTPUT:
    HISTORY:
        2018-Mar-17 - Written - Henry Leung (University of Toronto)
    """

    def __init__(self, center=(0,0)*u.kpc, radius=90750*u.lyr, unit=u.kpc, coord='galactic', annotation=True, rot180=False):
        """
        :param center: Coordinates of the center of the plot with astropy units
        :type center: u.quantity.Quantity
        :param radius: Radius of the plot with astropy units
        :type radius: u.quantity.Quantity
        :param unit: astropy units
        :type unit: u.quantity.Quantity
        :param coord: 'galactocentric' or 'galactic'
        :type coord: str
        :param annotation: whether use a milkyway background with annotation
        :type annotation: bool
        :param rot180: whether rotate the image by 180 deg
        :type rot180: bool
        """
        self.fontsize = 25
        self.s = 1.0
        self.figsize = (20, 20)
        self.dpi = 200
        self.cmap = "viridis"
        self.imalpha = 0.85
        self.tight_layout = True

        # user should not change these values anyway
        self.__center = center
        self.__radius = radius
        self._unit = unit
        self.__coord = coord
        self.__annotation = annotation
        self.__rot180 = rot180

        self._unit_english = None
        self._coord_english = None
        self.__ext = None
        self.__img = None
        self.__aspect = None

        # Fixed value
        self.__pixels = 5600
        self.__resolution = 24.2 * u.lyr
        self.__fig = None

        # preprocessing procedure
        self._unit_english = self._unit.long_names[0]
        if self.__center.unit is not None and self.__radius.unit is not None:
            self.__center = self.__center.to(self._unit)
            self.__radius = self.__radius.to(self._unit)

        self.images_read()

    def plot(self, *args, **kwargs):
        plt.plot(*args, **kwargs)

    @staticmethod
    def show(*args, **kwargs):
        plt.show(*args, **kwargs)

    def savefig(self, file='MWPlot.png'):
        if self.tight_layout is True:
            plt.tight_layout()

        # this is a pylab method
        self.__fig.savefig(file)

    def images_read(self):
        image_filename = 'MW_bg_annotate.jpg'
        if self.__annotation is False:
            image_filename = 'MW_bg_unannotate.jpg'

        try:
            img = plt.imread(image_filename)
        except FileNotFoundError:
            import mw_plot
            path = os.path.join(os.path.dirname(mw_plot.__path__[0]), 'mw_plot', image_filename)
            img = plt.imread(path)

        if self.__coord.lower() == 'galactic':
            # shift the coord by 8 to the new coord system
            x_shift = 8. * u.kpc
            self.__center[0] += x_shift
            self._coord_english = 'Galactic Coordinates'
        elif self.__coord.lower() == 'galactocentric':
            x_shift = 0. * u.kpc
            self._coord_english = 'Galactocentric Coordinates'
        else:
            raise ValueError("Unknown coordinates, can only be 'galactic' or 'galactocentric'")

        if not type(self.__center) == u.quantity.Quantity and not type(self.__radius) == u.quantity.Quantity:
            print(f"You did not specify units for center and radius, assuming the unit is {self._unit.long_names[0]}")
            if not type(self.__center) == u.quantity.Quantity:
                self.__center = self.__center * self._unit
            if not type(self.__radius) == u.quantity.Quantity:
                self.__radius = self.__radius * self._unit

        self.__resolution = self.__resolution.to(self._unit)
        self.__center = self.__center.to(self._unit)
        self.__radius = self.__radius.to(self._unit)

        # convert physical unit to pixel unit
        pixel_radius = int((self.__radius / self.__resolution).value)
        pixel_center = [int((self.__pixels / 2 + self.__center[0] / self.__resolution).value),
                        int((self.__pixels / 2 - self.__center[1] / self.__resolution).value)]

        # get the pixel coordinates
        x_left_px = pixel_center[0] - pixel_radius
        x_right_px = pixel_center[0] + pixel_radius
        y_top_px = self.__pixels - pixel_center[1] - pixel_radius
        y_bottom_px = self.__pixels - pixel_center[1] + pixel_radius

        # decide whether it needs to fill black pixels because the range outside the pre-compiled images
        if np.all(np.array([x_left_px, self.__pixels - x_right_px, y_top_px, self.__pixels - y_bottom_px]) >= 0):
            img = img[y_top_px:y_bottom_px, x_left_px:x_right_px]
        else:
            # create a black image first with 3 channel with the same data type
            black_img = np.zeros((pixel_radius * 2, pixel_radius * 2, 3), dtype=img.dtype)

            # assign them to temp value
            # just in case the area is outside the images, will fill black pixel
            temp_x_left_px = max(x_left_px, 0)
            temp_x_right_px = min(x_right_px, self.__pixels)
            temp_y_top_px = max(y_top_px, 0)
            temp_y_bottom_px = min(y_bottom_px, self.__pixels)

            left_exceed_px = abs(min(x_left_px, 0))
            top_exceed_px = abs(min(y_top_px, 0))
            # Extract available area from pre-compiled first
            img = img[temp_y_top_px:temp_y_bottom_px, temp_x_left_px:temp_x_right_px]

            # fill the black image with the background image
            black_img[top_exceed_px:top_exceed_px + img.shape[0], left_exceed_px:left_exceed_px + img.shape[1], :] = img

            # Set the images as the filled black-background image
            img = np.array(black_img)

        if self.__rot180:
            img = np.rot90(img, 2)
            self.__ext = [(self.__center[0] + self.__radius - x_shift).value, (self.__center[0] - self.__radius - x_shift).value,
                   (self.__center[1] - self.__radius).value, (self.__center[1] + self.__radius).value]
        else:
            self.__ext = [(self.__center[0] - self.__radius - x_shift).value, (self.__center[0] + self.__radius - x_shift).value,
                   (self.__center[1] + self.__radius).value, (self.__center[1] - self.__radius).value]

        self.__img = img
        self.__aspect = img.shape[0] / float(img.shape[1]) * ((self.__ext[1] - self.__ext[0]) / (self.__ext[3] - self.__ext[2]))

        return None

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

        if not type(x) == u.quantity.Quantity or not type(y) == u.quantity.Quantity:
            raise TypeError("Both x and y must carry astropy's unit")
        else:
            if x.unit is not None and y.unit is not None:
                x = x.to(self._unit)
                y = y.to(self._unit)
            else:
                raise TypeError("Both x, y, center and radius must carry astropy's unit")

        # decide whether we need colorbar or not
        if isinstance(c, list):
            color = c[0]
            cbar_label = c[1]
            cbar_flag = True
            if type(color) == u.quantity.Quantity:
                color = color.to(self._unit)
        else:
            color = c

        self.__fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.title(title, fontsize=self.fontsize)
        plt.xlabel(f'{self._coord_english} ({self._unit_english})', fontsize=self.fontsize)
        plt.ylabel(f'{self._coord_english} ({self._unit_english})', fontsize=self.fontsize)
        ax = plt.gca()
        ax.set_aspect(self.__aspect)
        ax.set_facecolor('k')  # have a black color background for image with <1.0 alpha
        plt.scatter(x, y, zorder=1, s=self.s, c=color, cmap=plt.get_cmap(self.cmap))
        ax.imshow(self.__img, zorder=0, extent=self.__ext, alpha=self.imalpha)

        if cbar_flag is True:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(cax=cax)
            cbar.ax.tick_params(labelsize=self.fontsize)
            cbar.set_label(f"{cbar_label}", size=self.fontsize)

        ax.tick_params(labelsize=self.fontsize)

    def mw_density(self, x, y, c, title=None):
        """
        NAME: mw_density
        PURPOSE:
        INPUT:
        OUTPUT:
        HISTORY:
            2018-Mar-17 - Written - Henry Leung (University of Toronto)
        """
        import seaborn as sns
        import matplotlib

        def transparent_cmap(cmap, N=255):
            "Copy colormap and set alpha values"

            mycmap = cmap
            mycmap._init()
            space = np.logspace(0, 100., N + 4)
            space[0] = 0
            mycmap._lut[:, -1] = space
            return mycmap

        cbar_flag = False

        if not type(x) == u.quantity.Quantity or not type(y) == u.quantity.Quantity:
            raise TypeError("Both x and y must carry astropy's unit")
        else:
            if x.unit is not None and y.unit is not None:
                x = x.to(self._unit)
                y = y.to(self._unit)
            else:
                raise TypeError("Both x, y, center and radius must carry astropy's unit")

        # decide whether we need colorbar or not
        if isinstance(c, list):
            color = c[0]
            cbar_label = c[1]
            cbar_flag = True
            if type(color) == u.quantity.Quantity:
                color = color.to(self._unit)
        else:
            color = c

        self.__fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.title(title, fontsize=self.fontsize)
        plt.xlabel(f'{self._coord_english} ({self._unit_english})', fontsize=self.fontsize)
        plt.ylabel(f'{self._coord_english} ({self._unit_english})', fontsize=self.fontsize)
        ax = plt.gca()
        ax.set_aspect(self.__aspect)
        ax.set_facecolor('k')  # have a black color background for image with <1.0 alpha
        # sns.kdeplot(x.value, y.value, gridsize=1000)
        print([self.__ext[:2], self.__ext[3], self.__ext[2]])
        heatmap, xedges, yedges = np.histogram2d(x.value, y.value, bins=250, range=[self.__ext[:2], [self.__ext[3],
                                                                                    self.__ext[2]]])
        ax.imshow(self.__img, zorder=0, extent=self.__ext, alpha=self.imalpha)
        ax.imshow(heatmap.T, extent=self.__ext, cmap=transparent_cmap(plt.get_cmap('Reds')))

        if cbar_flag is True:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(cax=cax)
            cbar.ax.tick_params(labelsize=self.fontsize)
            cbar.set_label(f"{cbar_label}", size=self.fontsize)

        ax.tick_params(labelsize=self.fontsize)
