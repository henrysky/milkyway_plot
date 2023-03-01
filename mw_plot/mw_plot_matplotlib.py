import os
import warnings
from astropy.coordinates.calculation import HumanError
import numpy as np

import astropy.units as u
import astropy.coordinates as apycoords

import pylab as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mw_plot.mw_plot_masters import MWPlotMaster, MWSkyMapMaster, rgb2gray


__all__ = ["MWPlot", "MWSkyMap"]


class MWPlot(MWPlotMaster):
    """
    MWPlot class plotting with Matplotlib
    
    :param mode: whether plot edge-on or face-on milkyway
    :type mode: str, either 'face-on' or 'edge-on'
    :param center: Coordinates of the center of the plot with astropy units
    :type center: astropy.Quantity
    :param radius: Radius of the plot with astropy units
    :type radius: astropy.Quantity
    :param unit: astropy units
    :type unit: astropy.Quantity
    :param coord: 'galactocentric' or 'galactic'
    :type coord: str
    :param annotation: whether use a milkyway background with annotation
    :type annotation: bool
    :param rot90: number of 90 degree rotation
    :type rot90: int
    :param grayscale: whether to use grayscale background
    :type grayscale: bool
    :param r0: distance to galactic center in kpc
    :type r0: float
    
    :param figsize: Matplotlib figure size
    :type figsize: turple
    :param dpi: Matplotlib figure dpi
    :type dpi: int
    """
    def __init__(self, mode='face-on', center=(0, 0) * u.kpc, radius=90750 * u.lyr, unit=u.kpc, coord='galactic',
                 annotation=True, rot90=0, grayscale=False, r0=8.125, figsize=(7.5, 7.5), dpi=None):
        super().__init__(grayscale=grayscale, 
                         annotation=annotation, 
                         rot90=rot90, 
                         coord=coord, 
                         mode=mode, 
                         r0=r0, 
                         center=center, 
                         radius=radius, 
                         unit=unit, 
                         figsize=figsize, 
                         dpi=dpi)
        self.fontsize = 20
        self.s = 20
        self.cmap = "viridis"
        self.imalpha = 1.
        self.facecolor = 'k' if not grayscale else 'w'
        self.tight_layout = True

        self._unit_english = None
        self._coord_english = None
        self._aspect = None

        self.fig = None
        self.ax = None
        self.title = None
        self.cbar_flag = False
        self.clim = None

        # prepossessing procedure
        self._unit_english = self._unit.short_names[0]
        if self._center.unit is not None and self._radius.unit is not None:
            self._center = self._center.to(self._unit)
            self._radius = self._radius.to(self._unit)

        self.images_read()
        
    def transform(self, x):
        """
        Transform matplotlib figure or a single axes
        """
        if isinstance(x, Figure):
            if len(x.axes) > 1: 
                warnings.warn("More than 1 axes in the figure, mw-plot will populate to all axes")
            fig = x
            ax = fig.axes
        elif isinstance(x, Axes):
            fig = x.figure
            ax = [x]
        elif (isinstance(x, list) or isinstance(x, np.ndarray)) and isinstance(x[0], Axes):
            fig = x[0].figure
            ax = x
        else:
            raise TypeError(f"Your input type {type(x)} is unsupported, can only be matplotlib figure or axes")
        for _ax in ax:
            self.initialize_mwplot(fig, _ax)
    
    def plot(self, x, y, *args, **kwargs):
        x, y = self.xy_unit_check(x, y)
        self.initialize_mwplot()
        self.ax.plot(x, y, zorder=3, *args, **kwargs)
        # just want to set the location right, we dont need image again
        self.ax.imshow(self._img, zorder=0, extent=self._ext, alpha=0., rasterized=True)
        if kwargs.get('label') is not None:
            self.ax.legend(loc='best', fontsize=self.fontsize)

    def scatter(self, x, y, c="r", *args, **kwargs):
        x, y = self.xy_unit_check(x, y)
        self.initialize_mwplot()
        if kwargs.get('s') is None:
            kwargs['s'] = self.s
        self.ax.scatter(x, y, c=c, rasterized=True, *args, **kwargs)
        # just want to set the location right, we dont need image again
        self.ax.imshow(self._img, zorder=0, extent=self._ext, alpha=0., rasterized=True)
        if kwargs.get('label') is not None:
            self.ax.legend(loc='best', fontsize=self.fontsize, markerscale=kwargs['s'])

    def hist2d(self, x, y, *args, **kwargs):
        x, y = self.xy_unit_check(x, y)
        self.initialize_mwplot()
        if kwargs.get('cmap') is None:
            kwargs['cmap'] = self.cmap
        kwargs['cmap'] = self.transparent_cmap(kwargs['cmap'])
        if kwargs.get('range') is None:
            kwargs['range'] = np.array([[self._ext[0], self._ext[1]],[self._ext[2], self._ext[3]]])
        self.ax.hist2d(x, y, zorder=3, *args, **kwargs)
        # just want to set the location right, we dont need image again
        self.ax.imshow(self._img, zorder=0, extent=self._ext, alpha=0., rasterized=True)
        if kwargs.get('label') is not None:
            self.ax.legend(loc='best', fontsize=self.fontsize)

    def show(self, *args, **kwargs):
        if self.fig is None:
            raise AttributeError('Nothing to show, please plot some data first')
        else:
            if self.tight_layout is True:
                self.fig.tight_layout()
            self.fig.show(*args, **kwargs)

    def savefig(self, file='MWPlot.png', dpi="figure", **kwargs):
        if self.tight_layout is True:
            self.fig.tight_layout()
        # this is a pylab method
        self.fig.savefig(file, dpi=dpi, **kwargs)

    @staticmethod
    def transparent_cmap(cmap, N=255):
        """
        Copy colormap and set alpha values

        :param cmap: Color map to covert to transparent color map
        :type cmap: Union[matplotlib.colors.ListedColormap, str]
        :param N: Color map to covert to transparent color map
        :type N: int
        :return: Transparent color map
        :rtype cmap: matplotlib.colors.ListedColormap
        """
        if type(cmap) == str:
            mycmap = plt.get_cmap(cmap)
        else:
            mycmap = cmap
        mycmap._init()
        mycmap._lut[0, -1] = 0
        return mycmap

    def initialize_mwplot(self, fig=None, ax=None):
        """
        Internal method to initial mw_plot images and plot

        :return: None
        """
        if not self._initialized:
            if self.fig is None and fig is None:
                fig, ax = plt.subplots(1, figsize=self.figsize, dpi=self.dpi)
            elif fig is not None:
                pass
            else:
                raise HumanError("Something is wrong duh")
            if self.title is not None:
                ax.set_title(self.title, fontsize=self.fontsize)
            ax.set_xlabel(f'{self._coord_english} ({self._unit_english})', fontsize=self.fontsize)
            ax.set_ylabel(f'{self._coord_english} ({self._unit_english})', fontsize=self.fontsize)
            ax.set_aspect(self._aspect)
            ax.set_facecolor(self.facecolor)  # have a black color background for image with <1.0 alpha
            if not self._grayscale:
                ax.imshow(self._img, extent=self._ext, zorder=0, alpha=self.imalpha, rasterized=True)
            else:
                ax.imshow(self._img[:, :, 0], extent=self._ext, zorder=0, alpha=self.imalpha, rasterized=True, cmap='gray')
            ax.tick_params(labelsize=self.fontsize * 0.8, width=self.fontsize / 10, length=self.fontsize / 2)
            self.fig, self.ax = fig, ax
            
            self._initialized = True

    def mw_scatter(self, x, y, c, **kwargs):
        """
        Plot scatter points with colorbar

        :param x: Scatter points x-coordinates on the plot
        :type x: astropy.Quantity
        :param y: Scatter points y-coordinates on the plot
        :type y: astropy.Quantity
        :param c: Scatter points color
        :type c: Union[str, list, ndarry]
        :History: 2018-Mar-17 - Written - Henry Leung (University of Toronto)
        """
        x, y = self.xy_unit_check(x, y)
        self.initialize_mwplot()

        # decide whether we need colorbar or not
        if isinstance(c, list):
            color = c[0]
            cbar_label = c[1]
            self.cbar_flag = True
            if type(color) == u.quantity.Quantity:
                color = color.to(self._unit).value
        else:
            color = c
            
        if kwargs.get('s') is None:
            kwargs['s'] = self.s

        mappable = self.ax.scatter(x, y, zorder=3, c=color, cmap=plt.get_cmap(self.cmap), rasterized=True,
                                   **kwargs)
        self.ax.imshow(self._img, zorder=0, extent=self._ext, alpha=0., rasterized=True)

        if self.cbar_flag is True:
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = self.fig.colorbar(mappable, cax=cax)
            cbar.ax.tick_params(labelsize=self.fontsize * 0.8, width=self.fontsize / 10, length=self.fontsize / 2)
            cbar.set_label(f"{cbar_label}", size=self.fontsize)
            if self.clim is not None:
                cbar.set_clim(self.clim)

    def mw_density(self, x, y, c, **kwargs):
        """
        Plot desnity with colorbar

        :param x: Scatter points x-coordinates on the plot
        :type x: astropy.Quantity
        :param y: Scatter points y-coordinates on the plot
        :type y: astropy.Quantity
        :param c: Scatter points color
        :type c: Union[str, list, ndarry]
        :param title: Plot title
        :type title: str
        :History: 2018-Mar-17 - Written - Henry Leung (University of Toronto)
        """
        x, y = self.xy_unit_check(x, y)
        self.initialize_mwplot()

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
            self.cbar_flag = True
            if type(color) == u.quantity.Quantity:
                color = color.to(self._unit)
        else:
            color = c

        heatmap, xedges, yedges = np.histogram2d(x.value, y.value, bins=250, range=[self._ext[:2], [self._ext[3],
                                                                                                     self._ext[2]]])
        mappable = self.ax.imshow(heatmap.T, extent=self._ext, cmap=self.transparent_cmap(plt.get_cmap('Reds')),
                                  rasterized=True)
        self.ax.imshow(self._img, zorder=0, extent=self._ext, alpha=0.0, rasterized=True)

        if self.cbar_flag is True:
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = self.fig.colorbar(mappable, cax=cax)
            cbar.ax.tick_params(labelsize=self.fontsize * 0.8, width=self.fontsize / 10, length=self.fontsize / 2)
            cbar.set_label(f"{cbar_label}", size=self.fontsize)
            if self.clim is not None:
                cbar.set_clim(self.clim)

    def scatter_annotate(self, text, position, arrowprops=dict(facecolor="black", width=1., headwidth=6., headlength=6.), fontsize=15, bbox=dict(pad=2), **kwargs):
        """
        Plot annotation with scatter

        :History: 2022-Jan-02 - Written - Henry Leung (University of Toronto)
        """
        if isinstance(position, apycoords.SkyCoord):
            position = self.skycoord_xy(position)
        position_wo_unit = self.xy_unit_check(position[0], position[1])
        position_text = np.add(position_wo_unit, 1.5)
        if isinstance(text, list):
            for t, p, pou, pt in zip(text, position, position_wo_unit, position_text):
                self.scatter(p[0], p[1])
                self.ax.annotate(t, xy=pou, xytext=pt, arrowprops=arrowprops, fontsize=fontsize, bbox=bbox, **kwargs)
        else:
            self.scatter(position[0], position[1])
            self.ax.annotate(text, xy=position_wo_unit, xytext=position_text, arrowprops=arrowprops, fontsize=fontsize, bbox=bbox, **kwargs)

    def annotate(self, *args, **kwargs):
        """
        Plot annotation

        :History: 2022-Jan-02 - Written - Henry Leung (University of Toronto)
        """
        return self.ax.annotate(*args, **kwargs)


class MWSkyMap(MWSkyMapMaster):
    """
    MWSkyMap class plotting with Matplotlib
    
    :param projection: projection system of the plot
    :type projection: string(["equirectangular", "aitoff", "hammer", "lambert", "mollweide"])
    :param center: Coordinates of the center of the plot with astropy degree/radian units
    :type center: astropy.Quantity
    :param radius: Radius of the plot with astropy degree/radian units
    :type radius: astropy.Quantity
    :param grayscale: whether to use grayscale background
    :type grayscale: bool
    
    :param figsize: Matplotlib figure size
    :type figsize: turple
    :param dpi: Matplotlib figure dpi
    :type dpi: int
    """
    def __init__(self, projection='equirectangular', center=(0, 0) * u.deg, radius=(180, 90) * u.deg, grayscale=False, 
                 figsize=(10, 6.5), dpi=None):
        super().__init__(grayscale=grayscale, 
                         projection=projection, 
                         center=center, 
                         radius=radius, 
                         figsize=figsize, 
                         dpi=dpi)
        self._unit = u.degree
        self.fontsize = 20
        self.s = 20.
        self.cmap = "viridis"
        self.imalpha = 1.
        self.tight_layout = True

        self.fig = None
        self.ax = None
        self.title = None
        self.cbar_flag = False
        self.clim = None
        self.facecolor = 'k' if not grayscale else 'w'

        #preprocessing
        if self._projection != 'equirectangular':  # other projections do not support zoom in
            if not np.all(self._center == (0, 0) * u.deg) or not np.all(self._radius == (180, 90) * u.deg):
                print("Projections other than equirectangular does not support custom center and radius, using default!")
                self._center = (0, 0) * u.deg
                self._radius = (180, 90) * u.deg
        else:
            if self._center.unit is not None and self._radius.unit is not None:
                self._center = self._center.to(self._unit)
                self._radius = self._radius.to(self._unit)

        if (self._center[0] + self._radius[0]).value > 180 or (self._center[0] - self._radius[0]).value < -180:
            raise ValueError("The border of the width will be outside the range of -180 to 180 which is not allowed\n")
        if (self._center[1] + self._radius[1]).value > 90 or (self._center[1] - self._radius[1]).value < -90:
            raise ValueError("The border of the height will be outside the range of -90 to 90 which is not allowed")
        if self._radius[0] <= 0 or self._radius[0] <= 0:
            raise ValueError("Radius cannot be negative or 0")

        self.images_read()

    def transform(self, x):
        """
        Transform matplotlib figure or a single axes
        """
        if isinstance(x, Figure):
            if len(x.axes) > 1: 
                warnings.warn("More than 1 axes in the figure, mw-plot will populate to all axes")
            fig = x
            ax = fig.axes
        elif isinstance(x, Axes):
            fig = x.figure
            ax = [x]
        elif (isinstance(x, list) or isinstance(x, np.ndarray)) and isinstance(x[0], Axes):
            fig = x[0].figure
            ax = x
        else:
            raise TypeError(f"Your input type {type(x)} is unsupported, can only be matplotlib figure or axes")
        for _ax in ax:
            self.initialize_mwplot(fig, _ax)

    def initialize_mwplot(self, fig=None, ax=None):
        """
        Initial mw_plot images and plot

        :return: None
        """
        if not self._initialized:
            if self._projection == 'equirectangular':
                if self.fig is None and fig is None:
                    fig, ax = plt.subplots(1, figsize=self.figsize, dpi=self.dpi)
                elif fig is not None:
                    pass
                else:
                    raise HumanError("Something is wrong duh")
                ax.set_xlabel('Galactic Longitude (Degree)', fontsize=self.fontsize)
                ax.set_ylabel('Galactic Latitude (Degree)', fontsize=self.fontsize)
                self._ext = [(self._center[0] - self._radius[0]).value, (self._center[0] + self._radius[0]).value,
                                (self._center[1] - self._radius[1]).value, (self._center[1] + self._radius[1]).value]
                ax.imshow(self._img, zorder=2, extent=self._ext, alpha=self.imalpha, rasterized=True)
            else:  # those cases if there is non-trivial projection
                if self.fig is None and fig is None:
                    fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
                    ax = fig.add_subplot(111, projection=self._projection)
                elif fig is not None:
                    pass
                else:
                    raise HumanError("Something is wrong duh")
                # color
                cmap_red = LinearSegmentedColormap.from_list("", ["black", "red"], N=256)
                cmap_grn = LinearSegmentedColormap.from_list("", ["black", "green"], N=256)
                cmap_blue = LinearSegmentedColormap.from_list("", ["black", "blue"], N=256)

                # coordinates
                lon = np.linspace(-np.pi, np.pi, 6500+1)
                lat = np.linspace(np.pi / 2., -np.pi / 2., 3250+1)
                Lon, Lat = np.meshgrid(lon, lat)
                if self._grayscale:
                    im = ax.pcolormesh(Lon, Lat, np.dot(self._img, [0.2989, 0.5870, 0.1140]), zorder=2, cmap="gray", alpha=self.imalpha, rasterized=True)
                else:
                    rgb = np.array(self._img)
                    color_tuple = rgb.reshape((rgb.shape[0]*rgb.shape[1],rgb.shape[2]))/255.0
                    im = ax.pcolormesh(Lon, Lat, self._img[:, :, 0], zorder=2, color=color_tuple, alpha=self.imalpha, rasterized=True)

            ax.set_facecolor(self.facecolor)  # have a black color background for image with <1.0 alpha
            ax.tick_params(labelsize=self.fontsize * 0.8, width=self.fontsize / 10, length=self.fontsize / 2)
            if self.title is not None:
                ax.set_title(self.title, fontsize=self.fontsize)
            self.fig, self.ax = fig, ax
            
            self._initialized = True

    def show(self, *args, **kwargs):
        if self.fig is None:
            raise AttributeError('Nothing to show, please plot some data first')
        else:
            if self.tight_layout is True:
                self.fig.tight_layout()
            self.fig.show(*args, **kwargs)

    def savefig(self, file='MWSkyMap.png', dpi="figure", **kwargs):
        if self.tight_layout is True:
            self.fig.tight_layout()
        # this is a pylab method
        self.fig.savefig(file, dpi=dpi, **kwargs)

    def mw_scatter(self, ra, dec, c, **kwargs):
        """
        Plot scatter points with colorbar

        :param ra: Scatter points x-coordinates on the plot
        :type ra: astropy.Quantity
        :param dec: Scatter points y-coordinates on the plot
        :type dec: astropy.Quantity
        :param c: Scatter points color
        :type c: Union[str, list, ndarry]
        :History: 2018-Mar-17 - Written - Henry Leung (University of Toronto)
        """
        ra, dec = self.radec_unit_check(ra, dec)
        
        self.initialize_mwplot()

        # decide whether we need colorbar or not
        if isinstance(c, list):
            color = c[0]
            cbar_label = c[1]
            self.cbar_flag = True
            if type(color) == u.quantity.Quantity:
                color = color.to(self._unit).value
        else:
            color = c

        mappable = self.ax.scatter(ra, dec, zorder=3, s=self.s, c=color, cmap=plt.get_cmap(self.cmap), rasterized=True,
                                   **kwargs)
        if self._projection == 'equirectangular':
            self.ax.imshow(self._img, zorder=0, extent=self._ext, alpha=0., rasterized=True)
        else:
            self.ax.imshow(self._img, zorder=0, extent=self._ext, alpha=self.imalpha, rasterized=True,
                           aspect=self.ax.get_aspect(), transform=self.ax.transAxes)
        if self.cbar_flag is True:
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            if self._projection == 'equirectangular':
                cbar = self.fig.colorbar(mappable, cax=cax)
            else:
                cbar = self.fig.colorbar(mappable, ax=self.ax)
            cbar.ax.tick_params(labelsize=self.fontsize * 0.8, width=self.fontsize / 10, length=self.fontsize / 2)
            cbar.set_label(f"{cbar_label}", size=self.fontsize)
            if self.clim is not None:
                cbar.set_clim(self.clim)

    def scatter(self, ra, dec, c="r", *args, **kwargs):
        ra, dec = self.radec_unit_check(ra, dec)
        self.initialize_mwplot()
        if kwargs.get('s') is None:
            kwargs['s'] = self.s
        self.ax.scatter(ra, dec, c=c, zorder=3, rasterized=True, *args, **kwargs)
        # just want to set the location right, we dont need image again
        if self._projection == 'equirectangular':
            self.ax.imshow(self._img, zorder=0, extent=self._ext, alpha=0., rasterized=True)
        else:
            self.ax.imshow(self._img, zorder=0, extent=self._ext, alpha=self.imalpha, rasterized=True,
                           aspect=self.ax.get_aspect(), transform=self.ax.transAxes)
        if kwargs.get('label') is not None:
            self.ax.legend(loc='best', fontsize=self.fontsize, markerscale=kwargs['s'])

    def scatter_annotate(self, text, position, arrowprops=dict(facecolor="black", width=1., headwidth=6., headlength=6.), fontsize=15, bbox=dict(pad=2), **kwargs):
        """
        Plot annotation with scatter

        :History: 2022-Jan-02 - Written - Henry Leung (University of Toronto)
        """
        if isinstance(position, apycoords.SkyCoord):
            position = self.skycoord_radec(position)
        position_wo_unit = self.xy_unit_check(position[0], position[1], checkrot=False)
        position_text = np.add(position_wo_unit, 10)
        if isinstance(text, list):
            for t, p, pou, pt in zip(text, position, position_wo_unit, position_text):
                self.scatter(p[0], p[1])
                self.ax.annotate(t, xy=pou, xytext=pt, arrowprops=arrowprops, fontsize=fontsize, bbox=bbox, **kwargs)
        else:
            self.scatter(position[0], position[1])
            self.ax.annotate(text, xy=position_wo_unit, xytext=position_text, arrowprops=arrowprops, fontsize=fontsize, bbox=bbox, **kwargs)

    def annotate(self, *args, **kwargs):
        """
        Plot annotation

        :History: 2022-Jan-02 - Written - Henry Leung (University of Toronto)
        """
        return self.ax.annotate(*args, **kwargs)
