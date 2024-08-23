import warnings
from typing import Tuple, Union

import astropy.coordinates as apycoords
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mw_plot.base import MWPlotBase, MWSkyMapBase


class MWFaceOn(MWPlotBase):
    """
    MWPlot class plotting with Matplotlib

    Parameters
    ----------
    grayscale : bool, optional
        Whether to use grayscale background. The default is False.
    annotation : bool, optional
        Whether to show annotation. The default is False.
    angle : int, optional
        Angle of the plot. The default is 90.
    r0 : astropy.Quantity, optional
        Distance to the Galactic center. The default is 8.125*u.kpc.
    coord : str, optional
        Coordinate system. The default is "galactic".
    center : tuple, optional
        Center of the plot. The default is (0.0, 0.0).
    radius : astropy.Quantity, optional
        Radius of the plot. The default is 20.0*u.kpc.
    unit : astropy.Unit, optional
        Unit of the plot. The default is u.kpc.
    figsize : tuple, optional
        Matplotlib figure size. The default is (5, 5).
    """

    def __init__(
        self,
        grayscale: bool = False,
        annotation: bool = False,
        angle: int = 90,
        r0: u.Quantity = 8.125 * u.kpc,
        coord: str = "galactic",
        center: tuple = (0.0, 0.0),
        radius: u.Quantity = 20.0 * u.kpc,
        unit: u.Unit = u.kpc,
        figsize: tuple = (5, 5),
    ):
        super().__init__(
            grayscale=grayscale,
            annotation=annotation,
            angle=angle,
            r0=r0,
            coord=coord,
            center=center,
            radius=radius,
            unit=unit,
            figsize=figsize,
        )
        self.s = 20
        self.cmap = "viridis"
        self.imalpha = 1.0
        self.tight_layout = True

        self.unit_english = None
        self.coord_english = None
        self._aspect = None

        self.fig = None
        self.ax = None
        self.title = None
        self.cbar_flag = False
        self.clim = None

        # prepossessing procedure
        self.unit_english = self.unit.short_names[0]
        self.unit_check(self.center, self.unit)
        self.unit_check(self.radius, self.unit)

        self.read_bg_img()

    def transform(self, x):
        """
        Transform matplotlib figure or a single axes
        """
        if isinstance(x, Figure):
            if len(x.axes) > 1:
                warnings.warn(
                    "More than 1 axes in the figure, mw-plot will populate to all axes"
                )
            fig = x
            ax = fig.axes
        elif isinstance(x, Axes):
            fig = x.figure
            ax = [x]
        elif (isinstance(x, list) or isinstance(x, np.ndarray)) and isinstance(
            x[0], Axes
        ):
            fig = x[0].figure
            ax = x
        else:
            raise TypeError(
                f"Your input type {type(x)} is unsupported, can only be matplotlib figure or axes"
            )
        for _ax in ax:
            self.initialize_mwplot(fig, _ax, _multi=True)

    def plot(self, x, y, *args, **kwargs):
        x, y = self.xy_unit_check(x, y)
        self.initialize_mwplot()
        self.ax.plot(x, y, zorder=3, *args, **kwargs)
        # just want to set the location right, we dont need image again
        self.ax.imshow(
            self.bg_img, zorder=0, extent=self._ext, alpha=0.0, rasterized=True
        )
        if kwargs.get("label") is not None:
            self.ax.legend(loc="best")

    def scatter(self, x, y, c="r", *args, **kwargs):
        x, y = self.xy_unit_check(x, y)
        self.initialize_mwplot()
        if kwargs.get("s") is None:
            kwargs["s"] = self.s
        self.ax.scatter(x, y, c=c, rasterized=True, *args, **kwargs)
        # just want to set the location right, we dont need image again
        self.ax.imshow(
            self.bg_img, zorder=0, extent=self._ext, alpha=0.0, rasterized=True
        )
        if kwargs.get("label") is not None:
            self.ax.legend(loc="best", markerscale=kwargs["s"])

    def show(self, *args, **kwargs):
        if self.fig is None:
            raise AttributeError("Nothing to show, please plot some data first")
        else:
            if self.tight_layout is True:
                self.fig.tight_layout()
            self.fig.show(*args, **kwargs)

    def savefig(self, file="MWPlot.png", **kwargs):
        if self.tight_layout is True:
            self.fig.tight_layout()
        # this is a pylab method
        self.fig.savefig(file, **kwargs)

    def initialize_mwplot(self, fig=None, ax=None, _multi=False):
        """
        Internal method to initial mw_plot images and plot

        :return: None
        """
        if not self._built or _multi:
            if self.fig is None and fig is None:
                fig, ax = plt.subplots(1, figsize=self.figsize)
            elif fig is not None:
                pass
            else:
                raise NotImplementedError("I think no one will ever reach here")
            if self.title is not None:
                ax.set_title(self.title)
            ax.set_xlabel(
                f"{self.coord_english} ({self.unit_english})"
            )
            ax.set_ylabel(
                f"{self.coord_english} ({self.unit_english})"
            )
            ax.set_aspect(self._aspect)
            ax.set_facecolor(
                self.facecolor
            )  # have a black color background for image with <1.0 alpha
            if not self.grayscale:
                ax.imshow(
                    self.bg_img,
                    extent=self.bg_img_ext,
                    zorder=0,
                    alpha=self.imalpha,
                    rasterized=True,
                )
                ax.set_xlim(self._ext[0], self._ext[1])
                ax.set_ylim(self._ext[2], self._ext[3])
            else:
                ax.imshow(
                    self.bg_img[:, :, 0],
                    extent=self.bg_img_ext,
                    zorder=0,
                    alpha=self.imalpha,
                    rasterized=True,
                    cmap="gray",
                )
                ax.set_xlim(self._ext[0], self._ext[1])
                ax.set_ylim(self._ext[2], self._ext[3])
            self.fig, self.ax = fig, ax

            self._built = True

    def mw_scatter(self, x, y, c="r", **kwargs):
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
            if hasattr(c[0], "__len__"):
                color = c[0]
                cbar_label = c[1]
                self.cbar_flag = True
                if isinstance(color, u.quantity.Quantity):
                    color = color.value
            else:
                color = c
        else:
            color = c

        if kwargs.get("s") is None:
            kwargs["s"] = self.s

        mappable = self.ax.scatter(
            x,
            y,
            zorder=3,
            c=color,
            cmap=plt.get_cmap(self.cmap) if self.cbar_flag else None,
            rasterized=True,
            **kwargs,
        )
        self.ax.imshow(
            self.bg_img, zorder=0, extent=self.bg_img_ext, alpha=0.0, rasterized=True
        )

        if self.cbar_flag is True:
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = self.fig.colorbar(mappable, cax=cax)
            cbar.set_label(f"{cbar_label}")
            if self.clim is not None:
                cbar.set_clim(self.clim)

    def scatter_annotate(
        self,
        text,
        position,
        arrowprops=dict(facecolor="black", width=1.0, headwidth=6.0, headlength=6.0),
        fontsize=15,
        bbox=dict(pad=2),
        **kwargs,
    ):
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
                self.ax.annotate(
                    t,
                    xy=pou,
                    xytext=pt,
                    arrowprops=arrowprops,
                    fontsize=fontsize,
                    bbox=bbox,
                    **kwargs,
                )
        else:
            self.scatter(position[0], position[1])
            self.ax.annotate(
                text,
                xy=position_wo_unit,
                xytext=position_text,
                arrowprops=arrowprops,
                fontsize=fontsize,
                bbox=bbox,
                **kwargs,
            )

    def annotate(self, *args, **kwargs):
        """
        Plot annotation

        :History: 2022-Jan-02 - Written - Henry Leung (University of Toronto)
        """
        return self.ax.annotate(*args, **kwargs)


class MWSkyMap(MWSkyMapBase):
    """
    MWSkyMap class plotting with Matplotlib

    Parameters
    ----------
    grayscale : bool, optional
        Whether to use grayscale background. The default is False.
    projection : str, optional
        Projection of the plot. The default is "equirectangular".
    background : str, optional
        Background image of the plot. The default is "optical".
        You can use ``MWSkyMap.search_background(keyword=None)`` to search for available background images.
    center : Union[Tuple[float, float], str], optional
        Center of the plot. The default is (0.0, 0.0) * u.deg.
    radius : tuple, optional
        Radius of the plot. The default is (180.0, 90.0).
    grid : str, optional
        Grid of the plot. The default is None.
    figsize : Tuple[float, float], optional
        Matplotlib figure size. The default is (5, 5).
    """

    def __init__(
        self,
        grayscale: bool = False,
        projection: str = "equirectangular",
        background: str = "optical",
        center: Union[Tuple[float, float], str] = (0.0, 0.0) * u.deg,
        radius: tuple = (180.0, 90.0) * u.deg,
        grid: str = None,
        figsize: Tuple[float, float] = (6, 4),
    ):
        super().__init__(
            grayscale=grayscale,
            projection=projection,
            background=background,
            center=center,
            radius=radius,
            figsize=figsize,
        )
        self.unit = u.degree
        self.s = 20.0
        self.cmap = "viridis"
        self.imalpha = 1.0
        self.tight_layout = True

        self.grid = False
        self.radecgrid = False
        self.eclgrid = False
        if grid is None:
            self.grid = False
        elif grid == "galactic":
            self.grid = True
        elif grid == "equatorial":
            self.radecgrid = True
        elif grid == "ecliptic":
            self.eclgrid = True
        else:
            raise ValueError(
                "grid must be either 'galactic', 'equatorial' or 'ecliptic'"
            )

        self.fig = None
        self.ax = None
        self.title = None
        self.cbar_flag = False
        self.clim = None
        self.facecolor = "k" if not grayscale else "w"

        # preprocessing
        if (
            self.projection != "equirectangular"
        ):  # other projections do not support zoom in
            if not np.all(self.center == (0, 0) * u.deg) or not np.all(
                self.radius == (180, 90) * u.deg
            ):
                warnings.warn(
                    "Projections other than equirectangular does not support custom center and radius, "
                    "using default center=(0, 0) degree and radius=(180, 90) degree"
                )
                self.center = (0, 0) * u.deg
                self.radius = (180, 90) * u.deg
        else:
            self.unit_check(self.center, self.unit)
            self.unit_check(self.radius, self.unit)

        if (self.center[0] + self.radius[0]).value > 180 or (
            self.center[0] - self.radius[0]
        ).value < -180:
            raise ValueError(
                "The border of the width will be outside the range of -180 to 180 which is not allowed\n"
            )
        if (self.center[1] + self.radius[1]).value > 90 or (
            self.center[1] - self.radius[1]
        ).value < -90:
            raise ValueError(
                "The border of the height will be outside the range of -90 to 90 which is not allowed"
            )
        if self.radius[0] <= 0 or self.radius[0] <= 0:
            raise ValueError("Radius cannot be negative or 0")

        self.read_bg_img()

    def transform(self, x):
        """
        Transform matplotlib figure or a single axes
        """
        if self.projection == "equirectangular":
            projection_name = "rectilinear"
        else:
            projection_name = self.projection
        if isinstance(x, Figure):
            if len(x.axes) > 1:
                warnings.warn(
                    "More than 1 axes in the figure, mw-plot will populate to all axes"
                )
            fig = x
            ax = fig.axes
        elif isinstance(x, Axes):
            fig = x.figure
            ax = [x]
        elif (isinstance(x, list) or isinstance(x, np.ndarray)) and isinstance(
            x[0], Axes
        ):
            fig = x[0].figure
            ax = x
        else:
            raise TypeError(
                f"Your input type {type(x)} is unsupported, can only be matplotlib figure or axes"
            )
        for _ax in ax:
            if _ax.name != projection_name:
                raise TypeError(
                    f"You can not transform a figure with different projection, you want to transform to '{projection_name}' but your figure is '{x.name}'"
                )
            self.initialize_mwplot(fig, _ax, _multi=True)

    def initialize_mwplot(self, fig=None, ax=None, _multi=False):
        """
        Initial mw_plot images and plot

        :return: None
        """
        if self.projection == "equirectangular":
            self._fake_rad2deg = np.rad2deg
        else:
            self._fake_rad2deg = lambda x: x
        if not self._built or _multi:
            if self.projection == "equirectangular":
                if self.fig is None and fig is None:
                    fig, ax = plt.subplots(1, figsize=self.figsize)
                elif fig is not None:
                    pass
                else:
                    raise NotImplementedError("I think no one will ever reach here")
                ax.set_xlabel(r"$l$ (deg)")
                ax.set_ylabel(r"$b$ (deg)")
                self._ext = [
                    (self.center[0] - self.radius[0]).value,
                    (self.center[0] + self.radius[0]).value,
                    (self.center[1] - self.radius[1]).value,
                    (self.center[1] + self.radius[1]).value,
                ]
                ax.imshow(
                    self.bg_img,
                    zorder=2,
                    extent=self._ext,
                    alpha=self.imalpha,
                    rasterized=True,
                )
            else:  # those cases if there is non-trivial projection
                if self.fig is None and fig is None:
                    fig = plt.figure(figsize=self.figsize)
                    ax = fig.add_subplot(111, projection=self.projection)
                elif fig is not None:
                    pass
                else:
                    raise NotImplementedError("I think no one will ever reach here")

                # coordinates
                lon = np.linspace(-np.pi, np.pi, self.bg_img.shape[1] + 1)
                lat = np.linspace(np.pi / 2.0, -np.pi / 2.0, self.bg_img.shape[0] + 1)
                Lon, Lat = np.meshgrid(lon, lat)
                if self.grayscale:
                    mappable = ax.pcolormesh(
                        Lon,
                        Lat,
                        np.dot(self.bg_img, [0.2989, 0.5870, 0.1140]),
                        zorder=2,
                        cmap="gray",
                        alpha=self.imalpha,
                        rasterized=True,
                    )
                else:
                    mappable = ax.pcolormesh(
                        Lon,
                        Lat,
                        self.bg_img,
                        zorder=2,
                        alpha=self.imalpha,
                        rasterized=True,
                    )

            # ax.set_facecolor(
            #     self.facecolor
            # )  # have a black color background for image with <1.0 alpha
            if self.title is not None:
                ax.set_title(self.title, y=1.05)
            self.fig, self.ax = fig, ax

            grad_alpha = 0.5
            grid_width = 0.5
            grid_style = "--"

            self._built = True
            if self.projection != "equirectangular":
                self.ax.set_xticklabels([])
            if self.grid is True:
                for i in [0, -15, 15, -30, 30, -45, 45, -60, 60, -75, 75]:
                    self.ax.plot(
                        self._fake_rad2deg(np.deg2rad([-180, 180])),
                        self._fake_rad2deg(np.deg2rad([i, i])),
                        c=self._opposite_color,
                        lw=grid_width,
                        ls=grid_style,
                        zorder=3,
                        alpha=grad_alpha,
                    )
                for i in [-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150]:
                    self.ax.plot(
                        self._fake_rad2deg(np.deg2rad([i, i])),
                        self._fake_rad2deg(np.deg2rad([-75, 75])),
                        c=self._opposite_color,
                        lw=grid_width,
                        ls=grid_style,
                        zorder=3,
                        alpha=grad_alpha,
                    )
            elif self.projection == "equirectangular":
                pass
            else:
                # disable ticks if not galactic grid
                self.ax.set_yticklabels([])

            epoch = "J2000"

            def radec_to_lb(ra, dec, degree=False):
                if degree is True:
                    unit = u.deg
                else:
                    unit = u.rad
                c = apycoords.SkyCoord(
                    ra * unit, dec * unit, equinox=epoch, frame="icrs"
                )
                c = c.transform_to(apycoords.Galactic)
                return c.l.to(unit).value, c.b.to(unit).value

            if self.radecgrid is True:
                for i in [-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75]:
                    ras = np.linspace(0, 360, 360)
                    des = np.linspace(i, i, 360)
                    l, b = radec_to_lb(ras, des, degree=True)
                    l = -(l + 180) % (2 * 180) - 180
                    if np.max(np.diff(l)) > 100.0:
                        idx = np.argmax(np.diff(l)) + 1
                        l = np.concatenate([l[idx:], l[:idx]])
                        b = np.concatenate([b[idx:], b[:idx]])
                    self.ax.plot(
                        self._fake_rad2deg(np.deg2rad(l)),
                        self._fake_rad2deg(np.deg2rad(b)),
                        c=self._opposite_color,
                        lw=grid_width,
                        ls=grid_style,
                        zorder=3,
                        alpha=grad_alpha,
                    )

                for i in [30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]:
                    ras = np.linspace(i, i, 360)
                    des = np.linspace(-75, 75, 360)
                    l, b = radec_to_lb(ras, des, degree=True)
                    l = -(l + 180) % (2 * 180) - 180
                    if np.max(np.diff(l)) > 0.25:
                        idx = np.argmax(np.diff(l))
                        idx = np.argmax(l) + 1
                        l = np.concatenate([l[idx:], l[:idx]])
                        b = np.concatenate([b[idx:], b[:idx]])
                        idx_break = np.argmax(np.diff(l))
                        self.ax.plot(
                            self._fake_rad2deg(np.deg2rad(l[:idx_break])),
                            self._fake_rad2deg(np.deg2rad(b[:idx_break])),
                            c=self._opposite_color,
                            lw=grid_width,
                            ls=grid_style,
                            zorder=3,
                            alpha=grad_alpha,
                        )
                        self.ax.plot(
                            self._fake_rad2deg(np.deg2rad(l[idx_break + 1 :])),
                            self._fake_rad2deg(np.deg2rad(b[idx_break + 1 :])),
                            c=self._opposite_color,
                            lw=grid_width,
                            ls=grid_style,
                            zorder=3,
                            alpha=grad_alpha,
                        )
                    else:
                        self.ax.plot(
                            self._fake_rad2deg(np.deg2rad(l)),
                            self._fake_rad2deg(np.deg2rad(b)),
                            c=self._opposite_color,
                            lw=grid_width,
                            ls=grid_style,
                            zorder=3,
                            alpha=grad_alpha,
                        )

            if self.eclgrid is True:

                def ecl_to_lb(elon, elat):
                    """
                    elon and elat in radian
                    """
                    e = 23.43928083333333 / 180 * np.pi
                    atan_top = np.sin(elon) * np.cos(e) - np.tan(elat) * np.sin(e)
                    atan_bottom = np.cos(elon)
                    ra = np.arctan(atan_top / atan_bottom)
                    dec = np.arcsin(
                        np.sin(elat) * np.cos(e)
                        + np.cos(elat) * np.sin(e) * np.sin(elon)
                    )
                    case_1_idx = (atan_top > 0) & (atan_bottom < 0)
                    case_2_idx = (atan_top < 0) & (atan_bottom > 0)
                    case_3_idx = (atan_top < 0) & (atan_bottom < 0)
                    ra[case_1_idx] += np.pi
                    ra[case_2_idx] += 2 * np.pi
                    ra[case_3_idx] += 3 * np.pi
                    l, b = radec_to_lb(ra, dec)
                    l = -(l + np.pi) % (2 * np.pi) - np.pi
                    return l, b

                for i in [-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75]:
                    elon = np.linspace(0, 360, 360)
                    elat = np.linspace(i, i, 360)
                    l, b = ecl_to_lb(np.deg2rad(elon), np.deg2rad(elat))
                    if np.max(np.diff(l)) > 2.0:
                        idx = np.argmax(np.diff(l)) + 1
                        l = np.concatenate([l[idx:], l[:idx]])
                        b = np.concatenate([b[idx:], b[:idx]])
                    self.ax.plot(
                        self._fake_rad2deg(l),
                        self._fake_rad2deg(b),
                        c=self._opposite_color,
                        lw=grid_width,
                        ls=grid_style,
                        zorder=3,
                        alpha=grad_alpha,
                    )

                for i in [30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]:
                    elon = np.linspace(i, i, 360)
                    elat = np.linspace(-75, 75, 360)
                    l, b = ecl_to_lb(np.deg2rad(elon), np.deg2rad(elat))
                    if np.max(np.diff(l)) > 0.004:
                        idx = np.argmax(np.diff(l))
                        idx = np.argmax(l) + 1
                        l = np.concatenate([l[idx:], l[:idx]])
                        b = np.concatenate([b[idx:], b[:idx]])
                        idx_break = np.argmax(np.diff(l))
                        self.ax.plot(
                            self._fake_rad2deg(l[:idx_break]),
                            self._fake_rad2deg(b[:idx_break]),
                            c=self._opposite_color,
                            lw=grid_width,
                            ls=grid_style,
                            zorder=3,
                            alpha=grad_alpha,
                        )
                        self.ax.plot(
                            self._fake_rad2deg(l[idx_break + 1 :]),
                            self._fake_rad2deg(b[idx_break + 1 :]),
                            c=self._opposite_color,
                            lw=grid_width,
                            ls=grid_style,
                            zorder=3,
                            alpha=grad_alpha,
                        )
                    else:
                        self.ax.plot(
                            self._fake_rad2deg(l),
                            self._fake_rad2deg(b),
                            c=self._opposite_color,
                            lw=grid_width,
                            ls=grid_style,
                            zorder=3,
                            alpha=grad_alpha,
                        )

    def show(self, *args, **kwargs):
        if self.fig is None:
            raise AttributeError("Nothing to show, please plot some data first")
        else:
            if self.tight_layout is True:
                self.fig.tight_layout()
            self.fig.show(*args, **kwargs)

    def savefig(self, file="MWSkyMap.png", **kwargs):
        if self.tight_layout is True:
            self.fig.tight_layout()
        # this is a pylab method
        self.fig.savefig(file, **kwargs)

    def mw_scatter(self, ra, dec, c="r", **kwargs):
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
            if hasattr(c[0], "__len__"):
                color = c[0]
                cbar_label = c[1]
                self.cbar_flag = True
                if isinstance(color, u.quantity.Quantity):
                    color = color.to(self.unit).value
            else:
                color = c
        else:
            color = c

        mappable = self.ax.scatter(
            ra,
            dec,
            zorder=3,
            s=self.s,
            c=color,
            cmap=plt.get_cmap(self.cmap) if self.cbar_flag else None,
            rasterized=True,
            **kwargs,
        )
        if self.projection == "equirectangular":
            self.ax.imshow(
                self.bg_img, zorder=0, extent=self._ext, alpha=0.0, rasterized=True
            )
            self.ax.set_aspect("equal")
        else:
            self.ax.imshow(
                self.bg_img,
                zorder=0,
                extent=self._ext,
                alpha=self.imalpha,
                rasterized=True,
                aspect=self.ax.get_aspect(),
                transform=self.ax.transAxes,
            )
        if self.cbar_flag is True:
            if self.projection == "equirectangular":
                divider = make_axes_locatable(self.ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = self.fig.colorbar(mappable, cax=cax)
            else:
                cbar = self.fig.colorbar(mappable, ax=self.ax)
            if self.clim is not None:
                cbar.set_clim(self.clim)

    def scatter(self, ra, dec, c="r", *args, **kwargs):
        ra, dec = self.radec_unit_check(ra, dec)
        self.initialize_mwplot()
        if kwargs.get("s") is None:
            kwargs["s"] = self.s
        self.ax.scatter(ra, dec, c=c, zorder=3, rasterized=True, *args, **kwargs)
        # just want to set the location right, we dont need image again
        if self.projection == "equirectang ular":
            self.ax.imshow(
                self.bg_img, zorder=0, extent=self._ext, alpha=0.0, rasterized=True
            )
        else:
            self.ax.imshow(
                self.bg_img,
                zorder=0,
                extent=self._ext,
                alpha=self.imalpha,
                rasterized=True,
                aspect=self.ax.get_aspect(),
                transform=self.ax.transAxes,
            )
        if kwargs.get("label") is not None:
            self.ax.legend(loc="best", markerscale=kwargs["s"])

    def scatter_annotate(
        self,
        text,
        position,
        arrowprops=dict(facecolor="black", width=1.0, headwidth=6.0, headlength=6.0),
        fontsize=15,
        bbox=dict(pad=2),
        **kwargs,
    ):
        """
        Plot annotation with scatter

        :History: 2022-Jan-02 - Written - Henry Leung (University of Toronto)
        """
        if isinstance(position, apycoords.SkyCoord):
            position = self.skycoord_radec(position)
        position_wo_unit = self.xy_unit_check(position[0], position[1])
        position_text = np.add(position_wo_unit, 10)
        if isinstance(text, list):
            for t, p, pou, pt in zip(text, position, position_wo_unit, position_text):
                self.scatter(p[0], p[1])
                self.ax.annotate(
                    t,
                    xy=pou,
                    xytext=pt,
                    arrowprops=arrowprops,
                    fontsize=fontsize,
                    bbox=bbox,
                    **kwargs,
                )
        else:
            self.scatter(position[0], position[1])
            self.ax.annotate(
                text,
                xy=position_wo_unit,
                xytext=position_text,
                arrowprops=arrowprops,
                fontsize=fontsize,
                bbox=bbox,
                **kwargs,
            )

    def annotate(self, *args, **kwargs):
        """
        Plot annotation

        :History: 2022-Jan-02 - Written - Henry Leung (University of Toronto)
        """
        return self.ax.annotate(*args, **kwargs)
