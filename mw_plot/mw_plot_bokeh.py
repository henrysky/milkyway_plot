try:
    import bokeh
except ImportError:
    raise ImportError(
        "Bokeh is not installed. Please install Bokeh to use this feature"
    )

import requests
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
from mw_plot.mw_plot_masters import MWPlotMaster, MWSkyMapMaster

from bokeh.plotting import figure, show
from bokeh.io import output_file, save, output_notebook
from bokeh.models import Range1d


def to_bokeh_img(imgarray):
    M, N, _ = imgarray.shape
    img = np.empty((M, N), dtype=np.uint32)
    view = img.view(dtype=np.uint8).reshape((M, N, 4))
    view[:, :, 0] = imgarray[:, :, 0]  # copy red channel
    view[:, :, 1] = imgarray[:, :, 1]  # copy blue channel
    view[:, :, 2] = imgarray[:, :, 2]  # copy green channel
    view[:, :, 3] = 255

    img = img[::-1]  # flip for Bokeh
    return img


class MWPlotBokeh(MWPlotMaster):
    """
    MWPlot Brokeh class plotting with Bokeh

    :param mode: whether plot edge-on or face-on milkyway
    :type mode: string, either 'face-on' or 'edge-on'
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
    """

    def __init__(
        self,
        mode="face-on",
        center=(0, 0) * u.kpc,
        radius=90750 * u.lyr,
        unit=u.kpc,
        coord="galactic",
        annotation=True,
        rot90=0,
        grayscale=False,
        r0=8.125,
    ):
        super().__init__(
            grayscale=grayscale,
            annotation=annotation,
            rot90=rot90,
            coord=coord,
            mode=mode,
            r0=r0,
            center=center,
            radius=radius,
            unit=unit,
            figsize=None,
            dpi=None,
        )

        # prepossessing procedure
        self._unit_english = self._unit.short_names[0]
        if self._center.unit is not None and self._radius.unit is not None:
            self._center = self._center.to(self._unit)
            self._radius = self._radius.to(self._unit)

        self.images_read()
        self.s = 1.0

        TOOLS = "pan, wheel_zoom, box_zoom, reset, save, box_select"

        self.bokeh_fig = figure(
            title="",
            tools=TOOLS,
            x_range=Range1d(
                self._ext[0],
                self._ext[1],
                bounds=[
                    min(self._ext[0], self._ext[1]),
                    max(self._ext[0], self._ext[1]),
                ],
            ),
            y_range=Range1d(
                self._ext[2],
                self._ext[3],
                bounds=[
                    min(self._ext[2], self._ext[3]),
                    max(self._ext[2], self._ext[3]),
                ],
            ),
            width=1000,
            height=1000,
        )
        if (
            requests.head(self._gh_img_url, allow_redirects=True).status_code == -9999
        ):  # connection successful
            # disabled currently because rotation and grayscale wont work
            self.bokeh_fig.image_url(
                url=[self._gh_img_url],
                x=self._ext[0],
                y=self._ext[2],
                w=abs(self._ext[1] - self._ext[0]),
                h=abs(self._ext[3] - self._ext[2]),
                anchor="bottom_left",
            )
        else:
            self._img = to_bokeh_img(self._img)
            self.bokeh_fig.image_rgba(
                image=[self._img],
                x=self._ext[0],
                y=self._ext[2],
                dw=abs(self._ext[1] - self._ext[0]),
                dh=abs(self._ext[3] - self._ext[2]),
            )

        self.bokeh_fig.xaxis.axis_label = (
            f"{self._coord_english} ({self._unit_english})"
        )
        self.bokeh_fig.yaxis.axis_label = (
            f"{self._coord_english} ({self._unit_english})"
        )

    def scatter(self, x, y, *args, **kwargs):
        x, y = self.xy_unit_check(x, y)
        if kwargs.get("s") is None:
            kwargs["s"] = self.s
        self.bokeh_fig.circle(x, y, size=self.s)

    def show(self, notebook=True):
        if self._in_jupyter and notebook:
            output_notebook()
        else:
            pass
        show(self.bokeh_fig)

    def savefig(self, file="MWPlot.html"):
        output_file(file)
        save(self.bokeh_fig)


class MWSkyMapBokeh(MWSkyMapMaster):
    """
    MWSkyMapBokeh class plotting with Bokeh

    :param center: Coordinates of the center of the plot with astropy degree/radian units
    :type center: astropy.Quantity
    :param radius: Radius of the plot with astropy degree/radian units
    :type radius: astropy.Quantity
    :param grayscale: whether to use grayscale background
    :type grayscale: bool
    """

    def __init__(
        self, center=(0, 0) * u.deg, radius=(180, 90) * u.deg, grayscale=False
    ):
        super().__init__(
            grayscale=grayscale,
            projection="equirectangular",
            center=center,
            radius=radius,
            figsize=None,
            dpi=None,
        )
        self._unit = u.degree
        self.s = 1.0

        # preprocessing
        if self._center.unit is not None and self._radius.unit is not None:
            self._center = self._center.to(self._unit)
            self._radius = self._radius.to(self._unit)

        if (self._center[0] + self._radius[0]).value > 180 or (
            self._center[0] - self._radius[0]
        ).value < -180:
            raise ValueError(
                "The border of the width will be outside the range of -180 to 180 which is not allowed\n"
            )
        if (self._center[1] + self._radius[1]).value > 90 or (
            self._center[1] - self._radius[1]
        ).value < -90:
            raise ValueError(
                "The border of the height will be outside the range of -90 to 90 which is not allowed"
            )
        if self._radius[0] <= 0 or self._radius[0] <= 0:
            raise ValueError("Radius cannot be negative or 0")

        self.images_read()
        self.s = 1.0

        TOOLS = "pan, wheel_zoom, box_zoom, reset, save, box_select"

        self.bokeh_fig = figure(
            title="",
            tools=TOOLS,
            x_range=Range1d(
                self._ext[0],
                self._ext[1],
                bounds=[
                    min(self._ext[0], self._ext[1]),
                    max(self._ext[0], self._ext[1]),
                ],
            ),
            y_range=Range1d(
                self._ext[2],
                self._ext[3],
                bounds=[
                    min(self._ext[2], self._ext[3]),
                    max(self._ext[2], self._ext[3]),
                ],
            ),
            width=1000,
            height=500,
        )

        if (
            requests.head(self._gh_img_url, allow_redirects=True).status_code == -9999
        ):  # connection successful
            # disabled currently because rotation and grayscale wont work
            self.bokeh_fig.image_url(
                url=[self._gh_img_url],
                x=self._ext[0],
                y=self._ext[2],
                w=abs(self._ext[1] - self._ext[0]),
                h=abs(self._ext[3] - self._ext[2]),
                anchor="bottom_left",
            )
        else:
            self._img = to_bokeh_img(self._img)
            self.bokeh_fig.image_rgba(
                image=[self._img],
                x=self._ext[0],
                y=self._ext[2],
                dw=abs(self._ext[1] - self._ext[0]),
                dh=abs(self._ext[3] - self._ext[2]),
            )

        self.bokeh_fig.xaxis.axis_label = "Galactic Longitude (Degree)"
        self.bokeh_fig.yaxis.axis_label = "Galactic Latitude (Degree)"

    def scatter(self, ra, dec, *args, **kwargs):
        ra, dec = self.radec_unit_check(ra, dec)
        if kwargs.get("s") is None:
            kwargs["s"] = self.s
        self.bokeh_fig.circle(ra, dec, size=self.s)

    def show(self, notebook=True):
        if self._in_jupyter and notebook:
            output_notebook()
        else:
            pass
        show(self.bokeh_fig)

    def savefig(self, file="MWSkyMap.html"):
        output_file(file)
        save(self.bokeh_fig)
