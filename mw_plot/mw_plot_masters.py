import os
import numpy as np
import pylab as plt
from abc import ABC, abstractmethod
import astropy.units as u
import astropy.coordinates as coord
import warnings


def rgb2gray(rgb):
    """
    Change RGB color image into grayscale in RGB representation using colorimetric (perceptual luminance-preserving) conversion

    :param rgb: NumPy array of the RGB image
    :return: NumPy array of grayscale image, same shape as input
    """
    r, g, b = rgb[:, :, 0] / 255, rgb[:, :, 1] / 255, rgb[:, :, 2] / 255
    r = np.where(r <= 0.04045, r / 12.92, ((r + 0.055) / 1.055) ** 2.4)
    g = np.where(g <= 0.04045, g / 12.92, ((g + 0.055) / 1.055) ** 2.4)
    b = np.where(b <= 0.04045, b / 12.92, ((b + 0.055) / 1.055) ** 2.4)
    gray = np.atleast_3d(1 - (0.2126 * r + 0.7152 * g + 0.0722 * b))
    return np.repeat(255 * gray, 3, axis=2).astype(int)


class MWPlotMaster(ABC):
    """
    MWPlot master class
    """

    def __init__(
        self,
        grayscale,
        annotation,
        rot90,
        coord,
        mode,
        r0,
        center,
        radius,
        unit,
        figsize,
        dpi,
    ):
        self.__coord = coord
        self.__annotation = annotation
        self.__rot90 = rot90
        self._grayscale = grayscale
        self.r0 = r0 * u.kpc
        self.mode = mode
        self._initialized = False
        self.figsize = figsize
        self.dpi = dpi

        # user should not change these values anyway
        self._ext = None
        self._img = None
        self._img_fname = None
        self._gh_img_url = None
        self._center = center
        self._radius = radius
        self._unit = unit
        self._gh_imgbase_url = (
            "https://github.com/henrysky/milkyway_plot/raw/master/mw_plot/"
        )

        # Fixed value
        if self.mode == "face-on":
            self.__pixels = 5600
            self.__resolution = (self.r0 / 1078).to(u.lyr)
        elif self.mode == "edge-on":
            self.__pixels = 6500
            self.__resolution = 15.384615846 * u.lyr
        else:
            raise LookupError(
                f"Unknown mode '{self.mode}', can only be 'edge-on' or 'face-on'"
            )

        # check if running in browser-based ipython (aka jupyter)
        self._in_jupyter = False
        try:
            from IPython import get_ipython

            ip = get_ipython()
        except ImportError:
            pass
        else:
            if hasattr(ip, "has_trait"):
                if ip.has_trait("kernel"):
                    self._in_jupyter = True

    def xy_unit_check(self, x, y, checkrot=True):
        if not type(x) == u.quantity.Quantity or not type(y) == u.quantity.Quantity:
            raise TypeError("All numbers must be astropy Quantity")
        if x.unit is None and y.unit is None:
            raise TypeError("All numbers must carry astropy unit")
        else:
            x = x.to(self._unit).value
            y = y.to(self._unit).value

        # check if rotation is 90deg or 270deg
        if checkrot:  # nested, do not unnest
            if self.__rot90 % 2 == 1:
                x, y = y, x
        return x, y

    def lrbt_rot(self):
        """This function rotate matplolti's extent ordered LRBT"""
        l, r, b, t = self._ext[0], self._ext[1], self._ext[2], self._ext[3]
        if self.__rot90 % 4 == 1:  # -90deg
            self._ext = [b, t, l, r]
        elif self.__rot90 % 4 == 2:  # -180deg
            self._ext = [r, l, t, b]
        elif self.__rot90 % 4 == 3:  # -270deg
            self._ext = [t, b, r, l]

    def images_read(self):
        if self.mode == "edge-on":
            image_filename = "MW_edgeon_edr3_unannotate.jpg"
            path = os.path.join(os.path.dirname(__file__), image_filename)
            img = np.zeros((6500, 6500, 3), dtype=int)
            img[1625:4875, :, :] = plt.imread(path)
        elif self.__annotation is False:
            image_filename = "MW_bg_unannotate.jpg"
            path = os.path.join(os.path.dirname(__file__), image_filename)
            img = plt.imread(path)
        else:
            image_filename = "MW_bg_annotate.jpg"
            path = os.path.join(os.path.dirname(__file__), image_filename)
            img = plt.imread(path)
        self._img_fname = image_filename
        self._gh_img_url = self._gh_imgbase_url + self._img_fname

        if self._grayscale:
            img = rgb2gray(img)

        if self.__coord.lower() == "galactic":
            # shift the coord by r0 to the new coord system
            x_shift = self.r0
            self._center[0] += x_shift
            self._coord_english = "Galactic Coordinates"
        elif self.__coord.lower() == "galactocentric":
            x_shift = 0.0 * u.kpc
            self._coord_english = "Galactocentric Coordinates"
        else:
            raise ValueError(
                "Unknown coordinates, can only be 'galactic' or 'galactocentric'"
            )

        if (
            not type(self._center) == u.quantity.Quantity
            and not type(self._radius) == u.quantity.Quantity
        ):
            warnings.warn(
                f"You did not specify units for center and radius, assuming the unit is {self._unit.long_names[0]}"
            )
            if not type(self._center) == u.quantity.Quantity:
                self._center = self._center * self._unit
            if not type(self._radius) == u.quantity.Quantity:
                self._radius = self._radius * self._unit

        self.__resolution = self.__resolution.to(self._unit)
        self._radius = self._radius.to(self._unit)

        # convert physical unit to pixel unit
        pixel_radius = int((self._radius / self.__resolution).value)
        pixel_center = [
            int((self.__pixels / 2 + self._center[0] / self.__resolution).value),
            int((self.__pixels / 2 - self._center[1] / self.__resolution).value),
        ]

        # get the pixel coordinates
        x_left_px = pixel_center[0] - pixel_radius
        x_right_px = pixel_center[0] + pixel_radius
        y_top_px = self.__pixels - pixel_center[1] - pixel_radius
        y_bottom_px = self.__pixels - pixel_center[1] + pixel_radius

        # decide whether it needs to fill black pixels because the range outside the pre-compiled images
        if np.all(
            np.array(
                [
                    x_left_px,
                    self.__pixels - x_right_px,
                    y_top_px,
                    self.__pixels - y_bottom_px,
                ]
            )
            >= 0
        ):
            img = img[y_top_px:y_bottom_px, x_left_px:x_right_px]
        else:
            # create a black/white image first with 3 channel with the same data type
            if self._grayscale:
                black_img = (
                    np.ones((pixel_radius * 2, pixel_radius * 2, 3), dtype=img.dtype)
                    * 255
                )
            else:
                black_img = np.zeros(
                    (pixel_radius * 2, pixel_radius * 2, 3), dtype=img.dtype
                )

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

            # fill the black/white image with the background image
            black_img[
                top_exceed_px : top_exceed_px + img.shape[0],
                left_exceed_px : left_exceed_px + img.shape[1],
                :,
            ] = img

            # Set the images as the filled black-background image
            img = np.array(black_img)

        img = np.rot90(img, self.__rot90)
        self._ext = [
            (self._center[0] - self._radius - x_shift).value,
            (self._center[0] + self._radius - x_shift).value,
            (self._center[1] - self._radius).value,
            (self._center[1] + self._radius).value,
        ]

        if self.mode == "edge-on":
            self._ext[2] *= -1
            self._ext[3] *= -1

        self._img = img
        self._aspect = (
            img.shape[0]
            / float(img.shape[1])
            * ((self._ext[1] - self._ext[0]) / (self._ext[3] - self._ext[2]))
        )
        self._aspect = np.abs(self._aspect)

        self.lrbt_rot()

        return None

    @staticmethod
    def skycoord_xy(skycoord):
        # convert astropy SkyCoord to cartesian x, y
        return [skycoord.cartesian.x, skycoord.cartesian.y]

    @staticmethod
    def skycoord_radec(skycoord):
        # convert astropy SkyCoord to list
        if not hasattr(skycoord, "ra"):
            skycoord = skycoord.icrs
        return [skycoord.ra.deg * u.deg, skycoord.dec.deg * u.deg]


class MWSkyMapMaster(MWPlotMaster):
    """
    MWSkyMap master class
    """

    def __init__(self, grayscale, projection, center, radius, figsize, dpi, grid=False):
        if projection in [
            "equirectangular",
            "aitoff",
            "hammer",
            "lambert",
            "mollweide",
        ]:
            self._projection = projection
        else:
            raise ValueError("Unknown projection")

        self._center = center
        self._radius = radius
        self._grayscale = grayscale
        self.figsize = figsize
        self.dpi = dpi
        self.grid = grid
        self._gh_imgbase_url = (
            "https://github.com/henrysky/milkyway_plot/raw/master/mw_plot/"
        )
        self._initialized = False

        self._opposite_color = "white"
        if self._grayscale:
            self._opposite_color = "black"

        # check if running in browser-based ipython (aka jupyter)
        self._in_jupyter = False
        try:
            from IPython import get_ipython

            ip = get_ipython()
        except ImportError:
            pass
        else:
            if hasattr(ip, "has_trait"):
                if ip.has_trait("kernel"):
                    self._in_jupyter = True

    def images_read(self):
        image_filename = "MW_edgeon_edr3_unannotate.jpg"
        path = os.path.join(os.path.dirname(__file__), image_filename)
        img = plt.imread(path)
        self._img = img

        # find center pixel and radius pixel
        y_img_center = 1625 - int((3250 / 180) * self._center[1].value)
        y_radious_px = int((3250 / 180) * self._radius[1].value)
        x_img_center = int((6500 / 360) * self._center[0].value) + 3250
        x_radious_px = int((6500 / 360) * self._radius[0].value)

        self._ext = [
            (self._center[0] - self._radius[0]).value,
            (self._center[0] + self._radius[0]).value,
            (self._center[1] - self._radius[1]).value,
            (self._center[1] + self._radius[1]).value,
        ]

        self._img = self._img[
            (y_img_center - y_radious_px) : (y_img_center + y_radious_px),
            (x_img_center - x_radious_px) : (x_img_center + x_radious_px),
            :,
        ]

        if self._grayscale:
            self._img = rgb2gray(self._img)

        self._img_fname = image_filename
        self._gh_img_url = self._gh_imgbase_url + self._img_fname

        return None

    def radec_unit_check(self, ra, dec):
        if not type(ra) == u.quantity.Quantity or not type(dec) == u.quantity.Quantity:
            raise TypeError("Both RA and DEC must carry astropy's unit")
        else:
            if ra.unit is not None and dec.unit is not None:
                ra = ra.to(self._unit)
                dec = dec.to(self._unit)
                c_icrs = coord.SkyCoord(ra=ra, dec=dec, frame="icrs")
                if self._projection == "equirectangular":
                    ra = coord.Angle(-c_icrs.galactic.l).wrap_at(180 * u.degree).value
                    dec = coord.Angle(c_icrs.galactic.b).value
                else:  # projection requires radian instead of degree
                    ra = (
                        coord.Angle(-c_icrs.galactic.l)
                        .wrap_at(180 * u.degree)
                        .to(u.radian)
                        .value
                    )
                    dec = coord.Angle(c_icrs.galactic.b).to(u.radian).value
            else:
                raise TypeError(
                    "Both x, y, center and radius must carry astropy's unit"
                )

        return ra, dec
