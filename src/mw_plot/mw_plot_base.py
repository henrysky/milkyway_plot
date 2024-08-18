import importlib.util
import pathlib
import warnings
from abc import ABC, abstractmethod
from dataclasses import dataclass

import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from mw_plot.utils import rgb2gray


@dataclass
class MWImage:
    """
    Dataclass to keep track of background images
    """

    filename: str
    citation: str

    @property
    def img_path(self) -> pathlib.Path:
        # get path to the mw_plot directory
        return pathlib.Path(importlib.util.find_spec("mw_plot").origin).parent.joinpath(
            self.filename
        )


class MWPlotCommon(ABC):
    """
    Common class for MWPlotBase and MWSkyMapBase
    """

    def __init__(self):
        self._gh_imgbase_url = (
            "https://github.com/henrysky/milkyway_plot/raw/master/mw_plot/"
        )

        # check if running in browser-based ipython (aka jupyter)
        self._in_jupyter = False
        try:
            from IPython import get_ipython
        except ImportError:
            # the case where ipython is not installed
            pass
        else:
            if (ip := get_ipython()) is None:
                # the case where ipython is installed but not in ipython runtime
                pass
            else:
                # only in jupyter notebook/lab has the kernel trait
                self._in_jupyter = ip.has_trait("kernel")

        # all the background images
        self._MW_IMAGES = {
            "MW_bg_unannotate": MWImage(
                filename="MW_bg_unannotate.jpg",
                citation="NASA/JPL-Caltech/R. Hurt (SSC/Caltech)",
            ),
            "MW_bg_annotate": MWImage(
                filename="MW_bg_annotate.jpg",
                citation="NASA/JPL-Caltech/R. Hurt (SSC/Caltech)",
            ),
            "MW_edgeon_edr3_unannotate": MWImage(
                filename="MW_edgeon_edr3_unannotate.jpg", citation="ESA/Gaia/DPAC"
            ),
            "MW_fermi_gamma": MWImage(
                filename="MW_fermi_gamma.jpg",
                citation="NASA/DOE/Fermi LAT Collaboration",
            ),
            "MW_2mass": MWImage(
                filename="MW_2mass.jpg",
                citation="2MASS/IPAC/Caltech/University of Massachusetts",
            ),
            "MW_farinfrared": MWImage(
                filename="MW_farinfrared.jpg",
                citation="WISE/NASA/JPL-Caltech/UCLA & IRAS/NASA",
            ),
        }

    @abstractmethod
    def read_bg_img(self):
        # class to read images
        pass

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

    def xy_unit_check(self, x, y, checkrot=True):
        if not isinstance(x, u.quantity.Quantity) or not isinstance(
            y, u.quantity.Quantity
        ):
            raise TypeError("All numbers must be astropy Quantity")
        if x.unit is None and y.unit is None:
            raise TypeError("All numbers must carry astropy unit")
        else:
            x = x.to(self._unit).value
            y = y.to(self._unit).value

        # check if rotation is 90deg or 270deg
        if checkrot:  # nested, do not unnest
            if self._rot90 % 2 == 1:
                x, y = y, x
        return x, y

    @property
    def reference(self):
        return self.reference_str


class MWPlotBase(MWPlotCommon):
    """
    MWPlot base class to plot the face-on Milky Way

    Parameters
    ----------
    grayscale : bool
        Whether to use grayscale background
    annotation : bool
        Whether use a milkyway background with annotation
    rot90 : int
        Number of 90 degree rotation
    coord : str
        'galactocentric' or 'galactic'
    r0 : float
        Distance to galactic center in kpc
    center : tuple
        Coordinates of the center of the plot with astropy units
    radius : float
        Radius of the plot with astropy units
    unit : astropy units
        Astropy units
    figsize : tuple
        Figure size
    dpi : int
        Dots per inch
    """

    def __init__(
        self,
        grayscale,
        annotation,
        rot90,
        coord,
        r0,
        center,
        radius,
        unit,
        figsize,
        dpi,
    ):
        super().__init__()
        self._coord = coord
        self._annotation = annotation
        self._rot90 = rot90
        self._grayscale = grayscale
        self.r0 = r0 * u.kpc
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
        self.reference_str = None

        # properties of the images
        self._pixels = 5600  # number of pixels in x and y axis
        self._ly_per_pixel = (self.r0 / 1078).to(u.lyr)  # light years per pixel resolution

    def lrbt_rot(self):
        """This function rotate matplolti's extent ordered LRBT"""
        l, r, b, t = self._ext[0], self._ext[1], self._ext[2], self._ext[3]
        if self._rot90 % 4 == 1:  # -90deg
            self._ext = [b, t, l, r]
        elif self._rot90 % 4 == 2:  # -180deg
            self._ext = [r, l, t, b]
        elif self._rot90 % 4 == 3:  # -270deg
            self._ext = [t, b, r, l]

    def read_bg_img(self):
        if self._annotation:
            img_obj = self._MW_IMAGES["MW_bg_annotate"]
        else:
            img_obj = self._MW_IMAGES["MW_bg_unannotate"]
        self._img = plt.imread(img_obj.img_path)
        self.reference_str = img_obj.citation
        self._gh_img_url = self._gh_imgbase_url + img_obj.filename

        if self._grayscale:
            self._img = rgb2gray(self._img)

        if self._coord.lower() == "galactic":
            # shift the coord by r0 to the new coord system
            x_shift = self.r0
            self._center[0] += x_shift
            self._coord_english = "Galactic Coordinates"
        elif self._coord.lower() == "galactocentric":
            x_shift = 0.0 * u.kpc
            self._coord_english = "Galactocentric Coordinates"
        else:
            raise ValueError(
                "Unknown coordinates, can only be 'galactic' or 'galactocentric'"
            )

        if not isinstance(self._center, u.quantity.Quantity) and not isinstance(
            self._radius, u.quantity.Quantity
        ):
            warnings.warn(
                f"You did not specify units for center and radius, assuming the unit is {self._unit.long_names[0]}"
            )
            if not isinstance(self._center, u.quantity.Quantity):
                self._center = self._center * self._unit
            if not isinstance(self._radius, u.quantity.Quantity):
                self._radius = self._radius * self._unit

        self._resolution = self._ly_per_pixel.to(self._unit)
        self._radius = self._radius.to(self._unit)

        # convert physical unit to pixel unit
        pixel_radius = int((self._radius / self._resolution).value)
        pixel_center = [
            int((self._pixels / 2 + self._center[0] / self._resolution).value),
            int((self._pixels / 2 - self._center[1] / self._resolution).value),
        ]

        # get the pixel coordinates
        x_left_px = pixel_center[0] - pixel_radius
        x_right_px = pixel_center[0] + pixel_radius
        y_top_px = self._pixels - pixel_center[1] - pixel_radius
        y_bottom_px = self._pixels - pixel_center[1] + pixel_radius

        # decide whether it needs to fill black pixels because the range outside the pre-compiled images
        if np.all(
            np.array(
                [
                    x_left_px,
                    self._pixels - x_right_px,
                    y_top_px,
                    self._pixels - y_bottom_px,
                ]
            )
            >= 0
        ):
            self._img = self._img[y_top_px:y_bottom_px, x_left_px:x_right_px]
        else:
            # create a black/white image first with 3 channel with the same data type
            if self._grayscale:
                black_img = (
                    np.ones(
                        (pixel_radius * 2, pixel_radius * 2, 3), dtype=self._img.dtype
                    )
                    * 255
                )
            else:
                black_img = np.zeros(
                    (pixel_radius * 2, pixel_radius * 2, 3), dtype=self._img.dtype
                )

        self._img = np.rot90(self._img, self._rot90)
        self._ext = [
            (self._center[0] - self._radius - x_shift).value,
            (self._center[0] + self._radius - x_shift).value,
            (self._center[1] - self._radius).value,
            (self._center[1] + self._radius).value,
        ]

        self._img_ext = [
            np.max([self._ext[0], (self._center[0] - ((5600 / 2) * self._ly_per_pixel).to(self._unit) - x_shift).value]),
            np.min([self._ext[1], (self._center[0] + ((5600 / 2) * self._ly_per_pixel).to(self._unit)  - x_shift).value]),
            np.max([self._ext[2], (self._center[1] - ((5600 / 2) * self._ly_per_pixel).to(self._unit)).value]),
            np.min([self._ext[3], (self._center[1] + ((5600 / 2) * self._ly_per_pixel).to(self._unit)).value]),
        ]

        self._img = self._img
        self._aspect = (
            self._img.shape[0]
            / float(self._img.shape[1])
            * ((self._ext[1] - self._ext[0]) / (self._ext[3] - self._ext[2]))
        )
        self._aspect = np.abs(self._aspect)

        self.lrbt_rot()

        return None


class MWSkyMapBase(MWPlotCommon):
    """
    MWSkyMap base class to plot the sky map
    """

    def __init__(
        self,
        grayscale,
        projection,
        wavelength,
        center,
        radius,
        figsize,
        dpi,
    ):
        super().__init__()
        self._projection = projection

        if self._projection in (
            allowed_proj := [
                "equirectangular",
                "aitoff",
                "hammer",
                "mollweide",
            ]
        ):
            pass
        # projection that is not suitable for sky map but available from matplotlib
        elif self._projection in ["lambert", "polar", "rectilinear"]:
            raise NotImplementedError(
                f"`{self._projection}` projection is not implemented for sky map"
            )
        else:
            raise ValueError(
                f"Unknown projection `{self._projection}`, allowed values are: {allowed_proj}"
            )

        if wavelength in (
            allowed_wavelength := ["gamma", "optical", "infrared", "far-infrared"]
        ):
            self.wavlength = wavelength
        else:
            raise ValueError(
                f"Unknown wavelength, allowed values are: {allowed_wavelength}"
            )

        self._center = center
        self._radius = radius
        self._grayscale = grayscale
        self.figsize = figsize
        self.dpi = dpi
        self._initialized = False

        self._opposite_color = "white"
        if self._grayscale:
            self._opposite_color = "black"

    def read_bg_img(self):
        if self.wavlength == "optical":
            img_key = "MW_edgeon_edr3_unannotate"
        elif self.wavlength == "gamma":
            img_key = "MW_fermi_gamma"
        elif self.wavlength == "infrared":
            img_key = "MW_2mass"
        elif self.wavlength == "far-infrared":
            img_key = "MW_farinfrared"
        else:
            raise ValueError("Unknown wavelength")
        img_obj = self._MW_IMAGES[img_key]
        self._img = plt.imread(img_obj.img_path)
        self.reference_str = img_obj.citation

        # find center pixel and radius pixel
        y_img_center = self._img.shape[0] // 2 - int(
            (self._img.shape[0] / 180) * self._center[1].value
        )
        y_radious_px = int((self._img.shape[0] / 180) * self._radius[1].value)
        x_img_center = (
            int((self._img.shape[1] / 360) * self._center[0].value) + self._img.shape[0]
        )
        x_radious_px = int((self._img.shape[1] / 360) * self._radius[0].value)

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

        self._gh_img_url = self._gh_imgbase_url + img_obj.filename

        return None

    def radec_unit_check(self, ra, dec):
        if not isinstance(ra, u.quantity.Quantity) or not isinstance(
            dec, u.quantity.Quantity
        ):
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
