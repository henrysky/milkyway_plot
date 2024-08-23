import importlib.util
import pathlib
from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import reduce
from typing import Optional

import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import PIL
import requests
from astropy import wcs as astropy_wcs
from astroquery.hips2fits import hips2fits
from astroquery.simbad import Simbad

from mw_plot.utils import rgb2gray

# global variable to store the HiPS metadata response
_HiPS_metadata_response = None

# increase respond timeout to 120s
hips2fits.timeout = 120


@dataclass
class MWImage:
    """
    Dataclass to keep track of pre-packaged background images
    """

    filename: str
    citation: str

    @property
    def img_path(self) -> pathlib.Path:
        # get path to the mw_plot directory
        return pathlib.Path(importlib.util.find_spec("mw_plot").origin).parent.joinpath(
            self.filename
        )

    @property
    def pillow_img(self):
        # get the image using PIL
        return PIL.Image.open(self.img_path)


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
        else:  # pragma: no cover
            if (ip := get_ipython()) is None:
                # the case where ipython is installed but not in ipython runtime
                pass
            else:
                # only in jupyter notebook/lab has the kernel trait
                self._in_jupyter = ip.has_trait("kernel")

        # all the packaged background images
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
    def read_bg_img(self):  # pragma: no cover
        # class to read images and set appropriate attributes like extent, aspect ratio etc
        pass

    @staticmethod
    def skycoord_xy(skycoord):
        # convert astropy SkyCoord to cartesian x, y
        return [skycoord.cartesian.x, skycoord.cartesian.y]

    @staticmethod
    def objname_to_coord(objname: str) -> coord.SkyCoord:
        """
        Convert object name to astropy coordinates

        Parameters
        ----------
        objname : str
            Name of an object

        Returns
        -------
        coord.SkyCoord
            Astropy SkyCoord object
        """
        simbad = Simbad()
        simbad.add_votable_fields("ra(d)", "dec(d)", "plx", "distance")
        result = simbad.query_object(objname)

        if result is None:
            raise ValueError(f"Object `{objname}` not found in Simbad")
        else:
            result = result[
            "RA_d", "DEC_d", "PLX_VALUE", "Distance_distance", "Distance_unit"
        ].filled(np.nan)  # result is a single object
            
        if len(result) != 1:
            raise ValueError(
                f"Multiple objects found for `{objname}` in Simbad but expected only one"
            )
        else:
            result = result[0]

        if np.isnan(result["PLX_VALUE"]) and np.isnan(result["Distance_distance"]):
            # in case the distance is not available
            distance = None
        elif np.isnan(result["PLX_VALUE"]):
            # only use distance if parallax is not available
            distance = (
                result["Distance_distance"] * u.Unit(result["Distance_unit"])
            ).to(u.kpc)
        else:  # use parallax to calculate distance
            distance = (1 / result["PLX_VALUE"]) * u.kpc

        return coord.SkyCoord(
            result["RA_d"] * u.deg,
            result["DEC_d"] * u.deg,
            distance=distance,
        )

    @staticmethod
    def skycoord_radec(skycoord):
        # convert astropy SkyCoord to list
        if not hasattr(skycoord, "ra"):
            skycoord = skycoord.icrs
        return [skycoord.ra.deg * u.deg, skycoord.dec.deg * u.deg]

    @staticmethod
    def parse_hips_background():
        global _HiPS_metadata_response
        # cache the response if it is not already cached
        if _HiPS_metadata_response is None:
            response = requests.get(
                "https://alasky.u-strasbg.fr/MocServer/query?*/P/*&get=record"
            )
            # only store the response if it is successful
            if response.status_code != 200:
                raise ConnectionError(
                    f"Failed to retrieve data. Status code: {response.status_code}"
                )
            # Split the content into chunks based on empty lines
            chunks = response.text.split("\n\n")

            # Store the chunks in a list
            chunks_list = [chunk.strip() for chunk in chunks if chunk.strip()]

            dist_ls = []
            for chunk in chunks_list:
                # parse the chunk into a dictionary
                data = {}
                for line in chunk.splitlines():
                    key, value = line.split("=", 1)
                    # only if key not start with
                    if not key.startswith("addendum_id") and not key.startswith("hipsgen_"):
                        data[key.strip()] = value.strip()
                dist_ls.append(data)
            df = pd.DataFrame(dist_ls)
            df["moc_sky_fraction"] = df["moc_sky_fraction"].astype(float)
            _HiPS_metadata_response = df[
                (~df["client_category"].str.contains("solar", case=False, na=False))
                & (df["moc_sky_fraction"].astype(float) > 0.7)
            ].reset_index(drop=True)
        return _HiPS_metadata_response

    @classmethod
    def search_sky_background(cls, keywords: Optional[str] = None):
        """
        Search for HiPS background images based on keywords

        Parameters
        ----------
        keywords : str
            Keywords to search for, multiple keywords can be separated by space

        Returns
        -------
        list
            List of HiPS background images
        """
        df = cls.parse_hips_background()
        if keywords is None:
            return df["obs_title"].tolist()
        else:
            # replace some common spelling variations
            keywords = keywords.replace("colour", "color").replace("gray", "grey").replace("tmass", "2mass")
            # search for all keywords in the title and description, only return the corresponding title only if all keywords are found
            keywords = keywords.lower()
            # use \b to match the whole word
            return df["obs_title"][
                reduce(
                    # make sure all keywords are found in the title or description
                    lambda x, y: x & y,
                    [
                        df[["obs_title", "obs_description", "obs_regime"]]
                        .apply(
                            lambda x: x.str.contains(
                                f"\\b{kw}\\b", case=False, na=False
                            )
                        )
                        .any(axis=1)
                        for kw in keywords.split()
                    ],
                )
            ].tolist()

    def get_hips_images(self, hips_id: str):
        df = self.parse_hips_background()
        obs_title = df["obs_title"].tolist()
        obs_copyright = df["obs_copyright"].tolist()
        allowed_id = df["ID"].tolist()

        # need to use & operator to avoid short-circuiting
        if (cond1 := hips_id not in allowed_id) & (cond2 := hips_id not in obs_title):
            raise ValueError(
                f"Unknown HiPS ID `{hips_id}`, allowed IDs are: {allowed_id}"
            )
        if cond2:
            obs_copyright = obs_copyright[allowed_id.index(hips_id)]
        if cond1:
            obs_copyright = obs_copyright[obs_title.index(hips_id)]

        # Create a new WCS astropy object
        horizontal_pix = 2000
        vertical_pix = horizontal_pix // (self.radius[0].value / self.radius[1].value)
        w = astropy_wcs.WCS(
            header={
                "NAXIS1": horizontal_pix,  # Width of the output fits/image
                "NAXIS2": vertical_pix,  # Height of the output fits/image
                "WCSAXES": 2,  # Number of coordinate axes
                "CRPIX1": horizontal_pix / 2,  # Pixel coordinate of reference point
                "CRPIX2": vertical_pix / 2,  # Pixel coordinate of reference point
                "CDELT1": 2.0
                * self.radius[0].value
                / horizontal_pix,  # [deg] Coordinate increment at reference point
                "CDELT2": 2.0
                * self.radius[1].value
                / vertical_pix,  # [deg] Coordinate increment at reference point
                "CUNIT1": "deg",  # Units of coordinate increment and value
                "CUNIT2": "deg",  # Units of coordinate increment and value
                # https://docs.astropy.org/en/stable/wcs/supported_projections.html
                "CTYPE1": "GLON-CAR",  # galactic longitude
                "CTYPE2": "GLAT-CAR",  # galactic latitude
                "CRVAL1": -self.center[
                    0
                ].value,  # [deg] Coordinate value at reference point
                "CRVAL2": self.center[
                    1
                ].value,  # [deg] Coordinate value at reference point
            }
        )
        result_image = hips2fits.query_with_wcs(
            hips=hips_id,
            wcs=w,
            get_query_payload=False,
            format="jpg",
        )
        return np.flip(result_image, axis=1), obs_copyright

    def unit_check(self, x, unit):
        """
        Check if the unit is the same as the default unit and convert to the default unit
        """
        # force the unit to be the same as the default unit
        x <<= unit
        return x

    def xy_unit_check(self, x, y):
        if not isinstance(x, u.quantity.Quantity) or not isinstance(
            y, u.quantity.Quantity
        ):
            raise TypeError("Both x and y must carry astropy's unit")
        else:
            if x.unit is not None and y.unit is not None:
                x = x.to(self.unit).value
                y = y.to(self.unit).value
            else:
                raise TypeError(
                    "Both x, y, center and radius must carry astropy's unit"
                )

        return x, y

    @property
    def citation(self):
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
    angle : int
        Where the Sun is from the center of the galaxy in degrees. 0 is north, 90 is east...
    r0 : float
        Distance to galactic center in kpc
    coord : str
        'galactocentric' or 'galactic'
    center : tuple
        Coordinates as in the system `coord` of the center of the plot with astropy units
    radius : float
        Radius of the plot with astropy units
    unit : astropy units
        Astropy units for the default if none provided and the axes
    figsize : tuple
        Figure size
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
        super().__init__()
        self.grayscale = grayscale
        self.annotation = annotation
        if angle != 90:
            raise NotImplementedError("Only 90 degrees is implemented")
        self.angle = -(angle - 90.0)
        self.r0 = self.unit_check(r0, unit)
        self.coord = coord
        self.center = self.unit_check(center, unit)
        self.radius = self.unit_check(radius, unit)
        self.unit = unit
        self.figsize = figsize
        self._built = False
        self.facecolor = (0, 0, 0) if not self.grayscale else (1, 1, 1)

        # properties of the images
        self.img_pixels = 5600  # number of pixels in x and y axis
        # light years per pixel resolution given
        # 1078 is the number of pixels in the image from galactic center to the Sun
        self.img_resolution = self.unit_check((self.r0 / 1078), unit)

    def read_bg_img(self):
        if self.annotation:
            self.img_obj = self._MW_IMAGES["MW_bg_annotate"]
        else:
            self.img_obj = self._MW_IMAGES["MW_bg_unannotate"]
        self.reference_str = self.img_obj.citation
        self._gh_img_url = self._gh_imgbase_url + self.img_obj.filename

        if self.coord.lower() == "galactic":
            # shift the coord by r0 backt to the galactocentric
            x_shift = self.r0
            # example: if the center is 0,0 and the angle is 90, the Sun is at 0,r0
            self.center[0] += x_shift
            self.coord_english = "Galactic"
        elif self.coord.lower() == "galactocentric":
            x_shift = 0.0 * self.unit
            self.coord_english = "Galactocentric"
        else:
            raise ValueError(
                "Unknown coordinates, can only be `galactic` or `galactocentric`"
            )
        if self.angle % 90 != 0:  # in case the angle is not 90, 180, 270, 360
            self.coord_english = ""

        # calculate the pixel radius and center from physical units
        pixel_radius = int((self.radius / self.img_resolution))
        pixel_center = [  # where the center of the image is in pixel
            int((self.img_pixels / 2 + self.center[0] / self.img_resolution)),
            int((self.img_pixels / 2 - self.center[1] / self.img_resolution)),
        ]

        # get the pixel coordinates of four corners
        x_left_px = pixel_center[0] - pixel_radius
        x_right_px = pixel_center[0] + pixel_radius
        y_top_px = pixel_center[1] + pixel_radius
        y_bottom_px = pixel_center[1] - pixel_radius

        # decide whether boundary is within the image
        crop_arg = [0, 0, self.img_pixels, self.img_pixels]
        # image extent in physical coordinates
        img_ext_phy = [
            -self.img_pixels / 2 * self.img_resolution.value,
            self.img_pixels / 2 * self.img_resolution.value,
            -self.img_pixels / 2 * self.img_resolution.value,
            self.img_pixels / 2 * self.img_resolution.value,
        ]
        if x_left_px >= 0:
            crop_arg[0] = x_left_px
        if x_right_px <= self.img_pixels:
            crop_arg[2] = x_right_px
        if y_bottom_px >= 0:
            crop_arg[1] = y_bottom_px
        if y_top_px <= self.img_pixels:
            crop_arg[3] = y_top_px

        # crop the image before rotation
        self.bg_img = self.img_obj.pillow_img.rotate(
            self.angle, fillcolor=self.facecolor
        )

        # rotate the four corners
        rotation_matrix = np.array(
            [
                [np.cos(np.deg2rad(self.angle)), -np.sin(np.deg2rad(self.angle))],
                [np.sin(np.deg2rad(self.angle)), np.cos(np.deg2rad(self.angle))],
            ]
        )
        bottom_left_phy = np.dot(rotation_matrix, [img_ext_phy[0], img_ext_phy[2]])
        top_left_phy = np.dot(rotation_matrix, [img_ext_phy[0], img_ext_phy[3]])
        bottom_right_phy = np.dot(rotation_matrix, [img_ext_phy[1], img_ext_phy[2]])
        top_right_phy = np.dot(rotation_matrix, [img_ext_phy[1], img_ext_phy[3]])

        # extent=[horizontal_min,horizontal_max,vertical_min,vertical_max]
        self.bg_img_ext = [
            np.min(
                [
                    bottom_left_phy[0],
                    top_left_phy[0],
                    bottom_right_phy[0],
                    top_right_phy[0],
                ]
            )
            - x_shift.value,
            np.max(
                [
                    bottom_left_phy[0],
                    top_left_phy[0],
                    bottom_right_phy[0],
                    top_right_phy[0],
                ]
            )
            - x_shift.value,
            np.min(
                [
                    bottom_left_phy[1],
                    top_left_phy[1],
                    bottom_right_phy[1],
                    top_right_phy[1],
                ]
            ),
            np.max(
                [
                    bottom_left_phy[1],
                    top_left_phy[1],
                    bottom_right_phy[1],
                    top_right_phy[1],
                ]
            ),
        ]

        # actual extent of the image
        self._ext = [
            (self.center[0] - self.radius - x_shift).value,
            (self.center[0] + self.radius - x_shift).value,
            (self.center[1] - self.radius).value,
            (self.center[1] + self.radius).value,
        ]

        self.bg_img = np.asanyarray(self.bg_img)
        if self.grayscale:
            self.bg_img = rgb2gray(self.bg_img)

        self._aspect = (
            self.bg_img.shape[0]
            / float(self.bg_img.shape[1])
            * ((self._ext[1] - self._ext[0]) / (self._ext[3] - self._ext[2]))
        )
        self._aspect = np.abs(self._aspect)

        return None


class MWSkyMapBase(MWPlotCommon):
    """
    MWSkyMap base class to plot the sky map
    """

    def __init__(
        self,
        grayscale: bool = False,
        projection: str = "equirectangular",
        background: str = "optical",
        center: tuple = (0.0, 0.0),
        radius: tuple = (180.0, 90.0),
        figsize: tuple = (5, 5),
    ):
        super().__init__()
        self.projection = projection

        if self.projection in (
            allowed_proj := [
                "equirectangular",
                "aitoff",
                "hammer",
                "mollweide",
            ]
        ):
            pass
        # projection that is not suitable for sky map but available from matplotlib
        elif self.projection in ["lambert", "polar", "rectilinear"]:
            raise NotImplementedError(
                f"`{self.projection}` projection is not implemented for sky map"
            )
        else:
            raise ValueError(
                f"Unknown projection `{self.projection}`, allowed values are: {allowed_proj}"
            )

        self.wavlength = background

        # turn object name to coordinates
        if isinstance(center, str):
            c = self.objname_to_coord(center)
            c = c.transform_to(coord.Galactic)
            center = (-c.l.wrap_at(180 * u.degree).value, c.b.value) * u.deg

        self.center = center
        self.radius = radius
        self.grayscale = grayscale
        self.figsize = figsize
        self._built = False

        self._opposite_color = "white"
        if self.grayscale:
            self._opposite_color = "black"

    def read_bg_img(self):
        img_key = None
        need_copping = True
        if self.wavlength == "optical":
            img_key = "MW_edgeon_edr3_unannotate"
        elif self.wavlength == "gamma":
            img_key = "MW_fermi_gamma"
        elif self.wavlength == "infrared":
            img_key = "MW_2mass"
        elif self.wavlength == "far-infrared":
            img_key = "MW_farinfrared"
        else:
            self.bg_img, self.reference_str = self.get_hips_images(self.wavlength)
            img_obj = None
            need_copping = False
        if img_key:
            img_obj = self._MW_IMAGES[img_key]
            self.bg_img = plt.imread(img_obj.img_path)
            self.reference_str = img_obj.citation

        # find center pixel and radius pixel
        y_img_center = self.bg_img.shape[0] // 2 - int(
            (self.bg_img.shape[0] / 180) * self.center[1].value
        )
        y_radious_px = int((self.bg_img.shape[0] / 180) * self.radius[1].value)
        x_img_center = (
            int((self.bg_img.shape[1] / 360) * self.center[0].value)
            + self.bg_img.shape[0]
        )
        x_radious_px = int((self.bg_img.shape[1] / 360) * self.radius[0].value)

        self._ext = [
            (self.center[0] - self.radius[0]).value,
            (self.center[0] + self.radius[0]).value,
            (self.center[1] - self.radius[1]).value,
            (self.center[1] + self.radius[1]).value,
        ]

        if need_copping:
            self.bg_img = self.bg_img[
                (y_img_center - y_radious_px) : (y_img_center + y_radious_px),
                (x_img_center - x_radious_px) : (x_img_center + x_radious_px),
                :,
            ]

        if self.grayscale:
            self.bg_img = rgb2gray(self.bg_img)

        if img_obj:
            self._gh_img_url = self._gh_imgbase_url + img_obj.filename

        return None

    def radec_unit_check(self, ra, dec):
        if not isinstance(ra, u.quantity.Quantity) or not isinstance(
            dec, u.quantity.Quantity
        ):
            raise TypeError("Both RA and DEC must carry astropy's unit")
        else:
            if ra.unit is not None and dec.unit is not None:
                ra = ra.to(self.unit)
                dec = dec.to(self.unit)
                c_icrs = coord.SkyCoord(ra=ra, dec=dec, frame="icrs")
                if self.projection == "equirectangular":
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
