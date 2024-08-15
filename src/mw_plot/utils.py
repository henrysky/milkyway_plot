from numpy.typing import NDArray
from typing import Tuple
from astropy.coordinates import SkyCoord, ICRS
import astropy.units as u
import numpy as np

# RA and DEC of galactic center, galactic anti-center, galactic north and south pole in degree
center_radec = [266.4167, -29.0078]
anti_center_radec = [86.4167, 28.0078]
northpole_radec = [192.7667, 27.1167]
southpole_radec = [12.7667, -27.1167]


def mw_radec(deg: bool = True, size: int = 3600) -> Tuple[NDArray, NDArray]:
    """
    Get RA DEC coordinates of the milkyway in the ICRS system

    Parameters
    ----------
    deg : bool, optional
        Return RA and DEC in degree if True else in radian, by default True
    size : int, optional
        Number of points to generate, by default 3600

    Returns
    -------
    Tuple
        Tuple of RA and DEC coordinates
    """
    c = SkyCoord(
        np.linspace(0.0, 360.0, size) * u.deg, np.zeros(size) * u.deg, frame="galactic"
    )
    c = c.transform_to(ICRS)
    idx = np.argsort(c.ra.to(u.deg).value)
    if deg:
        ang_unit = u.deg
    else:
        ang_unit = u.rad

    return c.ra.to(ang_unit).value[idx], c.dec.to(ang_unit).value[idx]


def rgb2gray(rgb: NDArray) -> NDArray:
    """
    Change 8-bit RGB color image [0...255] into grayscale in RGB 8-bit representation [0...255]
    using colorimetric (perceptual luminance-preserving) conversion

    Parameters
    ----------
    rgb : NDArray
        RGB 2D image array

    Returns
    -------
    NDArray
        Grayscale 2D image in RGB array
    """
    rgb_norm = rgb / 255.0
    # need to seperate the color channel to avoid using excessive memory
    rgb_norm[:, :, 0] = np.where(
        rgb_norm[:, :, 0] <= 0.04045,
        rgb_norm[:, :, 0] / 12.92,
        ((rgb_norm[:, :, 0] + 0.055) / 1.055) ** 2.4,
    )
    rgb_norm[:, :, 1] = np.where(
        rgb_norm[:, :, 1] <= 0.04045,
        rgb_norm[:, :, 1] / 12.92,
        ((rgb_norm[:, :, 1] + 0.055) / 1.055) ** 2.4,
    )
    rgb_norm[:, :, 2] = np.where(
        rgb_norm[:, :, 2] <= 0.04045,
        rgb_norm[:, :, 2] / 12.92,
        ((rgb_norm[:, :, 2] + 0.055) / 1.055) ** 2.4,
    )
    gray = np.rint(
        255.0
        * (
            1
            - (
                0.2126 * rgb_norm[:, :, 0]
                + 0.7152 * rgb_norm[:, :, 1]
                + 0.0722 * rgb_norm[:, :, 2]
            )
        )
    ).astype(np.uint8)
    return np.concatenate([gray[:, :, None]] * 3, axis=-1)
