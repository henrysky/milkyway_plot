from importlib.metadata import version
from importlib.util import find_spec

from mw_plot.mw_plot_matplotlib import MWFaceOn, MWSkyMap
from mw_plot.utils import (
    anti_center_radec,
    center_radec,
    mw_radec,
    northpole_radec,
    rgb2gray,
    southpole_radec,
)

if find_spec("bokeh") is not None:
    _HAS_BOKEH = True
else:
    _HAS_BOKEH = False


__all__ = [
    "MWPlot",
    "MWSkyMap",
    "center_radec",
    "anti_center_radec",
    "northpole_radec",
    "southpole_radec",
    "mw_radec",
    "rgb2gray",
    "_HAS_BOKEH",
]

version = __version__ = version("mw_plot")
