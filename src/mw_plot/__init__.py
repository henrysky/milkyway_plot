from importlib.metadata import version
from importlib.util import find_spec

from mw_plot.matplotlib_backend import MWFaceOn, MWSkyMap
try:
    import bokeh
    _HAS_BOKEH = True
    from mw_plot.bokeh_backend import MWFaceOnBokeh, MWSkyMapBokeh
except ImportError:
    _HAS_BOKEH = False

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
    "MWFaceOn",
    "MWSkyMap",
    "MWFaceOnBokeh", 
    "MWSkyMapBokeh",
    "center_radec",
    "anti_center_radec",
    "northpole_radec",
    "southpole_radec",
    "mw_radec",
    "rgb2gray",
    "_HAS_BOKEH",
]

version = __version__ = version("mw_plot")
