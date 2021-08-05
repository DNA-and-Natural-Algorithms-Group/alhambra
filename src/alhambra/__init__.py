__all__ = ["TileSet", "Tile", "TileList", "Glue", "GlueList"]

from . import classes

from .tilesets import TileSet
from .tiles import Tile, TileList
from .glues import Glue, GlueList
from importlib.metadata import PackageNotFoundError, version  # pragma: no cover

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = "alhambra"
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
