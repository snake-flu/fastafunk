from pkg_resources import get_distribution

try:
    __version__ = get_distribution("fastafunk").version
except:
    __version__ = "local"

__all__ = ["consensus", "extract", "merge", "remove", "split", "count", "subsample"]

from fastafunk import *
