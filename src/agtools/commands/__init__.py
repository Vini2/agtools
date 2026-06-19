"""agtools: A Software Framework to Manipulate Assembly Graphs"""

from agtools import __version__

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = __version__
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


from .asqg2gfa import asqg2gfa
from .clean import clean
from .component import component
from .concat import concat
from .fastg2gfa import fastg2gfa
from .filter import filter
from .gfa2adj import gfa2adj
from .gfa2asqg import gfa2asqg
from .gfa2dot import gfa2dot
from .gfa2fasta import gfa2fasta
from .gfa2fastg import gfa2fastg
from .rename import rename
from .stats import stats

__all__ = [
    "__version__",
    "asqg2gfa",
    "clean",
    "component",
    "concat",
    "fastg2gfa",
    "filter",
    "gfa2adj",
    "gfa2asqg",
    "gfa2dot",
    "gfa2fasta",
    "gfa2fastg",
    "rename",
    "stats",
]
