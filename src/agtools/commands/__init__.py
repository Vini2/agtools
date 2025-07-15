"""agtools: Tools for manipulating assembly graphs"""

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Alpha"


from .component import component
from .fastg2gfa import fastg2gfa
from .filter import filter
from .gfa2fasta import gfa2fasta
from .gfa2fastg import gfa2fastg
from .merge import merge
from .rename import rename
from .stats import stats
