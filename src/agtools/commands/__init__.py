"""agtools: Tools for manipulating assembly graphs"""

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Alpha"


from .stats import stats
from .rename import rename
from .merge import merge
from .fastg2gfa import fastg2gfa
from .gfa2fasta import gfa2fasta
from .gfa2fastg import gfa2fastg
