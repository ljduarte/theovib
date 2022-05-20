
from .input import *
from .internal import *
from .molecule import Molecule
from .ir import *
from .ptable import * 
from .matrices import *

# read version from installed package
from importlib.metadata import version  
__version__ = version("theovib")

