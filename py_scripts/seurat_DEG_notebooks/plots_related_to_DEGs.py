######### PLOTTING DEGs ##############

# pylint: disable=all

import os
import struct
import sys
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd
import scipy.stats as stats
from memory_profiler import profile
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA

import warnings

# add the parent directory of 'src' to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

############ fold change plots code ############
