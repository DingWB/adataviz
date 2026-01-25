#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
# from .clustermap import (heatmap, ClusterMapPlotter, composite,
#                          DendrogramPlotter)
# from .oncoPrint import oncoprint, oncoPrintPlotter
# __all__=['*']
import sys
from ._version import version as __version__
from . import tools as tl
from . import preprocessing as pp
from . import plotting as pl
from loguru import logger as logger

_ROOT = os.path.abspath(os.path.dirname(__file__))
logger.remove()
logger.add(sys.stderr, level="DEBUG")