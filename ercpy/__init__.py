# -*- coding: utf-8 -*-
# Copyright 2014 by Forschungszentrum Juelich GmbH
# Author: J. Caron
#
"""# TODO: Global description (one line)!

Modules
-------
emd
    # TODO: Add description!
holography
    # TODO: Add description!
eelsedx
    # TODO: Add description!
utils
    # TODO: Add description!

"""


from .emd import *  # analysis:ignore
import holography
import eelsedx
import utils
import config
from .version import version as __version__

import logging
_log = logging.getLogger(__name__)
_log.info("Starting PYRAMID V{}".format(__version__))
del logging

__all__ = ['emd', 'holography', 'eelsedx', 'utils', 'config']
__all__.extend(emd.__all__)
