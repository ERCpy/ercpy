# -*- coding: utf-8 -*-
# Copyright 2014 by Forschungszentrum Juelich GmbH
# Author: J. Caron
#
"""# TODO: Global description (one line)!

Modules
-------
formats
    # TODO: Add description!
holography
    # TODO: Add description!
eelsedx
    # TODO: Add description!
utils
    # TODO: Add description!

"""


from .formats import *  # analysis:ignore
from .holography import *  # analysis:ignore
from .eelsedx import *  # analysis:ignore
from .mtools import *  # analysis:ignore
import utils
import config
from .version import version as __version__

import logging
_log = logging.getLogger(__name__)
_log.info("Starting ERCpy V{}".format(__version__))
del logging

__all__ = ['utils', 'config']
__all__.extend(formats.__all__)
__all__.extend(holography.__all__)
__all__.extend(eelsedx.__all__)
