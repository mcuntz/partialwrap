#!/usr/bin/env python
"""
Purpose
=======

partialwrap provides wrappers for Python functions and external
executables so that they can easily be partialised with
functools.partial.

:copyright: Copyright 2016-2020 Matthias Cuntz, see AUTHORS.md for details.
:license: MIT License, see LICENSE for details.

Subpackages
===========
.. autosummary::
    wrappers
    std_io
    version
"""
from __future__ import division, absolute_import, print_function

from .version import __version__, __author__

# Function wrappers to be used with partial from functools
from .wrappers import exe_wrapper, exe_mask_wrapper
from .wrappers import function_wrapper, function_mask_wrapper

# Standard parameter reader and writer functions as well as output reader
# functions
from .std_io import sub_params_ja
from .std_io import sub_params_names, sub_params_names_case
from .std_io import sub_params_names_ignorecase
from .std_io import standard_output_reader
from .std_io import standard_parameter_reader, standard_parameter_writer
from .std_io import standard_parameter_reader_bounds_mask
from .std_io import standard_parameter_writer_bounds_mask
from .std_io import standard_time_series_reader, standard_timeseries_reader
