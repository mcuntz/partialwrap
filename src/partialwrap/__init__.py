#!/usr/bin/env python
"""
Purpose
=======

partialwrap provides wrappers for Python functions and external executables so
that they can easily be partialised with functools.partial.

:copyright: Copyright 2016-2022 Matthias Cuntz, see AUTHORS.md for details.
:license: MIT License, see LICENSE for details.

Subpackages
===========
.. autosummary::
    wrappers
    std_io

History
    * Written Mar 2020 by Matthias Cuntz (mc (at) macu (dot) de)
    * v1.0, initial PyPI commit, Jun 2020, Matthias Cuntz
    * v1.0.1, trigger zenodo, Jun 2020, Matthias Cuntz
    * v1.1, make flake8 compliant, Dec 2020, Matthias Cuntz
    * v1.2, remove scipy.optimize from tests, Dec 2020, Matthias Cuntz
    * v1.3, allow non-numeric parameters, May 2021, Matthias Cuntz
    * v1.3.1, more general treatment of right-hand side in sub_params_names,
      May 2021, Matthias Cuntz
    * v1.3.2, protect saved groups in replacement pattern in std_io,
      May 2021, Matthias Cuntz
    * v1.3.3, deleted trailing print statement, May 2021, Matthias Cuntz
    * v1.3.4, escape backslash in params in std_io, May 2021, Matthias Cuntz
    * v1.3.5, escape backslash also in sub_params_ja, May 2021, Matthias Cuntz
    * v1.4, move to pyproject.toml setup with Github actions,
      Oct 2021, Matthias Cuntz
    * v1.5, Allow empty right-hand sides in files, such as `name1 =`,
      Nov 2021, Matthias Cuntz
    * v1.6, Use subprocess.run for Python > v3.6, Aug 2022, Matthias Cuntz
    * v1.7, Added error function in `exe_wrapper`, Aug 2022, Matthias Cuntz
    * v2.0, Major update of documentation and docstrings,
      Nov 2023, Matthias Cuntz

"""
# version, author
try:
    from ._version import __version__
except ImportError:  # pragma: nocover
    # package is not installed
    __version__ = "0.0.0.dev0"
__author__  = "Matthias Cuntz"

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


__all__ = ['exe_wrapper', 'exe_mask_wrapper',
           'function_wrapper', 'function_mask_wrapper',
           'sub_params_ja',
           'sub_params_names', 'sub_params_names_case',
           'sub_params_names_ignorecase',
           'standard_output_reader',
           'standard_parameter_reader', 'standard_parameter_writer',
           'standard_parameter_reader_bounds_mask',
           'standard_parameter_writer_bounds_mask',
           'standard_time_series_reader', 'standard_timeseries_reader']
