#!/usr/bin/env python
"""
Standard parameter reader and writer functions as well as output reader
functions, including substitution of tags #JA????# or names in files.

Parameter files can be written in some arbitrary standard formats.
Or a prepared parameter file can be read and parameter values replaced.

This module was written by Matthias Cuntz while at Institut National de
Recherche pour l'Agriculture, l'Alimentation et l'Environnement (INRAE), Nancy,
France.

Copyright (c) 2016-2020 Matthias Cuntz - mc (at) macu (dot) de

Released under the MIT License; see LICENSE file for details.

History:

* Written Nov 2016 by Matthias Cuntz (mc (at) macu (dot) de),
  including standard readers and writers
* Added substitution function for #JA????# in prepared parameter files
  (sub_ja_params_files), Jan 2018, Matthias Cuntz
* Use .pid as suffix for multiprocessing, Feb 2018, Matthias Cuntz
* Return ndarrays instead of lists and line IDs in standard_parameter_reader,
  Mar 2018, Matthias Cuntz
* Added substitution functions for arbitrary names in parameter files
  (sub_names_params_files*), Mar 2018, Matthias Cuntz
* Allow substitution in several parameter files (_msub_files), Apr 2018,
  Matthias Cuntz
* Changed to Sphinx docstring and numpydoc, Dec 2019, Matthias Cuntz
* Renamed standard_parameter_reader/_writer to
  standard_parameter_reader_bounds_mask/standard_parameter_writer_bounds_mask,
  Jan 2020, Matthias Cuntz
* New standard_parameter_reader/_writer that read/write single values
  per parameter, Jan 2020, Matthias Cuntz
* Swapped names and params in call to sub_names_params_files* to be compatible
  with new generic exe_wrapper, Jan 2020, Matthias Cuntz
* Call standard_parameter_writer with 2 or 3 arguments, i.e. pid given or not,
  Jan 2020, Matthias Cuntz
* Make all strings raw strings in sub_names_params_files_* routines to deal
  with regular expressions, Jan 2020, Matthias Cuntz
* Keep formatting of names and spaces with sub_names_params functions;
  close input file before raising error, Feb 2020, Matthias Cuntz
* Prepend underscore to local function, Jun 2020, Matthias Cuntz
* Use with statement for reading file so that they are properly closed,
  Jun 2020, Matthias Cuntz
* Overwrite input files by default if no pid given, Jun 2020, Matthias Cuntz
* Remove _files from function names, Jun 2020, Matthias Cuntz
* pid keyword after arguments, Jun 2020, Matthias Cuntz
* Rename all substitution functions sub_* to start with sub_params_*, Jun 2020,
  Matthias Cuntz
* Add pid keyword to all routines, Jun 2020, Matthias Cuntz
* Make flake8 compliant, Dec 2020, Matthias Cuntz
* Allow non-numeric parameters, Feb 2021, Matthias Cuntz
* More generic right-hand in sub_params_names, May 2021, Matthias Cuntz
* Escape backslash in names in sub_params_names, May 2021, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided:

.. autosummary::
   sub_params_ja
   sub_params_names
   sub_params_names_case
   sub_params_names_ignorecase
   standard_output_reader
   standard_parameter_reader
   standard_parameter_writer
   standard_parameter_reader_bounds_mask
   standard_parameter_writer_bounds_mask
   standard_timeseries_reader
   standard_time_series_reader
"""
from __future__ import division, absolute_import, print_function
import re
import numpy as np


__all__ = ['sub_params_ja',
           'sub_params_names', 'sub_params_names_case',
           'sub_params_names_ignorecase',
           'standard_output_reader',
           'standard_parameter_reader', 'standard_parameter_writer',
           'standard_parameter_reader_bounds_mask',
           'standard_parameter_writer_bounds_mask',
           'standard_time_series_reader', 'standard_timeseries_reader']


# ------------------------------------------------------------------------------


def _msub(dic, text, flags=0):
    """
    Helper function for substituting several patterns in one string.

    pattern/replacement are given as dictionary: d[pattern] = replacement

    Parameters
    ----------
    dic : dict
        Pattern/replacement dictionary: `dic[pattern] = replacement`
    text : string
        String on which patterns will be replaced
    flags : int, optional
        Flags will be passed to `re.sub()` (default: 0)

    Returns
    -------
    text : string
        Input string text with replaced patterns

    Notes
    -----
    Compiled version as in
        http://code.activestate.com/recipes/81330-single-pass-multiple-replace/
    and
        https://www.safaribooksonline.com/library/view/python-cookbook-2nd/0596007973/ch01s19.h
    do not work because match.string returns the string and not the pattern so
    that the key for the dictionary does not work anymore.

    History
    -------
    Written,  Matthias Cuntz, Mar 2018
    Modified, Matthias Cuntz, Dec 2019 - Sphinx docstring
    """
    for d in dic:
        text = re.sub(d, dic[d], text, flags=flags)
    return text


# ------------------------------------------------------------------------------


def _msub_files(files, dic, pid=None, flags=0):
    """
    Helper function for applying replacement dictionary on several files.
    pattern/replacement are given as dictionary: d[pattern] = replacement

    Parameters
    ----------
    files : list
        List with file names in which pattern replacement will be applied
        on each line.
    dic : dict
        Pattern/replacement dictionary: dic[pattern] = replacement
    pid : int, optional
        Output files will be input files suffixed by .pid
        (default: overwrite input file)
    flags : int, optional
        Flags will be passed to re.sub() (default: 0)

    Returns
    -------
    files :
        No return value but either changed input files or output files with
        names of the input files suffixed by .pid, in which all occurences of
        all patterns were replaced.

    History
    -------
    Written,  Matthias Cuntz, Apr 2018
    Modified, Matthias Cuntz, Dec 2019 - Sphinx docstring
              Matthias Cuntz, Jun 2020 - use with statement for reading files
              Matthias Cuntz, Jun 2020 - overwrite input files by default
                                         if pid not given
    """
    from os.path import exists

    for f in files:
        if not exists(f):
            raise IOError('File does not exist: '+f)
        with open(f, 'r') as fi:
            tt = fi.read()
        tt = _msub(dic, tt, flags=flags)
        if pid:
            fname = f+'.'+str(pid)
        else:
            fname = f
        with open(fname, 'w') as ff:
            ff.write(tt)

    return


# ------------------------------------------------------------------------------


def sub_params_ja(files, params, pid=None):
    """
    Substitute #JA????# with parameter value in several files, i.e.

        #JA0000# with params[0]

        #JA0001# with params[1]

        ...

    Parameters
    ----------
    files : list
        List with file names in which #JA????# will be replaced.

    params : iterable
        Parameter values to replace #JA????# patterns.

        params[0] will replace #JA0000#

        params[1] will replace #JA0001#

        ...

    pid : int, optional
        Output files will be input files suffixed by .pid
        (default: overwrite input file)

        Note that the default `pid=None` might not be useful if you have
        carefully prepared input file(s).

    Returns
    -------
    None
        No return value but either changed input files or output files with
        names of the input files suffixed by .pid, in which all #JA????#
        patterns were replaced by params elements.

    Examples
    --------
    >>> sub_params_ja([file1, file2], [0, 1, 2, 3], pid=1234)


    History
    -------
    Written,  Matthias Cuntz, Jan 2018
    Modified, Matthias Cuntz, Feb 2018 - pid
              Matthias Cuntz, Mar 2018 - use _msub
              Matthias Cuntz, Apr 2018 - use _msub_files
              Matthias Cuntz, Dec 2019 - Sphinx docstring
              Matthias Cuntz, Jun 2020 - overwrite input files by default
                                         if pid not given
              Matthias Cuntz, Jun 2020 - pid keyword after arguments
              Matthias Cuntz, Feb 2021 - allow non-numeric parameter
    """
    # assert list of files
    if isinstance(files, str):
        files = [files]

    # make dict for replacement
    dd = {}
    for i, p in enumerate(params):
        k = "#JA{:04d}#".format(i)
        try:
            repl = "{:.14e}".format(p)
        except ValueError:
            repl = "{}".format(p)
        dd[k] = repl

    # replace in each file
    _msub_files(files, dd, pid)

    return


# ------------------------------------------------------------------------------


def sub_params_names_case(files, params, names, pid=None):
    """
    Substitute `name = .*` with `name = parameter` value in several files, i.e.

        `names[i] = params[i]`

    Note, `names` are case sensitive.

    Parameters
    ----------
    files : list
        List with file in which values of `names` will be replaced with
        `params`.
    params : iterable
        Parameter values to be given to variables on the right of = sign

        Variable in names[0] will be assigned value in params[0]

        Variable in names[1] will be assigned value in params[1]

        ...

    names : iterable
        Variable names on left of = sign in files
    pid : int, optional
        Output files will be input files suffixed by .pid
        (default: overwrite input file)

    Returns
    -------
    None
        No return value but either changed input files or output files with
        names of the input files suffixed by .pid, in which all variables given
        in names are assigned the values given in params.

    Examples
    --------
    >>> sub_params_names_case([file1, file2], [0, 1, 2, 3],
    ...                       ['param1', 'param2', 'param3', 'param4'],
    ...                       pid=1234)


    History
    -------
    Written,  Matthias Cuntz, Mar 2018
    Modified, Matthias Cuntz, Apr 2018 - use _msub_files
              Matthias Cuntz, Dec 2019 - Sphinx docstring
              Matthias Cuntz, Jan 2020 - swap names and params in argument list
              Matthias Cuntz, Jan 2020 - make all raw strings for regular
                                         expressions
              Matthias Cuntz, Feb 2020 - keep formatting of names and spaces
              Matthias Cuntz, Jun 2020 - pid keyword after arguments
              Matthias Cuntz, Feb 2021 - allow non-numeric parameter
              Matthias Cuntz, May 2021 - respect space after = sign
                                       - replace everything after space after
                                         = sign up to next space
              Matthias Cuntz, May 2021 - protect saved group with \g<> in
                                         replacement pattern
              Matthias Cuntz, May 2021 - escape backslash in names
    """
    # assert list of files
    if isinstance(files, str):
        files = [files]

    # make dict for _msub with dict[pattern] = replacement
    dd = {}
    for i, p in enumerate(params):
        try:
            repl = r"\g<1>\g<2>=\g<3>{:.14e}".format(p)
        except ValueError:
            p = p.replace('\\', '\\\\')
            repl = r"\g<1>\g<2>=\g<3>{}".format(p)
        nep = r"(" + names[i] + r"\s*)=(\s*).*"  # name = value
        k = r"^(\s*)" + nep                      # beginning of line
        dd[k] = repl                             # replacement using
                                                 # substitutions \1, \2, and \3
        k = r"(\n+\s*)" + nep                    # after newline
        dd[k] = repl
        print(k, dd[k])

    # replace in each file
    _msub_files(files, dd, pid)

    return


def sub_params_names_ignorecase(files, params, names, pid=None):
    """
    Substitute `name = .*` with `name = parameter` value in several files, i.e.

        `names[i] = params[i]`

    Note, `names` are case insensitive.

    Parameters
    ----------
    files : list
        List with file in which values of `names` will be replaced with
        `params`.
    params : iterable
        Parameter values to be given to variables on the right of = sign

        Variable in names[0] will be assigned value in params[0]

        Variable in names[1] will be assigned value in params[1]

        ...

    names : iterable
        Variable names on left of = sign in files
    pid : int, optional
        Output files will be input files suffixed by .pid
        (default: overwrite input file)

    Returns
    -------
    None
        No return value but either changed input files or output files with
        names of the input files suffixed by .pid, in which all variables
        given in names are assigned the values given in params.

    Examples
    --------
    >>> sub_params_names_ignorecase([file1, file2], [0, 1, 2, 3],
    ...                             ['param1', 'param2', 'param3', 'param4'],
    ...                             pid=1234)


    History
    -------
    Written,  Matthias Cuntz, Mar 2018
    Modified, Matthias Cuntz, Apr 2018 - use _msub_files
              Matthias Cuntz, Dec 2019 - Sphinx docstring
              Matthias Cuntz, Jan 2020 - swap names and params in argument list
              Matthias Cuntz, Jan 2020 - make all raw strings for regular
                                         expressions
              Matthias Cuntz, Feb 2020 - keep formatting of names and spaces
              Matthias Cuntz, Jun 2020 - pid keyword after arguments
              Matthias Cuntz, Feb 2021 - allow non-numeric parameter
              Matthias Cuntz, May 2021 - respect space after = sign
                                       - replace everything after space after
                                         = sign up to next space
              Matthias Cuntz, May 2021 - protect saved group with \g<> in
                                         replacement pattern
              Matthias Cuntz, May 2021 - escape backslash in names
    """
    # assert list of files
    if isinstance(files, str):
        files = [files]

    # make dict for _msub with dict[pattern] = replacement
    dd = {}
    for i, p in enumerate(params):
        try:
            repl = r"\g<1>\g<2>=\g<3>{:.14e}".format(p)
        except ValueError:
            p = p.replace('\\', '\\\\')
            repl = r"\g<1>\g<2>=\g<3>{}".format(p)
        nep = r"("+names[i]+r"\s*)=(\s*)\S*"  # name = value
        k = r"^(\s*)"+nep                     # beginning of line
        dd[k] = repl                          # replacement using
                                              # substitutions \1, \2, and \3
        k = r"(\n+\s*)"+nep                   # after newline
        dd[k] = repl

    # replace in each file
    _msub_files(files, dd, pid, flags=re.I)

    return


def sub_params_names(*args, **kwargs):
    """
    Wrapper for :any:`sub_params_names_ignorecase`.
    """
    return sub_params_names_ignorecase(*args, **kwargs)


# ------------------------------------------------------------------------------


def standard_output_reader(filename, pid=None):
    """
    Standard output reader.

    The standard output reader (if outputreader=None) reads a single value
    from a file without header, comment line or similar.

    That means for example:

        0.0123456789e-02

    Parameters
    ----------
    filename : string
        Filename of file with output value
    pid : int, optional
        If present, `filename` will be suffixed by .pid

    Returns
    -------
    float
        Single number read from filename

    Examples
    --------
    >>> subprocess.call(model)
    >>> obj = standard_output_reader(filename, pid=1234)


    History
    -------
    Written,  Matthias Cuntz, Nov 2016
    Modified, Matthias Cuntz, Dec 2019 - Sphinx docstring
              Matthias Cuntz, Jun 2020 - pid keyword
    """
    # read output value
    if pid:
        fname = filename+'.'+str(pid)
    else:
        fname = filename
    with open(fname, 'r') as ff:
        obj = ff.readline()

    # return float
    return np.float(obj)


# ------------------------------------------------------------------------------


def standard_parameter_reader(filename, pid=None):
    """
    Read standard parameter file.

    The standard parameter file is a file containing
    1 line per parameter with the parameter value:

    Lines starting with # will be excluded.

    That means a standard parameter file might look like:

        #par

        3.000000000000000e-01

        2.300000000000000e-01

        1.440000000000000e+01

        3.000000000000000e-01

        ...

    Parameters
    ----------
    filename : string
        Filename with parameter values
    pid : int, optional
        If present, `filename` will be suffixed by .pid

    Returns
    -------
    ndarray
        Parameter values

    Examples
    --------
    >>> params = standard_parameter_reader(paramfile)


    History
    -------
    Written,  Matthias Cuntz, Jan 2020
    Modified, Matthias Cuntz, Jun 2020 - pid keyword
    """
    if pid:
        fname = filename+'.'+str(pid)
    else:
        fname = filename
    params = []
    with open(fname, 'r') as ff:
        for line in ff:
            if line.startswith('#'):
                continue
            params.append(np.float(line.strip()))

    return np.array(params, dtype=np.float)


# ------------------------------------------------------------------------------


def standard_parameter_writer(filename, params, pid=None):
    """
    Standard parameter writer.

    The standard parameter writer writes a file containing
    1 line per parameter with the parameter value.

    All values will be written in IEEE double precision: {:.14e}.

    That means:

        3.000000000000000e-01

        2.300000000000000e-01

        1.440000000000000e+01

        3.000000000000000e-01

        ...

    Parameters
    ----------
    filename : string
        Output filename with parameter values
    params : iterable
        Parameter values
    pid : int, optional
        If given, output file will be `filename` suffixed by .pid

    Returns
    -------
    None
        No return value but output file written: filename or filename.pid

    Examples
    --------
    >>> randst = np.random.RandomState()
    >>> pid = str(randst.randint(2147483647))
    >>> params = sample_parameter(pis, pmin, pmax, pmask)
    >>> standard_parameter_writer(paramfile, params, pid)


    History
    -------
    Written,  Matthias Cuntz, Jan 2020
    Modified, Matthias Cuntz, Jan 2020 - call with 2 or 3 arguments,
                                         i.e. pid given or not
              Matthias Cuntz, Jun 2020 - pid keyword after arguments
              Matthias Cuntz, Feb 2021 - allow non-numeric parameter
    """
    # Existing file will be overwritten
    if pid:
        ofile = filename+'.'+str(pid)
    else:
        ofile = filename
    with open(ofile, 'w') as ff:
        for pp in params:
            try:
                repl = '{:.14e}'.format(pp)
            except ValueError:
                repl = '{}'.format(pp)
            dstr = repl
            print(dstr, file=ff)

    return


# ------------------------------------------------------------------------------


def standard_parameter_reader_bounds_mask(filename, pid=None):
    """
    Read standard parameter file with parameter bounds and mask.

    The standard parameter file is a space separated file containing
    1 line per parameter with the following columns:

    identifier, current parameter value, minimum parameter value,
    maximum parameter value, parameter mask (1: include, 0: exclude).

    Lines starting with # will be excluded.

    That means a standard parameter file might look like:

        # value min max mask

        1 3.000000000000000e-01 0.000000000000000e+00 1.000000000000000e+00 1

        2 2.300000000000000e-01 -1.000000000000000e+00 1.000000000000000e+00 1

        3 1.440000000000000e+01 9.000000000000000e+00 2.000000000000000e+01 1

        4 3.000000000000000e-01 0.000000000000000e+00 1.000000000000000e+00 0

        ...

    Parameters
    ----------
    filename : string
        Filename with parameter values
    pid : int, optional
        If present, `filename` will be suffixed by .pid

    Returns
    -------
    list
        List with ndarrays of

            ids - identifier

            params - parameter values

            pmin - minimum parameter value

            pmax - maximum parameter value

            mask - parameter mask (1: include, 0: exclude from optimisation)

    Examples
    --------
    >>> ids, params, pmin, pmax, pmask = \
    ...     standard_parameter_reader_bounds_mask(paramfile)


    History
    -------
    Written,  Matthias Cuntz, Nov 2016
    Modified, Matthias Cuntz, Mar 2018 - return ids
                                       - return numpy.arrays
              Matthias Cuntz, Dec 2019 - Sphinx docstring
              Matthias Cuntz, Jan 2020 - renamed from standard_parameter_reader
                  to standard_parameter_reader_bounds_mask
              Matthias Cuntz, Feb 2020 - close file before raising error
              Matthias Cuntz, Jun 2020 - pid keyword
    """
    ids    = []
    params = []
    pmin   = []
    pmax   = []
    pmask  = []
    if pid:
        pfile = filename+'.'+str(pid)
    else:
        pfile = filename
    with open(pfile, 'r') as ff:
        for line in ff:
            ll = line.strip()
            if ll.startswith('#'):
                continue
            el = ll.split()
            if len(el) != 5:
                raise IOError('Line has no 5 columns for parameter: '+line)
            ids.append(el[0])
            params.append(np.float(el[1]))
            pmin.append(np.float(el[2]))
            pmax.append(np.float(el[3]))
            pmask.append(int(el[4]))

    return [ids,
            np.array(params, dtype=np.float),
            np.array(pmin, dtype=np.float),
            np.array(pmax, dtype=np.float),
            np.array(pmask, dtype=np.bool)]


# ------------------------------------------------------------------------------


def standard_parameter_writer_bounds_mask(filename, params, pmin, pmax, mask,
                                          pid=None):
    """
    Standard parameter writer with parameter bounds and mask.

    The standard parameter writer writes a space separated file containing
    1 header line (# value min max mask) plus 1 line per parameter with the
    following columns:

    consecutive parameter number, current parameter value,
    minimum parameter value, maximum parameter value,
    parameter mask (1: include, 0: exclude).

    All values will be written in IEEE double precision: {:.14e}.

    That means:

        # value min max mask

        1 3.000000000000000e-01 0.000000000000000e+00 1.000000000000000e+00 1

        2 2.300000000000000e-01 -1.000000000000000e+00 1.000000000000000e+00 1

        3 1.440000000000000e+01 9.000000000000000e+00 2.000000000000000e+01 1

        4 3.000000000000000e-01 0.000000000000000e+00 1.000000000000000e+00 0

        ...

    Parameters
    ----------
    filename : string
        Output filename with parameter values
    params : iterable
        Parameter values
    pmin : iterable
        Minimum parameter values
    pmax : iterable
        Maximum parameter values
    mask : iterable
       Parameter mask (1: include, 0: exclude from optimisation)
    pid : int, optional
        If given, output file will be `filename` suffixed by .pid

    Returns
    -------
    None
        No return value but output file written: filename or filename.pid

    Examples
    --------
    >>> randst = np.random.RandomState()
    >>> pid = str(randst.randint(2147483647))
    >>> params = sample_parameter(pis, pmin, pmax, pmask)
    >>> standard_parameter_writer_bounds_mask(paramfile, params, pmin, pmax,
    ...                                       pmask, pid)


    History
    -------
    Written,  Matthias Cuntz, Nov 2016
    Modified, Matthias Cuntz, Feb 2018 - pid
              Matthias Cuntz, Dec 2019 - Sphinx docstring
              Matthias Cuntz, Jan 2020 - renamed from standard_parameter_writer
                  to standard_parameter_writer_bounds_mask
                                       - no .pid to filename if pid is None
              Matthias Cuntz, Jun 2020 - pid keyword after arguments
    """
    # Assert correct call
    astr = 'Parameter and minima do not have the same length.'
    assert len(params) == len(pmin), astr
    astr = 'Parameter and maxima do not have the same length.'
    assert len(params) == len(pmax), astr
    astr = 'Parameter and mask do not have the same length.'
    assert len(params) == len(mask), astr

    # Convert mask to integer if boolean
    pmask = [ int(i) for i in mask ]

    # Existing file will be overwritten
    if pid:
        fname = filename+'.'+str(pid)
    else:
        fname = filename
    with open(fname, 'w') as ff:
        # header
        hstr = '# value min max mask'
        print(hstr, file=ff)
        # data
        for i in range(len(params)):
            dstr = '{:d} {:.14e} {:.14e} {:.14e} {:d}'.format(
                i+1, params[i], pmin[i], pmax[i], pmask[i])
            print(dstr, file=ff)

    return


# ------------------------------------------------------------------------------


def standard_time_series_reader(filename, pid=None):
    """
    Standard reader for time series.

    The standard time series reader reads a time series of arbitrary length
    from a file without header, comment line or similar.

    That means for example:

        0.0123456789e-02

        0.1234567890e-02

        0.2345678900e-02

        0.3456789000e-02

        0.4567890000e-02

        ...

    Parameters
    ----------
    filename : string
        Filename of with time series values
    pid : int, optional
        If present, `filename` will be suffixed by .pid

    Returns
    -------
    timeseries : ndarray
        ndarray with values of each line in filename

    Examples
    --------
    >>> subprocess.call(model)
    >>> ts = standard_time_series_reader(filename)


    History
    -------
    Written,  Matthias Cuntz, Jan 2018
    Modified, Matthias Cuntz, Dec 2019 - Sphinx docstring
              Matthias Cuntz, Jun 2020 - pid keyword
    """
    if pid:
        fname = filename+'.'+str(pid)
    else:
        fname = filename
    # read output value
    with open(fname, 'r') as ff:
        serie = ff.readlines()

    # return float array
    return np.array(serie, dtype=np.float)


def standard_timeseries_reader(*args, **kwargs):
    """
    Wrapper for :any:`standard_time_series_reader`
    """
    return standard_time_series_reader(*args, **kwargs)


# ------------------------------------------------------------------------------


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
