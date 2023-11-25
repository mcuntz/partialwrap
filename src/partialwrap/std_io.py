#!/usr/bin/env python
"""
Parameter writer and output reader functions

This module was written by Matthias Cuntz while at Institut National
de Recherche pour l'Agriculture, l'Alimentation et l'Environnement
(INRAE), Nancy, France.

Copyright (c) 2016-2023 Matthias Cuntz - mc (at) macu (dot) de

Released under the MIT License; see LICENSE file for details.

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

History
    * Written Nov 2016 by Matthias Cuntz (mc (at) macu (dot) de),
      including standard readers and writers
    * Added substitution function for #JA????# in prepared parameter files
      (sub_ja_params_files), Jan 2018, Matthias Cuntz
    * Use .pid as suffix for multiprocessing, Feb 2018, Matthias Cuntz
    * Return ndarrays instead of lists and line IDs in
      standard_parameter_reader, Mar 2018, Matthias Cuntz
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
    * Swapped names and params in call to sub_names_params_files* to be
      compatible with new generic exe_wrapper, Jan 2020, Matthias Cuntz
    * Call standard_parameter_writer with 2 or 3 arguments, i.e. pid given or
      not, Jan 2020, Matthias Cuntz
    * Make all strings raw strings in sub_names_params_files_* routines to deal
      with regular expressions, Jan 2020, Matthias Cuntz
    * Keep formatting of names and spaces with sub_names_params functions;
      close input file before raising error, Feb 2020, Matthias Cuntz
    * Prepend underscore to local function, Jun 2020, Matthias Cuntz
    * Use with statement for reading file so that they are properly closed,
      Jun 2020, Matthias Cuntz
    * Overwrite input files by default if no pid given,
      Jun 2020, Matthias Cuntz
    * Remove _files from function names, Jun 2020, Matthias Cuntz
    * pid keyword after arguments, Jun 2020, Matthias Cuntz
    * Rename all substitution functions sub_* to start with sub_params_*,
      Jun 2020, Matthias Cuntz
    * Add pid keyword to all routines, Jun 2020, Matthias Cuntz
    * Make flake8 compliant, Dec 2020, Matthias Cuntz
    * Allow non-numeric parameters, Feb 2021, Matthias Cuntz
    * More generic right-hand side in sub_params_names,
      May 2021, Matthias Cuntz
    * Escape backslash in params in sub_params routines,
      May 2021, Matthias Cuntz
    * np.float-> float, etc., Sep 2021, Matthias Cuntz
    * Use specific whitespace characters instead of all whitespace characters
      on right-hand side because the latter includes line breaks,
      Nov 2021, Matthias Cuntz
    * Helper function for constructing dict for `_msub`,
      Nov 2021, Matthias Cuntz
    * Added examples and reformatted docstrings, Aug 2022, Matthias Cuntz
    * Use formatted strings where appropriate, Nov 2023, Matthias Cuntz

"""
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


def _dict4msub(names, params):
    """
    Helper function for constructing the dict for `_msub`.

    pattern/replacement are given as dictionary: d[pattern] = replacement

    Parameters
    ----------
    names : iterable
        Variable names on left of equal (=) sign
    params : iterable
        Parameter values to be given to variables on the right of the equal (=)
        sign.
        Variable in *names[0]* will be assigned value in *params[0]*.
        Variable in *names[1]* will be assigned value in *params[1]*.
        etc.

    Returns
    -------
    dic : dict
        Pattern/replacement dictionary: `dic[pattern] = replacement`

    """
    dd = {}
    for i, p in enumerate(params):
        try:
            repl = r"\g<1>\g<2>=\g<3>{:.14e}".format(p)
        except ValueError:
            p1 = p.replace('\\', '\\\\')
            repl = r"\g<1>\g<2>=\g<3>{}".format(p1)
        nep = r"(" + names[i] + r"\s*)=([ \t\f\v]*).*"  # name = value
        k = r"^(\s*)" + nep     # beginning of line
        dd[k] = repl            # replacement using substitutions \\1, \\2, ...
        k = r"(\n+\s*)" + nep   # after newline
        dd[k] = repl

    return dd


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
    Compiled versions as in
        http://code.activestate.com/recipes/81330-single-pass-multiple-replace/
    and
        https://www.safaribooksonline.com/library/view/python-cookbook-2nd/0596007973/ch01s19.h
    do not work because match.string returns the string and not the pattern so
    that the key for the dictionary does not work anymore.

    """
    for d in dic:
        text = re.sub(d, dic[d], text, flags=flags)
    return text


# ------------------------------------------------------------------------------


def _msub_files(files, dic, pid=None, flags=0):
    """
    Helper function for applying replacement dictionary on several files

    pattern/replacement are given as dictionary: d[pattern] = replacement

    Parameters
    ----------
    files : list
        List with file names in which pattern replacement will be applied
        on each line.
    dic : dict
        Pattern/replacement dictionary: `dic[pattern] = replacement`
    pid : int, optional
        Output files will be input files suffixed by *.pid*
        (default: overwrite input file)
    flags : int, optional
        Flags will be passed to *re.sub()* (default: 0)

    Returns
    -------
    files
        No return value but either changed input files or output files with
        names of the input files suffixed by *.pid*, in which all occurences of
        all patterns were replaced.

    """
    from os.path import exists

    for f in files:
        if not exists(f):
            raise IOError('File does not exist: ' + f)
        with open(f, 'r') as fi:
            tt = fi.read()
        tt = _msub(dic, tt, flags=flags)
        if pid:
            fname = f'{f}.{pid}'
        else:
            fname = f
        with open(fname, 'w') as ff:
            ff.write(tt)

    return


# ------------------------------------------------------------------------------


def sub_params_ja(files, params, pid=None):
    """
    Substitute #JA????# with parameter values in several files

    This means replacing

       | #JA0000# with params[0]
       | #JA0001# with params[1]
       | ...

    Parameters
    ----------
    files : list
        List with file names in which #JA????# will be replaced.
    params : iterable
        | Parameter values to replace #JA????# patterns.
        | params[0] will replace #JA0000#
        | params[1] will replace #JA0001#
        | ...
    pid : int, optional
        Output files will be input files suffixed by *.pid*
        (default: overwrite input file)

        Note that the default *pid=None* might not be useful if you have
        carefully prepared input files.

    Returns
    -------
    None
        No return value but either changed input files or output files with
        names of the input files suffixed by *.pid*, in which all #JA????#
        patterns were replaced by *params*.

    Examples
    --------
    Replaces #JA0000# with 0, #JA0001# with 1, #JA0002# with 2, and
    #JA0003# with 3 in the two input files 'file1.txt' and 'file2.txt',
    writing the two output files 'file1.txt.1234' and 'file2.txt.1234':

    >>> sub_params_ja(['file1.txt', 'file1.txt'], [0, 1, 2, 3], pid=1234)

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
            p1 = p.replace('\\', '\\\\')
            repl = "{}".format(p1)
        dd[k] = repl

    # replace in each file
    _msub_files(files, dd, pid)

    return


# ------------------------------------------------------------------------------


def sub_params_names_case(files, params, names, pid=None):
    """
    Substitute `name = .*` with `name = parameter` value in several files

    This means writing `for i in range(len(params))`:

       `names[i] = params[i]`

    Note, *names* are case sensitive.

    Parameters
    ----------
    files : list
        List with files in which values of *names* will be replaced with
        *params*.
    params : iterable
        | Parameter values to be given to variables on the right of the
          equal (=) sign.
        | Variable in *names[0]* will be assigned value in *params[0]*
        | Variable in *names[1]* will be assigned value in *params[1]*
        | ...
    names : iterable
        Variable names on left of the equal (=) sign in *files*
    pid : int, optional
        Output files will be input files suffixed by *.pid*
        (default: overwrite input file)

    Returns
    -------
    None
        No return value but either changed input files or output files with
        names of the input files suffixed by *.pid*, in which all variables
        given in *names* are assigned the values given in *params*.

    Examples
    --------
    Replaces `param1 = ...` with `param1 = 0`, `param2 = ...` with
    `param2 = 1`, `param3 = ...` with `param3 = 2`, and `param4 = ...`
    with `param4 = 3` in the two input files 'file1.txt' and 'file2.txt',
    writing the two output files 'file1.txt.1234' and 'file2.txt.1234'.

    This version is case sensitive so it does not replace `PARAM1 = ...`
    with `PARAM1 = 0`, for example.

    >>> sub_params_names_case(['file1.txt', 'file1.txt'], [0, 1, 2, 3],
    ...                       ['param1', 'param2', 'param3', 'param4'],
    ...                       pid=1234)

    """
    # assert list of files
    if isinstance(files, str):
        files = [files]

    # make dict for _msub with dict[pattern] = replacement
    dd = _dict4msub(names, params)

    # replace in each file
    _msub_files(files, dd, pid)

    return


def sub_params_names_ignorecase(files, params, names, pid=None):
    """
    Substitute `name = .*` with `name = parameter` value in several files

    This means writing `for i in range(len(params))`:

       `names[i] = params[i]`

    Note, *names* are case insensitive.

    Parameters
    ----------
    files : list
        List with file in which values of *names* will be replaced with
        *params*.
    params : iterable
        | Parameter values to be given to variables on the right of the
          equal (=) sign
        | Variable in names[0] will be assigned value in params[0]
        | Variable in names[1] will be assigned value in params[1]
        | ...
    names : iterable
        Variable names on left of the equal (=) sign in *files*
    pid : int, optional
        Output files will be input files suffixed by *.pid*
        (default: overwrite input file)

    Returns
    -------
    None
        No return value but either changed input files or output files with
        names of the input files suffixed by *.pid*, in which all variables
        given in *names* are assigned the values given in *params*.

    Examples
    --------
    Replaces `param1 = ...` with `param1 = 0`, `param2 = ...` with
    `param2 = 1`, `param3 = ...` with `param3 = 2`, and `param4 = ...`
    with `param4 = 3` in the two input files 'file1.txt' and 'file2.txt',
    writing the two output files 'file1.txt.1234' and 'file2.txt.1234'.

    This version is case insensitive so it replaces also `PARAM1 = ...`,
    `Param1 = ...`, or `ParAm1 = ...` with `PARAM1 = 0`, `Param1 = 0`, or
    `ParAm1 = 0`, for example.

    >>> sub_params_names_ignorecase(['file1.txt', 'file1.txt'], [0, 1, 2, 3],
    ...                             ['param1', 'param2', 'param3', 'param4'],
    ...                             pid=1234)

    """
    # assert list of files
    if isinstance(files, str):
        files = [files]

    # make dict for _msub with dict[pattern] = replacement
    dd = _dict4msub(names, params)

    # replace in each file
    _msub_files(files, dd, pid, flags=re.I)

    return


def sub_params_names(*args, **kwargs):
    """
    Wrapper for :any:`sub_params_names_ignorecase`.

    Examples
    --------
    Replaces `param1 = ...` with `param1 = 0`, `param2 = ...` with
    `param2 = 1`, `param3 = ...` with `param3 = 2`, and `param4 = ...`
    with `param4 = 3` in the two input files 'file1.txt' and 'file2.txt',
    writing the two output files 'file1.txt.1234' and 'file2.txt.1234'.

    This version is case insensitive so it replaces also `PARAM1 = ...`,
    `Param1 = ...`, or `ParAm1 = ...` with `PARAM1 = 0`, `Param1 = 0`, or
    `ParAm1 = 0`, for example.

    >>> sub_params_names(['file1.txt', 'file1.txt'], [0, 1, 2, 3],
    ...                  ['param1', 'param2', 'param3', 'param4'],
    ...                  pid=1234)

    """
    return sub_params_names_ignorecase(*args, **kwargs)


# ------------------------------------------------------------------------------


def standard_output_reader(filename, pid=None):
    """
    Reader of simple output files

    The standard output reader reads a single value from a file without header,
    comment lines or similar. That means for example:

        0.0123456789e-02

    Parameters
    ----------
    filename : string
        Filename of file with output value
    pid : int, optional
        If present, *filename* will be suffixed by *.pid*, i.e. the file
        *filename.pid* will be read.

    Returns
    -------
    float
        Single number read from filename

    Examples
    --------
    Imagine an external program *prog.exe* that takes a number *pid* as input
    and writes out its result in *output.txt.pid*

    >>> import subprocess
    >>> exe = 'prog.exe'
    >>> pid = 1234
    >>> err = subprocess.run([exe, str(pid)])

    The results can be read with:

    >>> obj = standard_output_reader('output.txt', pid=pid)

    """
    # read output value
    if pid:
        fname = f'{filename}.{pid}'
    else:
        fname = filename
    with open(fname, 'r') as ff:
        obj = ff.readline()

    # return float
    return float(obj)


# ------------------------------------------------------------------------------


def standard_parameter_reader(filename, pid=None):
    """
    Reader of simple parameter files

    The standard parameter file is a file containing
    1 line per parameter with the parameter value.
    Lines starting with '#' will be excluded.
    That means a standard parameter file might look like:

       | # parameter value
       | 3.000000000000000e-01
       | 2.300000000000000e-01
       | 1.440000000000000e+01
       | 3.000000000000000e-01
       | ...

    Parameters
    ----------
    filename : string
        Filename with parameter values
    pid : int, optional
        If present, *filename* will be suffixed by *.pid*,
        i.e. the file *filename.pid* will be read.

    Returns
    -------
    ndarray
        Parameter values

    Examples
    --------
    A most simple parameter file *paramfile* with one line per parameter value
    can be read with:

    >>> params = standard_parameter_reader(paramfile)

    """
    if pid:
        fname = f'{filename}.{pid}'
    else:
        fname = filename
    params = []
    with open(fname, 'r') as ff:
        for line in ff:
            if line.startswith('#'):
                continue
            params.append(float(line.strip()))

    return np.array(params, dtype=float)


# ------------------------------------------------------------------------------


def standard_parameter_writer(filename, params, pid=None):
    """
    Writer of simple parameter files

    The standard parameter writer writes a file containing
    1 line per parameter with the parameter value.
    All values will be written in IEEE double precision: {:.14e}.
    That means:

       | 3.000000000000000e-01
       | 2.300000000000000e-01
       | 1.440000000000000e+01
       | 3.000000000000000e-01
       | ...

    Parameters
    ----------
    filename : string
        Output filename with parameter values
    params : iterable
        Parameter values
    pid : int, optional
        If given, output file will be *filename* suffixed by *.pid*,
        i.e. the file *filename.pid* will be written.

    Returns
    -------
    None
        No return value but output file written: *filename* or *filename.pid*

    Examples
    --------
    Writes the parameters *params* that were sampled between *bounds* (with
    the imaginary function *sample_parameters*) to file 'parameters.txt.1234'

    >>> pid = 1234
    >>> bounds = [(-5.12, 5.12)] * 5
    >>> params = sample_parameters(bounds)
    >>> standard_parameter_writer('parameters.txt', params, pid)

    """
    # Existing file will be overwritten
    if pid:
        ofile = f'{filename}.{pid}'
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
    Read a parameter file with parameter bounds and a mask

    The standard parameter reader reads from a space-separated file
    containing 1 line per parameter with the following columns:

    identifier, current parameter value, minimum parameter value,
    maximum parameter value, parameter mask (e.g., 1: include, 0: exclude).

    Lines starting with # will be excluded.

    That means a standard parameter file with bounds and mask might look like:

       | # value min max mask
       | 1 3.000000000000000e-01 0.000000000000000e+00 1.000000000000000e+00 1
       | 2 2.300000000000000e-01 -1.000000000000000e+00 1.000000000000000e+00 1
       | 3 1.440000000000000e+01 9.000000000000000e+00 2.000000000000000e+01 1
       | 4 3.000000000000000e-01 0.000000000000000e+00 1.000000000000000e+00 0
       | ...

    Parameters
    ----------
    filename : string
        Filename with parameter values
    pid : int, optional
        If present, *filename* will be suffixed by *.pid*,
        i.e. the file *filename.pid* will be read.

    Returns
    -------
    list
        List with ndarrays of

           | ids - identifier
           | params - parameter values
           | pmin - minimum parameter value
           | pmax - maximum parameter value
           | mask - parameter mask (1: include, 0: exclude)

    Examples
    --------
    Read identifiers, parameter values, lower and upper bounds, as well as
    the mask from the parameter file *paramfile*:

    >>> ids, params, pmin, pmax, pmask = (
    ...     standard_parameter_reader_bounds_mask(paramfile) )

    """
    ids    = []
    params = []
    pmin   = []
    pmax   = []
    pmask  = []
    if pid:
        pfile = f'{filename}.{pid}'
    else:
        pfile = filename
    with open(pfile, 'r') as ff:
        for line in ff:
            ll = line.strip()
            if ll.startswith('#'):
                continue
            el = ll.split()
            if len(el) != 5:
                raise IOError('Line has no 5 columns for parameter: ' + line)
            ids.append(el[0])
            params.append(float(el[1]))
            pmin.append(float(el[2]))
            pmax.append(float(el[3]))
            pmask.append(int(el[4]))

    return [ids,
            np.array(params, dtype=float),
            np.array(pmin, dtype=float),
            np.array(pmax, dtype=float),
            np.array(pmask, dtype=bool)]


# ------------------------------------------------------------------------------


def standard_parameter_writer_bounds_mask(filename, params, pmin, pmax, mask,
                                          pid=None):
    """
    Standard parameter writer with parameter bounds and a mask

    The standard parameter writer writes a space-separated file containing
    1 header line (# value min max mask) plus 1 line per parameter with the
    following columns:

    consecutive parameter number, current parameter value,
    minimum parameter value, maximum parameter value,
    parameter mask (e.g., 1: include, 0: exclude)

    All values will be written in IEEE double precision: {:.14e}.

    That means for example:

       | # value min max mask
       | 1 3.000000000000000e-01 0.000000000000000e+00 1.000000000000000e+00 1
       | 2 2.300000000000000e-01 -1.000000000000000e+00 1.000000000000000e+00 1
       | 3 1.440000000000000e+01 9.000000000000000e+00 2.000000000000000e+01 1
       | 4 3.000000000000000e-01 0.000000000000000e+00 1.000000000000000e+00 0
       | ...

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
       Parameter mask (1: include, 0: exclude)
    pid : int, optional
        If given, output file will be *filename* suffixed by *.pid*,
        i.e. the file *filename.pid* will be written.

    Returns
    -------
    None
        No return value but output file written: *filename* or *filename.pid*

    Examples
    --------
    Writes the parameters *params* that were sampled between *pmin* and
    *pmax* (with the imaginary function *sample_parameters*) to file
    'parameters.txt.1234'

    >>> pid   = 1234
    >>> p0    = np.zeros(5)
    >>> pmin  = np.full(5, -1.)
    >>> pmax  = np.full(5, 1.)
    >>> pmask = np.ones(5, dtype=int)
    >>> params = sample_parameters(p0, pmin, pmax, pmask)
    >>> standard_parameter_writer_bounds_mask(
    ...     'parameters.txt', params, pmin, pmax, pmask, pid)

    """
    # Assert correct call
    assert len(params) == len(pmin), (
        'Parameter and minima do not have the same length.')
    assert len(params) == len(pmax), (
        'Parameter and maxima do not have the same length.')
    assert len(params) == len(mask), (
        'Parameter and mask do not have the same length.')

    # Convert mask to integer if boolean
    pmask = [ int(i) for i in mask ]

    # Existing file will be overwritten
    if pid:
        fname = f'{filename}.{pid}'
    else:
        fname = filename
    with open(fname, 'w') as ff:
        # header
        hstr = '# value min max mask'
        print(hstr, file=ff)
        # data
        for i in range(len(params)):
            dstr = '{:d} {:.14e} {:.14e} {:.14e} {:d}'.format(
                i + 1, params[i], pmin[i], pmax[i], pmask[i])
            print(dstr, file=ff)

    return


# ------------------------------------------------------------------------------


def standard_time_series_reader(filename, pid=None):
    """
    Simple reader of time series without datetime

    The standard time series reader gets a time series of arbitrary length
    from a file without header, comment line or similar. The time series is
    given without any datetime information.

    That means for example:

       | 0.0123456789e-02
       | 0.1234567890e-02
       | 0.2345678900e-02
       | 0.3456789000e-02
       | 0.4567890000e-02
       | ...

    Parameters
    ----------
    filename : string
        Filename of with time series values
    pid : int, optional
        If present, *filename* will be suffixed by *.pid*,
        i.e. the file *filename.pid* will be read.

    Returns
    -------
    timeseries : ndarray
        Values of each line in filename

    Examples
    --------
    Imagine an external program *prog.exe* that takes a number *pid* as input
    and writes out a result time series in *output.txt.pid*

    >>> import subprocess
    >>> exe = 'prog.exe'
    >>> pid = 1234
    >>> err = subprocess.run([exe, str(pid)])

    The results can be read with:

    >>> ts = standard_time_series_reader('output.txt', pid=pid)

    """
    if pid:
        fname = f'{filename}.{pid}'
    else:
        fname = filename
    # read output value
    with open(fname, 'r') as ff:
        serie = ff.readlines()

    # return float array
    return np.array(serie, dtype=float)


def standard_timeseries_reader(*args, **kwargs):
    """
    Wrapper for :any:`standard_time_series_reader`

    Examples
    --------
    Imagine an external program *prog.exe* that takes a number *pid* as input
    and writes out a result time series in *output.txt.pid*

    >>> import subprocess
    >>> exe = 'prog.exe'
    >>> pid = 1234
    >>> err = subprocess.run([exe, str(pid)])

    The results can be read with:

    >>> ts = standard_timeseries_reader('output.txt', pid=pid)

    """
    return standard_time_series_reader(*args, **kwargs)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
