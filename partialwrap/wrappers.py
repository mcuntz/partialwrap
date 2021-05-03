#!/usr/bin/env python
"""
Module with wrappers to partialise functions so that they can be called simply
as func(x).

.. code-block:: python

    from functools import partial
    if isinstance(func, (str,list)):
        obj = partial(exe_wrapper, func,
                      parameterfile, parameterwriter, outputfile, outputreader,
                      {'shell':bool, 'debug':bool, 'pid':bool,
                       'pargs':list, 'pkwargs':dict, 'keepparameterfile':bool,
                       'oargs':list, 'okwargs':dict, 'keepoutputfile':bool})
    else:
        obj = partial(function_wrapper, func, arg, kwarg)
    fx = obj(x)

`func` can be a Python function or external executable.
External executables will be passed to `subprocess.check_output(func)`.
`func` can then be a string (e.g. './prog -arg') if `shell=True`
or a list (e.g. ['./prog', '-arg']) if `shell=False` (default).
Programs without arguments, pipes, etc. can simply be strings with
`shell=True` or `False`.

This module was written by Matthias Cuntz while at Institut National de
Recherche pour l'Agriculture, l'Alimentation et l'Environnement (INRAE), Nancy,
France.

Copyright (c) 2016-2020 Matthias Cuntz - mc (at) macu (dot) de

Released under the MIT License; see LICENSE file for details.

History:

* Written Nov 2016 by Matthias Cuntz (mc (at) macu (dot) de)
* Added x0 and mask to wrapper of external programs, Jan 2018, Matthias Cuntz
* Added that `pid` is passed to parameterwriter,
  and check parameterwriter (getargspec) for number or args, Feb 2018,
  Matthias Cuntz
* Removed check of number of parameters of parameterwriter (getargspec) but add
  separate wrappers for separate parmeterwriters with different number
  or arguments, Feb 2018, Matthias Cuntz
* Added `plotfile` and made docstring sphinx compatible option, Jan 2018,
  Matthias Cuntz
* Changed to Sphinx docstring and numpydoc, Nov 2019, Matthias Cuntz
* Remove that `exe_wrappers` support also Python functions. User should use
  function_wrappers, Nov 2019, Matthias Cuntz
* Make one `exe_wrapper`, passing bounds, mask, etc. via kwarg dictionary to
  parameterwriter; distinguish iterable and array_like parameter types,
  Jan 2020, Matthias Cuntz
* Replaced kwarg.pop mechanism because it removed the keywords from subsequent
  function calls, Feb 2020, Matthias Cuntz
* Change from ValueError to TypeError if function given to exe wrappers,
  Feb 2020, Matthias Cuntz
* Renamed func to function in calling names, objective to output, and renamed
  module name to wrappers, May 2020, Matthias Cuntz
* Add arguments for outputreader, its handling of pid, and the ability to have
  multiple output files, Jun 2020, Matthias Cuntz
* Correct removal of list of parameterfiles and/or outputfiles, Jun 2020,
  Matthias Cuntz
* Added ability to keep produced parameterfiles and outputfiles, Jun 2020,
  Matthias Cuntz
* Pass pid just before keyword arguments to parameterwriter and outputreader,
  Jun 2020, Matthias Cuntz
* Use `exe_wrapper` in `exe_mask_wrapper`, Jun 2020, Matthias Cuntz
* Make flake8 compliant, Dec 2020, Matthias Cuntz
* pass exe and pid as list if not shell in exe_wrapper,
  May 2021, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following wrappers are provided:

.. autosummary::
   exe_wrapper
   exe_mask_wrapper
   function_wrapper
   function_mask_wrapper
"""
from __future__ import division, absolute_import, print_function
import subprocess
import os
import numpy as np


__all__ = ['exe_wrapper', 'exe_mask_wrapper',
           'function_wrapper', 'function_mask_wrapper']


def _tolist(arg):
    """
    Assure that `arg` is a list, e.g. if string or None are given.

    Parameters
    ----------
    arg :
        argument to make list

    Returns
    -------
    list
        list(arg)

    Examples
    --------
    >>> _tolist('string')
    ['string']
    >>> _tolist([1,2,3])
    [1, 2, 3]
    >>> _tolist(None)
    [None]
    """
    if isinstance(arg, str):
        return [arg]
    try:
        return list(arg)
    except TypeError:
        return [arg]


def exe_wrapper(func,
                parameterfile, parameterwriter, outputfile, outputreader,
                kwarg, x):
    """
    Wrapper function for external programs using a `parameterwriter` and
    `outputreader` with the interfaces:

        `parameterwriter(parameterfile, x, *pargs, **pkwargs)`

        `outputreader(outputfile, *oargs, **okwargs)`

    or if `pid==True`:

        `parameterwriter(parameterfile, x, *pargs, pid=pid, **pkwargs)`

        `outputreader(outputfile, *oargs, pid=pid, **okwargs)`

    Examples of `parameterwriter` with `pid==True` are
    :any:`standard_parameter_writer` or :any:`sub_params_ja`.
    Examples of `outputreader` with or without `pid` is
    :any:`standard_output_reader`.

    To be used with :any:`functools.partial`:

        `obj = partial(exe_wrapper, func,
        parameterfile, parameterwriter, outputfile, outputreader,
        {'shell':bool, 'debug':bool, 'pid':bool,
        'pargs':list, 'pkwargs':dict, 'keepparameterfile':bool,
        'oargs':list, 'okwargs':dict, 'keepoutputfile':bool})`

    This allows then calling obj with only the non-masked parameters:

        `fx = obj(x)`

    which translates to:

        `parameterwriter(parameterfile, x, *pargs, **pkwargs)`

        `err = subprocess.check_output(
             func, stderr=subprocess.STDOUT, shell=shell)`

        `obj = outputreader(outputfile, *oargs, **okwargs)`

    or if `pid==True` to:

        `parameterwriter(parameterfile, x, *pargs, pid=pid, **pkwargs)`

        `err = subprocess.check_output(
            func, stderr=subprocess.STDOUT, shell=shell)`

        `obj = outputreader(outputfile, *oargs, pid=pid, **okwargs)`

    Note if `pid==True`, it is assumed that `parameterwriter` produces
    files in the form `parameterfile+'.'+processid`, which will be removed
    after the function call.
    Also only files in the form `outputfile+'.'+processid` will be removed
    automatically.

    Parameters
    ----------
    func : string or list of strings
        External program to launch by :any:`subprocess`
    parameterfile : string or iterable
        Filename(s) of parameter file(s)
    parameterwriter : callable
        Python function writing the `parameterfile`, called as:

            `parameterwriter(parameterfile, x, *pargs, **pkwargs)`

        or if `pid==True` as:

            `parameterwriter(parameterfile, x, *pargs, pid=pid, **pkwargs)`

    outputfile : string or iterable
        Filename(s) of file(s) with output values written by the external
        executable `func`
    outputreader : callable
        Python function for reading and processing output value(s) from
        `outputfile`, called as:

            `outputreader(ouputfile, x, *oargs, **okwargs)`

        or if `pid==True` as:

            `outputreader(ouputfile, x, *oargs, pid=pid, **okwargs)`

    kwarg : dict
        Dictionary with keyword arguments for `exe_wrapper`. Possible
        arguments are:

            **shell** (bool)

                If True, :any:`subprocess` opens shell for external executable

            **debug** (bool)

                If True, model output is displayed while executable is running

            **pid** (bool)

                If True, append '.RandomNumber' to `parameterfile` and
                `outputfile` for parallel calls of `func`

            **pargs** (iterable)

                List of arguments of `parameterwriter`.

            **pkwargs** (dict)

                Dictionary with keyword arguments of `parameterwriter`.

            **keepparameterfile** (bool)

                If True, `parameterfile` produced with `parameterwriter` will
                not be deleted after function call.

            **oargs** (iterable)

                List of arguments of `outputreader`.

            **okwargs** (dict)

                Dictionary with keyword arguments of `outputreader`.

            **keepoutputfile** (bool)

                If True, `outputfile` produced by `func` will not be deleted
                after function call.

    Returns
    -------
    float
        Output value calculated by the external executable `func` or via
        the `outputreader`


    History
    -------
    Written,  Matthias Cuntz, Mar 2018
    Modified, Matthias Cuntz, Dec 2019 - rm Python function
                                         -> use `function_wrapper`
                                       - Sphinx docstring
              Matthias Cuntz, Jan 2020 - renamed to `exe_wrapper`
                                       - shell, debug in kwarg dictionary
                                       - remove parameterfile if it exists
                                       - pid switch in kwarg dictionary
                                       - pargs, pkwargs for passing bounds,
                                         mask, etc. to parameterwriter
                                         replacing other `exe_wrappers`
                                       - distinguish iterable and `array_like`
                                         parameter types
              Matthias Cuntz, Feb 2020 - replaced kwarg.pop because it removed
                                         keywords from subsequent calls
                                       - ValueError -> TypeError
              Matthias Cuntz, May 2020 - renamed objective to output
              Matthias Cuntz, Jun 2020 - added arguments to outputreader
                                       - changed use of pid with outputfile
                                       - outputfile can be list of files
                                       - rm also list of parameterfiles and
                                         outputfiles
                                       - keekparameterfile, keepoutputfile
              Matthias Cuntz, Jun 2020 - pass pid just before keyword arguments
                                         to parameterwriter and outputreader
              Matthias Cuntz, Jun 2020 - called from `exe_mask_wrapper`
              Matthias Cuntz, May 2021 - pass func and pid as list if not shell
    """
    shell   = kwarg['shell']   if 'shell'   in kwarg else False
    debug   = kwarg['debug']   if 'debug'   in kwarg else False
    pid     = kwarg['pid']     if 'pid'     in kwarg else False
    pargs   = kwarg['pargs']   if 'pargs'   in kwarg else []
    pkwargs = kwarg['pkwargs'] if 'pkwargs' in kwarg else {}
    keepparameterfile = (kwarg['keepparameterfile']
                         if 'keepparameterfile' in kwarg else False)
    oargs   = kwarg['oargs']   if 'oargs'   in kwarg else []
    okwargs = kwarg['okwargs'] if 'okwargs' in kwarg else {}
    keepoutputfile = (kwarg['keepoutputfile']
                      if 'keepoutputfile' in kwarg else False)
    # For multiprocess but not MPI: pid = mp.current_process()._identity[0]
    # seed uses global variables so all processes will produce same random
    # numbers
    # -> use np.random.RandomState() for each processes for individual seeds in
    # each process
    if pid:
        randst = np.random.RandomState()
        ipid   = str(randst.randint(2147483647))
    if isinstance(func, (str, list)):
        if pid:
            parameterwriter(parameterfile, x, *pargs, pid=ipid, **pkwargs)
            if isinstance(func, str):
                if shell:
                    func1 = func + ' ' + ipid
                else:
                    func1 = [func, ipid]
            else:
                func1 = func + [ipid]
            if debug:
                err = subprocess.check_call(func1, stderr=subprocess.STDOUT,
                                            shell=shell)
            else:
                err = subprocess.check_output(func1, stderr=subprocess.STDOUT,
                                              shell=shell)
            obj = outputreader(outputfile, *oargs, pid=ipid, **okwargs)
            if not keepparameterfile:
                for ff in _tolist(parameterfile):
                    if os.path.exists(ff + '.' + ipid):
                        os.remove(ff + '.' + ipid)
            if not keepoutputfile:
                for ff in _tolist(outputfile):
                    if os.path.exists(ff + '.' + ipid):
                        os.remove(ff + '.' + ipid)
        else:
            parameterwriter(parameterfile, x, *pargs, **pkwargs)
            if debug:
                err = subprocess.check_call(func, stderr=subprocess.STDOUT,
                                            shell=shell)
            else:
                err = subprocess.check_output(func, stderr=subprocess.STDOUT,
                                              shell=shell)
            obj = outputreader(outputfile, *oargs, **okwargs)
            if not keepparameterfile:
                for ff in _tolist(parameterfile):
                    if os.path.exists(ff):
                        os.remove(ff)
            if not keepoutputfile:
                for ff in _tolist(outputfile):
                    if os.path.exists(ff):
                        os.remove(ff)
        return obj
    else:
        estr  = 'func must be string or list of strings for subprocess.'
        estr += ' Use function_wrapper for Python functions.'
        raise TypeError(estr)


def exe_mask_wrapper(func, x0, mask,
                     parameterfile, parameterwriter, outputfile, outputreader,
                     kwarg, x):
    """
    Same as `exe_wrapper` incl./excl. parameters with mask.

    Makes the transformation:

        `xx = np.copy(x0)`

        `xx[mask] = x`

    and calls `exe_wrapper` with `xx`.

    See `exe_wrapper` for details.
    """
    xx = np.copy(x0)
    xx[mask] = x
    return exe_wrapper(func,
                       parameterfile, parameterwriter, outputfile,
                       outputreader, kwarg, xx)


# Python function wrappers
def function_wrapper(func, arg, kwarg, x):
    """
    Wrapper function for Python function.
    To be used with partial:

        `obj = partial(function_wrapper, func, arg, kwarg)`

    This allows then calling obj with only the non-masked parameters:

        `fx = obj(x)`

    which translates to:

        `fx = func(x, *arg, **kwarg)`

    Parameters
    ----------
    func : callable
        Python function to be called `func(x, *arg, **kwarg)`
    arg : iterable
        Arguments passed to `func`
    kwarg : dictionary
        Keyword arguments passed to `func`

    Returns
    -------
    float
        Output value calculated by `func`


    History
    -------
    Written,  Matthias Cuntz, Nov 2016
    Modified, Matthias Cuntz, Nov 2019 - Sphinx docstring
    Matthias Cuntz, May 2020 - renamed func to function in name
    """
    return func(x, *arg, **kwarg)


def function_mask_wrapper(func, x0, mask, arg, kwarg, x):
    """
    Wrapper function for Python function using a mask.
    To be used with partial:

        `obj = partial(function_mask_wrapper, func, x0, mask, arg, kwarg)`

    This allows then calling obj with only the non-masked parameters:

        `fx = obj(x)`

    which translates to:

        `xx = np.copy(x0)`

        `xx[mask] = x`

        `fx = func(xx, *arg, **kwarg)`

    Parameters
    ----------
    func : callable
        Python function to be called `func(x, *arg, **kwarg)`
    x0 : array_like
        Fixed values of masked parameters
    mask : array_like
        Mask to use `x0` values ('mask[i]=1') or use new parameters
        `x` ('mask[i]=0') in call of function
    arg : iterable
        Arguments passed to `func`
    kwarg : dictionary
        Keyword arguments passed to `func`

    Returns
    -------
    float
        Output value calculated by `func`


    History
    -------
    Written,  Matthias Cuntz, Nov 2016
    Modified, Matthias Cuntz, Nov 2019 - Sphinx docstring
              Matthias Cuntz, Jan 2020 - distinguish iterable and `array_like`
                                         parameter types
              Matthias Cuntz, May 2020 - renamed func to function in name
    """
    xx       = np.copy(x0)
    xx[mask] = x
    return func(xx, *arg, **kwarg)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
