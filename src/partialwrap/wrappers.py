#!/usr/bin/env python
"""
Wrappers to partialise functions so that they can be simply called as func(x).

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

`func` can be a Python function or external executable. External executables
will be passed to `subprocess.run(func)`. `func` can then be a string (e.g.
'./prog -arg') if `shell=True` or a list (e.g. ['./prog', '-arg']) if
`shell=False` (default). Programs without arguments, pipes, etc. can simply be
strings with `shell=True` or `False`.

This module was written by Matthias Cuntz while at Institut National de
Recherche pour l'Agriculture, l'Alimentation et l'Environnement (INRAE), Nancy,
France.

Copyright (c) 2016-2022 Matthias Cuntz - mc (at) macu (dot) de

Released under the MIT License; see LICENSE file for details.

.. moduleauthor:: Matthias Cuntz

The following wrappers are provided:

.. autosummary::
   exe_wrapper
   exe_mask_wrapper
   function_wrapper
   function_mask_wrapper

History
    * Written Nov 2016 by Matthias Cuntz (mc (at) macu (dot) de)
    * Added x0 and mask to wrapper of external programs,
      Jan 2018, Matthias Cuntz
    * Added that `pid` is passed to parameterwriter,
      and check parameterwriter (getargspec) for number or args, Feb 2018,
      Matthias Cuntz
    * Removed check of number of parameters of parameterwriter (getargspec) but
      add separate wrappers for separate parmeterwriters with different number
      or arguments, Feb 2018, Matthias Cuntz
    * Added `plotfile` and made docstring sphinx compatible option, Jan 2018,
      Matthias Cuntz
    * Changed to Sphinx docstring and numpydoc, Nov 2019, Matthias Cuntz
    * Remove that `exe_wrappers` support also Python functions. User should use
      function_wrappers, Nov 2019, Matthias Cuntz
    * Make one `exe_wrapper`, passing bounds, mask, etc. via kwarg dictionary
      to parameterwriter; distinguish iterable and array_like parameter types,
      Jan 2020, Matthias Cuntz
    * Replaced kwarg.pop mechanism because it removed the keywords from
      subsequent function calls, Feb 2020, Matthias Cuntz
    * Change from ValueError to TypeError if function given to exe wrappers,
      Feb 2020, Matthias Cuntz
    * Renamed func to function in calling names, objective to output, and
      renamed module name to wrappers, May 2020, Matthias Cuntz
    * Add arguments for outputreader, its handling of pid, and the ability to
      have multiple output files, Jun 2020, Matthias Cuntz
    * Correct removal of list of parameterfiles and/or outputfiles, Jun 2020,
      Matthias Cuntz
    * Added ability to keep produced parameterfiles and outputfiles, Jun 2020,
      Matthias Cuntz
    * Pass pid just before keyword arguments to parameterwriter and
      outputreader, Jun 2020, Matthias Cuntz
    * Use `exe_wrapper` in `exe_mask_wrapper`, Jun 2020, Matthias Cuntz
    * Make flake8 compliant, Dec 2020, Matthias Cuntz
    * Pass exe and pid as list if not shell in exe_wrapper,
      May 2021, Matthias Cuntz
    * Change to subprocess.run for Python >= v3.5, Aug 2022, Matthias Cuntz
    * Added examples and reformatted docstrings, Aug 2022, Matthias Cuntz
    * Added error function as fallback option in case subprocess exits with
      error code > 0, Aug 2022, Matthias Cuntz

"""
import sys
import subprocess as sp
import os
import numpy as np


__all__ = ['exe_wrapper', 'exe_mask_wrapper',
           'function_wrapper', 'function_mask_wrapper']


def _tolist(arg):
    """
    Assure that *arg* is a list, e.g. if string or None are given.

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
    Wrapper function for external programs using a *parameterwriter* and
    *outputreader* with the interfaces:
    ``parameterwriter(parameterfile, x, *pargs, **pkwargs)``
    and
    ``outputreader(outputfile, *oargs, **okwargs)``
    or if *pid==True*:
    ``parameterwriter(parameterfile, x, *pargs, pid=pid, **pkwargs)``
    and
    ``outputreader(outputfile, *oargs, pid=pid, **okwargs)``

    Examples of *parameterwriter* with *pid==True* are
    :any:`standard_parameter_writer` or :any:`sub_params_ja`.
    An example of *outputreader* with or without *pid* is
    :any:`standard_output_reader`.

    To be used with :any:`functools.partial`:

    .. code-block:: python

       obj = partial(exe_wrapper, func,
                     parameterfile, parameterwriter,
                     outputfile, outputreader,
                     {'shell': bool, 'debug': bool, 'pid': bool,
                      'pargs': list, 'pkwargs': dict,
                      'keepparameterfile': bool,
                      'oargs': list, 'okwargs': dict,
                      'keepoutputfile': bool})

    This allows then calling *obj* simply with the parameters *x*:
    ``fx = obj(x)``

    which translates to:

    .. code-block:: python

       parameterwriter(parameterfile, x, *pargs, **pkwargs)
       err = subprocess.run(func)
       obj = outputreader(outputfile, *oargs, **okwargs)

    or if *pid==True* to:

    .. code-block:: python

       parameterwriter(parameterfile, x, *pargs, pid=pid, **pkwargs)
       err = subprocess.run(func)
       obj = outputreader(outputfile, *oargs, pid=pid, **okwargs)

    Note if *pid==True*, it is assumed that *parameterwriter* produces
    files in the form `parameterfile + '.' + processid`,
    which will be removed after the function call.
    Also files in the form `outputfile + '.' + processid`
    will be removed automatically.

    Parameters
    ----------
    func : string or list of strings
        External program to launch by :any:`subprocess`
    parameterfile : string or iterable
        Filename(s) of parameter file(s)
    parameterwriter : callable
        Python function writing the *parameterfile*, called as:
        ``parameterwriter(parameterfile, x, *pargs, **pkwargs)``

        or if *pid==True* as:
        ``parameterwriter(parameterfile, x, *pargs, pid=pid, **pkwargs)``
    outputfile : string or iterable
        Filename(s) of file(s) with output values written by the external
        executable *func*
    outputreader : callable
        Python function for reading and processing output value(s) from
        *outputfile*, called as:
        ``outputreader(ouputfile, x, *oargs, **okwargs)``

        or if *pid==True* as:
        ``outputreader(ouputfile, x, *oargs, pid=pid, **okwargs)``
    kwarg : dict
        Dictionary with keyword arguments for `exe_wrapper`. Possible
        arguments are:

        * ``shell`` (bool)
          If True, :any:`subprocess` opens shell for external executable
        * ``debug`` (bool)
          If True, model output is displayed while executable is running
        * ``pid`` (bool)
          If True, append '.RandomNumber' to *parameterfile* and
          *outputfile* for parallel calls of *func*
        * ``pargs`` (iterable)
          List of arguments of `parameterwriter`.
        * ``pkwargs`` (dict)
          Dictionary with keyword arguments of `parameterwriter`.
        * ``keepparameterfile`` (bool)
          If True, `parameterfile` produced with `parameterwriter` will
          not be deleted after function call.
        * ``oargs`` (iterable)
          List of arguments of `outputreader`.
        * ``okwargs`` (dict)
          Dictionary with keyword arguments of `outputreader`.
        * ``keepoutputfile`` (bool)
          If True, `outputfile` produced by `func` will not be deleted
          after function call.
        * ``error`` (function)
          Function to call in case subprocess exists with error code > 0. The
          function should take one argument, which will be the random number if
          `pid=True` or *None* otherwise. The return value will be used instead
          of the result of *outputreader*.
          This is not available for Python < v3.7

    Returns
    -------
    float
        Output value calculated by the external executable *func* or via
        the *outputreader*

    Examples
    --------
    Imagine an external program *prog.exe* written in C or Fortan, for example,
    that reads in two parameters from the file *params.txt* and writes its
    results to a file *out.txt*. The parameters can be written with
    `standard_parameter_writer` and read with `standard_output_reader`.

    Then this external program can be wrapped as:

    >>> from functools import partial
    >>> from partialwrap import exe_wrapper, standard_parameter_writer
    >>> from partialwrap import standard_output_reader
    >>> exe = 'prog.exe'
    >>> parameterfile = 'params.txt'
    >>> outputfile = 'out.txt'
    >>> exewrap = partial(exe_wrapper, exe,
    ...                   parameterfile, standard_parameter_writer,
    ...                   outputfile, standard_output_reader, {})

    And is run as:

    >>> x0 = [0.1, 0.2]
    >>> res = exewrap(x0)

    *res* then holds the output of *prog.exe*.

    Or one can find the minimum of the function in *prog.exe* with the global
    optimizer Differential Evolution (assuming here that *prog.exe* calculates
    the Rastrigin function).

    >>> import scipy.optimize as opt
    >>> bounds = [(-5.12, 5.12)] * 2
    >>> res = opt.differential_evolution(exewrap, bounds)
    >>> res.x
    array([-7.32336275e-09,  6.00261842e-09])
    >>> res.fun
    0.0

    If *prog.exe* might fail with some parameter combinations from the
    Differential Evolution algorithm, then one can give a fallback function for
    this cases that, for example, returns arbitrary large values.

    >>> def err(x):
    ...     return np.random.random() * 10000.
    >>> exewrap = partial(exe_wrapper, exe,
    ...                   parameterfile, standard_parameter_writer,
    ...                   outputfile, standard_output_reader,
    ...                   {'error': err})
    >>> res = opt.differential_evolution(exewrap, bounds)

    """
    if sys.version_info < (3, 7):
        return exe_wrapper_v34(
            func, parameterfile, parameterwriter, outputfile, outputreader,
            kwarg, x)
    shell   = kwarg['shell'] if 'shell' in kwarg else False
    debug   = kwarg['debug'] if 'debug' in kwarg else False
    pid     = kwarg['pid'] if 'pid' in kwarg else False
    pargs   = kwarg['pargs'] if 'pargs' in kwarg else []
    pkwargs = kwarg['pkwargs'] if 'pkwargs' in kwarg else {}
    keepparameterfile = (kwarg['keepparameterfile']
                         if 'keepparameterfile' in kwarg else False)
    oargs   = kwarg['oargs'] if 'oargs' in kwarg else []
    okwargs = kwarg['okwargs'] if 'okwargs' in kwarg else {}
    keepoutputfile = (kwarg['keepoutputfile']
                      if 'keepoutputfile' in kwarg else False)
    errorfunc = kwarg['error'] if 'error' in kwarg else ''
    # For multiprocess but not MPI: pid = mp.current_process()._identity[0]
    # seed uses global variables so all processes will produce the same
    # random numbers
    # -> use np.random.RandomState() for each process to have individual
    #    seeds in each process
    if pid:
        randst = np.random.RandomState()
        ipid = str(randst.randint(2147483647))
        pipid = '.' + ipid
    else:
        ipid = None
        pipid = ''
    if isinstance(func, (str, list)):
        if pid:
            parameterwriter(parameterfile, x, *pargs, pid=ipid, **pkwargs)
            if isinstance(func, str):
                if shell:
                    func1 = func + ' ' + ipid
                else:
                    func1 = [func, ipid]
            else:
                if shell:
                    func1 = ' '.join(func) + ' ' + ipid
                else:
                    func1 = func + [ipid]
        else:
            parameterwriter(parameterfile, x, *pargs, **pkwargs)
            if shell and (not isinstance(func, str)):
                func1 = ' '.join(func)
            else:
                func1 = func
        if debug:
            spkwarg = {'stderr': sp.PIPE}
        else:
            spkwarg = {'capture_output': True}
        try:
            err = sp.run(func1, check=True, shell=shell, **spkwarg)
            if pid:
                obj = outputreader(outputfile, *oargs, pid=ipid, **okwargs)
            else:
                obj = outputreader(outputfile, *oargs, **okwargs)
        except sp.CalledProcessError as e:
            if errorfunc:
                obj = errorfunc(ipid)
            else:
                print(e.cmd, 'returned with exit code', e.returncode)
                if not debug:
                    print('stdout was:', e.stdout)
                    print('stderr was:', e.stderr)
                raise ValueError('Subprocess error')
        if not keepparameterfile:
            for ff in _tolist(parameterfile):
                if os.path.exists(ff + pipid):
                    os.remove(ff + pipid)
        if not keepoutputfile:
            for ff in _tolist(outputfile):
                if os.path.exists(ff + pipid):
                    os.remove(ff + pipid)
        return obj
    else:
        raise TypeError('func must be string or list of strings for'
                        ' subprocess. Use function_wrapper for Python'
                        ' functions.')


def exe_mask_wrapper(func, x0, mask,
                     parameterfile, parameterwriter, outputfile, outputreader,
                     kwarg, x):
    """
    Same as `exe_wrapper` incl./excl. parameters with mask.

    Makes the transformation:

    .. code-block:: python

       xx = np.copy(x0)
       xx[mask] = x

    and calls `exe_wrapper` with *xx*.

    See `exe_wrapper` for details.

    Examples
    --------
    Imagine an external program *prog.exe* written in C or Fortan, for example,
    that reads in two parameters from the file *params.txt* and writes its
    results to a file *out.txt*. The parameters can be written with
    `standard_parameter_writer` and read with `standard_output_reader`.

    Then this external program can be wrapped as:

    >>> from functools import partial
    >>> from partialwrap import exe_wrapper, standard_parameter_writer
    >>> from partialwrap import standard_output_reader
    >>> exe = 'prog.exe'
    >>> parameterfile = 'params.txt'
    >>> outputfile = 'out.txt'
    >>> exewrap = partial(exe_wrapper, exe,
    ...                   parameterfile, standard_parameter_writer,
    ...                   outputfile, standard_output_reader, {})

    And is run as:

    >>> x0 = [0.1, 0.2]
    >>> res = exewrap(x0)

    *res* then holds the output of *prog.exe*.

    Or one can find the minimum of the function in *prog.exe* with the global
    optimizer Differential Evolution (assuming here that *prog.exe* calculates
    the Rastrigin function).

    >>> import scipy.optimize as opt
    >>> bounds = [(-5.12, 5.12)] * 2
    >>> res = opt.differential_evolution(exewrap, bounds)
    >>> res.x
    array([-7.32336275e-09,  6.00261842e-09])
    >>> res.fun
    0.0

    If one wants to fix the second parameter during optimisation, for example
    because the parameter is insensitive or one knows the parameter very well,
    then this parameter can be excluded from the optimisation by using
    `exe_mask_wrapper`, setting the mask to *False*, and providing *initial*
    values of the parameters.

    >>> from partialwrap import exe_mask_wrapper
    >>> mask = np.array([True, False])
    >>> x0 = np.array([1., 0.])
    >>> exemwrap = partial(exe_mask_wrapper, exe, x0, mask,
    ...                    parameterfile, standard_parameter_writer,
    ...                    outputfile, standard_output_reader, {})

    >>> bounds = [ bb for ii, bb in enumerate(bounds) if mask[ii] ]
    >>> res = opt.differential_evolution(exemwrap, bounds)
    >>> res.x
    array([-7.32336275e-09])
    >>> res.fun
    0.0

    The final output parameters are then:

    >>> xout = x0.copy()
    >>> xout[mask] = res.x
    >>> xout
    array([-7.32336275e-09,  0.00000000e+00])

    """
    xx = np.copy(x0)
    xx[mask] = x
    return exe_wrapper(func,
                       parameterfile, parameterwriter, outputfile,
                       outputreader, kwarg, xx)


def exe_wrapper_v34(func,
                    parameterfile, parameterwriter, outputfile, outputreader,
                    kwarg, x):
    """
    exe_wrapper for Python < v3.7 using subprocess.check_call and
    subprocess.check_output

    See `exe_wrapper` for details.

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
                err = sp.check_call(func1, stderr=sp.STDOUT,
                                    shell=shell)
            else:
                err = sp.check_output(func1, stderr=sp.STDOUT,
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
                err = sp.check_call(func, stderr=sp.STDOUT,
                                    shell=shell)
            else:
                err = sp.check_output(func, stderr=sp.STDOUT,
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


# Python function wrappers
def function_wrapper(func, arg, kwarg, x):
    """
    Wrapper for Python functions.

    To be used with partial:
    ``obj = partial(function_wrapper, func, arg, kwarg)``

    This allows then calling obj with only the non-masked parameters:
    ``fx = obj(x)``

    which translates to:
    ``fx = func(x, *arg, **kwarg)``

    Parameters
    ----------
    func : callable
        Python function to be called `func(x, *arg, **kwarg)`
    arg : iterable
        Arguments passed to *func*
    kwarg : dictionary
        Keyword arguments passed to *func*

    Returns
    -------
    float
        Output value calculated by *func*

    Examples
    --------
    Having the Python function *rastrigin* with the parameter *a* and
    the keyword argument *b*.

    >>> import numpy as np
    >>> def rastrigin(x, a, b=2.*np.pi):
    ...     return a * x.size + np.sum(x**2 - a * np.cos(b * x))

    We want to evaluate the 5-dimensional rastrigin function at
    `[-1, 0, 1, 2, 3]` using the parameters `a=20.`, `b=np.pi`:

    >>> from functools import partial
    >>> from partialwrap import function_wrapper
    >>> args   = [20.]
    >>> kwargs = {'b': 1.*np.pi}
    >>> rastra = partial(function_wrapper, rastrigin, args, kwargs)
    >>> x0 = np.array([-1, 0, 1, 2, 3])
    >>> res = rastra(x0)
    >>> res
    135.0

    Or find its minimum with the global optimizer Differential Evolution,
    which allows passing arguments such as `args=(20.)` but does not allow
    passing keyword arguments such `b=np.pi`.

    >>> import scipy.optimize as opt
    >>> bounds = [(-5.12, 5.12)] * 5
    >>> res = opt.differential_evolution(rastra, bounds)
    >>> res.x
    array([-7.32336275e-09,  6.00261842e-09,  1.26476721e-09, -1.61965090e-10,
           -7.35012009e-09])
    >>> res.fun
    0.0

    """
    return func(x, *arg, **kwarg)


def function_mask_wrapper(func, x0, mask, arg, kwarg, x):
    """
    Wrapper function for Python function using a mask.

    To be used with partial:
    ``obj = partial(function_mask_wrapper, func, x0, mask, arg, kwarg)``

    This allows then calling obj with only the non-masked parameters:
    ``fx = obj(x)``

    which translates to:

    .. code-block:: python

       xx = np.copy(x0)
       xx[mask] = x
       fx = func(xx, *arg, **kwarg)

    Parameters
    ----------
    func : callable
        Python function to be called `func(x, *arg, **kwarg)`
    x0 : array_like
        Fixed values of masked parameters
    mask : array_like
        Mask to use *x0* values (`mask[i]=1`) or use new parameters
        *x* (`ask[i]=0`) in call of the function
    arg : iterable
        Arguments passed to *func*
    kwarg : dictionary
        Keyword arguments passed to *func*

    Returns
    -------
    float
        Output value calculated by *func*

    Examples
    --------
    Having the Python function *rastrigin* with the parameter *a* and
    the keyword argument *b*.

    >>> import numpy as np
    >>> def rastrigin(x, a, b=2.*np.pi):
    ...     return a * x.size + np.sum(x**2 - a * np.cos(b * x))

    We want to evaluate the 5-dimensional rastrigin function at
    `[-1, 0, 1, 2, 3]` using the parameters `a=20.`, `b=np.pi`:

    >>> from functools import partial
    >>> from partialwrap import function_wrapper
    >>> args   = [20.]
    >>> kwargs = {'b': 1.*np.pi}
    >>> rastra = partial(function_wrapper, rastrigin, args, kwargs)
    >>> x0 = np.array([-1, 0, 1, 2, 3])
    >>> res = rastra(x0)
    >>> res
    135.0

    Or we want find its minimum with the global optimizer Differential
    Evolution, which allows passing arguments such as `args=(20.)` but does not
    allow passing keyword arguments such `b=np.pi`.

    >>> import scipy.optimize as opt
    >>> bounds = [(-5.12, 5.12)] * 5
    >>> res = opt.differential_evolution(rastra, bounds)
    >>> res.x
    array([-7.32336275e-09,  6.00261842e-09,  1.26476721e-09, -1.61965090e-10,
           -7.35012009e-09])
    >>> res.fun
    0.0

    If one wants to fix the second parameter during optimisation, for example
    because the parameter is insensitive or one knows the parameter very well
    (very hypothetical in the example of rastrigin), then this parameter can be
    excluded from the optimisation by using `function_mask_wrapper`, setting
    the mask to *False*, and providing *initial* values of the parameters.

    >>> from partialwrap import function_mask_wrapper
    >>> mask = np.array([True, False, True, True, True])
    >>> x0 = np.array([1., 0., 2., 3., 4.])
    >>> rastram = partial(function_mask_wrapper, rastrigin, x0, mask,
    ...                   args, kwargs)

    >>> bounds = [ bb for ii, bb in enumerate(bounds) if mask[ii] ]
    >>> res = opt.differential_evolution(rastram, bounds)
    >>> res.x
    array([-7.32336275e-09,  1.26476721e-09, -1.61965090e-10, -7.35012009e-09])
    >>> res.fun
    0.0

    The final output parameters are then:

    >>> xout = x0.copy()
    >>> xout[mask] = res.x
    >>> xout
    array([-7.32336275e-09,  0.00000000e+00,  1.26476721e-09, -1.61965090e-10,
           -7.35012009e-09])

    """
    xx       = np.copy(x0)
    xx[mask] = x
    return func(xx, *arg, **kwarg)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
