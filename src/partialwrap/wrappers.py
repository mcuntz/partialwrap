#!/usr/bin/env python
"""
Wrappers to partialize functions so that they can simply be called as func(x).

.. code-block:: python

   obj = partial(exe_wrapper, func,
                 parameterfile, parameterwriter,
                 outputfile, outputreader,
                 {'shell':bool, 'debug':bool, 'pid':bool,
                  'pargs':list, 'pkwargs':dict, 'keepparameterfile':bool,
                  'oargs':list, 'okwargs':dict, 'keepoutputfile':bool})
   fx = obj(x)

`func` can hence be an external executable or a Python function using
`exe_wrapper` or `function_wrapper`, respectively.

This module was written by Matthias Cuntz while at Institut National
de Recherche pour l'Agriculture, l'Alimentation et l'Environnement
(INRAE), Nancy, France.

Copyright (c) 2016-2023 Matthias Cuntz - mc (at) macu (dot) de

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
    * Rename arg to args and kwarg to kwargs, Nov 2023, Matthias Cuntz
    * Use dict.get() function for keywords, Nov 2023, Matthias Cuntz
    * Improved docstrings, Nov 2023, Matthias Cuntz
    * Removed exe_wrapper_v34 and hence support for < Python 3.7,
      Nov 2023, Matthias Cuntz

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
                kwargs, x):
    """
    Wrapper function for external programs

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
    `fx = obj(x)` by any python function such as the utilities in
    :any:`scipy`.

    Parameters
    ----------
    func : string or list of strings
        External program to be launched by :any:`subprocess`.

        External executables will be passed to
        :any:`subprocess.run`. Due to `subprocess.run()`,
        ``func`` is thus a string, e.g. ``func = './prog.exe -arg'``,
        if *shell=True*, or it is a list, e.g.
        ``func = ['./prog.exe', '-arg']``, if *shell=False*.
        If the external program is simply an executable without any arguments,
        pipes, or similar, ``func`` can simply be a string in all cases,
        e.g. ``func = './prog.exe'``.

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

    kwargs : dict
        Dictionary with keyword arguments for `exe_wrapper`. Possible
        arguments are:

        * ``shell`` (bool)
          If True, :any:`subprocess` will open a shell to run the external
          executable
        * ``debug`` (bool)
          If True, model output is displayed while the executable is running
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
          *pid=True* or *None* otherwise. The return value of ``error`` will
          then be used instead of the result of *outputreader*.

    Returns
    -------
    float
        Output value returned from *outputreader*

    Examples
    --------
    Imagine an external program *prog.exe* written in C or Fortan, for example,
    that reads in two parameters from the file *params.txt* and writes its
    results to a file *out.txt*. The parameters are written with
    `standard_parameter_writer` and output is read in with
    `standard_output_reader`.

    Then this external program can be wrapped as:

    >>> from functools import partial
    >>> from partialwrap import exe_wrapper, standard_parameter_writer
    >>> from partialwrap import standard_output_reader
    >>> exe = 'prog.exe'
    >>> parameterfile   = 'params.txt'
    >>> parameterwriter = standard_parameter_writer
    >>> outputfile      = 'out.txt'
    >>> outputreader    = standard_output_reader
    >>> exewrap = partial(exe_wrapper, exe,
    ...                   parameterfile, parameterwriter,
    ...                   outputfile, outputreader, {})

    And is run as:

    >>> x0 = [0.1, 0.2]
    >>> res = exewrap(x0)

    *res* then holds the output of *prog.exe*.

    Or one can find the minimum of the function in *prog.exe* with the global
    optimizer Differential Evolution (assuming here that *prog.exe* calculates
    the Rastrigin function with the minimum 0 at the origin).

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
    ...     return (1. + np.random.random()) * 10000.
    >>> exewrap = partial(exe_wrapper, exe,
    ...                   parameterfile, parameterwriter,
    ...                   outputfile, outputreader,
    ...                   {'error': err})
    >>> res = opt.differential_evolution(exewrap, bounds)

    """
    if sys.version_info < (3, 7):
        return exe_wrapper_v34(
            func, parameterfile, parameterwriter, outputfile, outputreader,
            kwargs, x)
    shell   = kwargs.get('shell', False)
    debug   = kwargs.get('debug', False)
    pid     = kwargs.get('pid', False)
    pargs   = kwargs.get('pargs', [])
    pkwargs = kwargs.get('pkwargs', {})
    keepparameterfile = kwargs.get('keepparameterfile', False)
    oargs   = kwargs.get('oargs', [])
    okwargs = kwargs.get('okwargs', {})
    keepoutputfile = kwargs.get('keepoutputfile', False)
    errorfunc = kwargs.get('error', '')
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

    if isinstance(func, (str, list, tuple)):
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
                     kwargs, x):
    """
    Same as :any:`exe_wrapper` incl./excl. parameters with mask.

    Makes the transformation:

    .. code-block:: python

       xx = np.copy(x0)
       xx[mask] = x

    and calls :any:`exe_wrapper` with *xx*.

    Parameters
    ----------
    func : string or list of strings
        External program to be launched by :any:`subprocess`.

        External executables will be passed to
        :any:`subprocess.run`. Due to `subprocess.run()`,
        ``func`` is thus a string, e.g. ``func = './prog.exe -arg'``,
        if *shell=True*, or it is a list, e.g.
        ``func = ['./prog.exe', '-arg']``, if *shell=False*.
        If the external program is simply an executable without any arguments,
        pipes, or similar, ``func`` can simply be a string in all cases,
        e.g. ``func = './prog.exe'``.

    x0 : array_like
        Initial values of all parameters. These will be used for parameters
        where `mask` is set to False
    mask : array_like
        Use default values *x0* in parameterfile where `mask` is False
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

    kwargs : dict
        Dictionary with keyword arguments for `exe_wrapper`. Possible
        arguments are:

        * ``shell`` (bool)
          If True, :any:`subprocess` will open a shell to run the external
          executable
        * ``debug`` (bool)
          If True, model output is displayed while the executable is running
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
          *pid=True* or *None* otherwise. The return value of ``error`` will
          then be used instead of the result of *outputreader*.

    Returns
    -------
    float
        Output value returned from *outputreader*

    Examples
    --------
    Imagine an external program *prog.exe* written in C or Fortan, for example,
    that reads in two parameters from the file *params.txt* and writes its
    results to a file *out.txt*. The parameters are written with
    `standard_parameter_writer` and the output is read in with
    `standard_output_reader`.

    Then this external program can be wrapped as:

    >>> from functools import partial
    >>> from partialwrap import exe_wrapper, standard_parameter_writer
    >>> from partialwrap import standard_output_reader
    >>> exe = 'prog.exe'
    >>> parameterfile   = 'params.txt'
    >>> parameterwriter = standard_parameter_writer
    >>> outputfile      = 'out.txt'
    >>> outputreader    = standard_output_reader
    >>> exewrap = partial(exe_wrapper, exe,
    ...                   parameterfile, parameterwriter,
    ...                   outputfile, outputreader, {})

    And is run as:

    >>> x0 = [0.1, 0.2]
    >>> res = exewrap(x0)

    *res* then holds the output of *prog.exe*.

    Or one can find the minimum of the function in *prog.exe* with the global
    optimizer Differential Evolution (assuming here that *prog.exe* calculates
    the Rastrigin function with the minimum 0 at the origin).

    >>> import scipy.optimize as opt
    >>> bounds = [(-5.12, 5.12)] * 2
    >>> res = opt.differential_evolution(exewrap, bounds)
    >>> res.x
    array([-7.32336275e-09,  6.00261842e-09])
    >>> res.fun
    0.0

    If one wants to fix the second parameter during optimization, for
    example because the parameter is insensitive or one knows the
    parameter very well, then this parameter can be excluded from the
    optimization by using `exe_mask_wrapper`. One provides then
    *initial* values for the parameters and sets a mask to *False* at
    the positions of the parameters to exclude.

    >>> from partialwrap import exe_mask_wrapper
    >>> x0   = np.array([1., 0.])
    >>> mask = np.array([True, False])
    >>> exemwrap = partial(exe_mask_wrapper, exe, x0, mask,
    ...                    parameterfile, standard_parameter_writer,
    ...                    outputfile, standard_output_reader, {})

    In case of Differential Evolution, one also has to gives bounds
    for all parameters. This would be the bounds of only the
    non-masked parameters.

    >>> bounds = [ bb for ii, bb in enumerate(bounds) if mask[ii] ]
    >>> res = opt.differential_evolution(exemwrap, bounds)
    >>> res.x
    array([-7.32336275e-09])
    >>> res.fun
    0.0

    The output has only 1 value because only 1 value was optimized
    by Differential Evolution. The final two output parameters are then:

    >>> xout = x0.copy()
    >>> xout[mask] = res.x
    >>> xout
    array([-7.32336275e-09,  0.00000000e+00])

    """
    xx = np.copy(x0)
    xx[mask] = x
    return exe_wrapper(func,
                       parameterfile, parameterwriter, outputfile,
                       outputreader, kwargs, xx)


# Python function wrappers
def function_wrapper(func, args, kwargs, x):
    """
    Wrapper for Python functions.

    To be used with :any:`functools.partial`:

    .. code-block:: python

       obj = partial(function_wrapper, func, args, kwargs)

    This allows then calling *obj* simply with the parameters *x*:
    `fx = obj(x)` by any python function such as the utilities in
    :any:`scipy`.

    Parameters
    ----------
    func : callable
        Python function to be called `func(x, *args, **kwargs)`
    args : iterable
        Arguments passed to *func*
    kwargs : dictionary
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
    return func(x, *args, **kwargs)


def function_mask_wrapper(func, x0, mask, args, kwargs, x):
    """
    Same as :any:`function_wrapper` incl./excl. parameters with mask.

    Makes the transformation:

    .. code-block:: python

       xx = np.copy(x0)
       xx[mask] = x

    and calls :any:`function_wrapper` with *xx*.

    Parameters
    ----------
    func : callable
        Python function to be called `func(x, *args, **kwargs)`
    x0 : array_like
        Initial values of all parameters. These will be used for parameters
        where `mask` is set to False
    mask : array_like
        Use default values *x0* in function call where `mask` is False
    args : iterable
        Arguments passed to *func*
    kwargs : dictionary
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

    If one wants to fix the second parameter during optimization, for example
    because the parameter is insensitive or one knows the parameter very well
    (very hypothetical in the example of rastrigin), then this parameter can be
    excluded from the optimization by using `function_mask_wrapper`.
    One then provides *initial* values of the parameters and sets a
    mask to *False* at the second index.

    >>> from partialwrap import function_mask_wrapper
    >>> x0   = np.array([1., 0., 2., 3., 4.])
    >>> mask = np.array([True, False, True, True, True])
    >>> rastram = partial(function_mask_wrapper, rastrigin, x0, mask,
    ...                   args, kwargs)

    In case of Differential Evolution, one also has to gives bounds
    for all parameters. This would be the bounds of only the
    non-masked parameters.

    >>> bounds = [ bb for ii, bb in enumerate(bounds) if mask[ii] ]
    >>> res = opt.differential_evolution(rastram, bounds)
    >>> res.x
    array([-7.32336275e-09,  1.26476721e-09, -1.61965090e-10, -7.35012009e-09])
    >>> res.fun
    0.0

    The output has only 4 values because only 4 values were optimized
    by Differential Evolution. The final five output parameters are then:

    >>> xout = x0.copy()
    >>> xout[mask] = res.x
    >>> xout
    array([-7.32336275e-09,  0.00000000e+00,  1.26476721e-09, -1.61965090e-10,
           -7.35012009e-09])

    """
    xx       = np.copy(x0)
    xx[mask] = x
    return func(xx, *args, **kwargs)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
