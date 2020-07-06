**********
User Guide
**********

``partialwrap`` provides wrappers for external executables and Python functions so that they can
easily be partialised with Python's :func:`functools.partial`.

Partial's ingenious mechanism allows to use even very complex functions with many arguments and
keyword arguments with routines that need functions in the simple form `func(x)`. This includes
Python's :func:`map` function, the minimizers in :mod:`scipy.optimize`, and plenty of third-party
modules such as :mod:`emcee` or :mod:`pyeee`. It also allows easy parallelization of code with, for
example, the parallel :func:`~multiprocessing.pool.Pool.map` function of the :mod:`multiprocessing`
module or via the Message Passing Interface (MPI) :mod:`mpi4py`.

``partialwrap`` provides two easy wrappers so that :func:`functools.partial` can be used with
external programs as well as Python functions. It provides further functions to deal with
input/output external programs. The user guide gives examples how to use pretty much every external
program within Python modules.


Simple Python functions
=======================

The easiest wrapper in ``partialwrap`` is for Python functions.

Consider the *Rastrigin function* (https://en.wikipedia.org/wiki/Rastrigin_function), which is a
popular function for performance testing of optimization algorithms: :math:`y = a \cdot n +
\sum_i^n (x_i^2 - a \cos(b \cdot x_i))`:

.. code-block:: python

   import numpy as np
   def rastrigin(x, a, b=2.*np.pi):
       return a*x.size + np.sum(x**2 - a*np.cos(b*x))

The Rastrigin function has a global minimum of :math:`0` at all :math:`x_i = 0`. :math:`a` influences
mainly the depths of the (local and global) minima, whereas :math:`b` influences mainly the size of
the minima.

A common form uses :math:`a = 10` and :math:`b = 2 \pi`. The parameters :math:`x_i` should then be
in the interval :math:`[-5.12,5.12]`.

The 5-dimensional Rastrigin function is hence called in Python as:

.. code-block:: python

   ndim = 5
   a    = 10.
   xmin = -5.12
   xmax =  5.12
   x = xmin + (xmax - xmin) * np.random.random_sample(ndim)
   y = rastrigin(x, a)

One can search the global minimum of the Rastrigin function using, for example, the differential
evolution algorithm of the :mod:`scipy.optimize` package:

.. code-block:: python

   import scipy.optimize as opt
   res = opt.differential_evolution(rastrigin, [(xmin,xmax),]*ndim, args=np.array([a]))

:mod:`scipy.optimize` does not allow to pass keyword arguments to the function. It is hence *a
priori* not possible to search the minimum of the Rastrigin function with another :math:`b`
parameter. (This is an illustrative example while it is, of course, possible to just pass `b` as a
second argument in this case.)

In this case, one can use Python's :func:`functools.partial` function:

.. code-block:: python

   from functools import partial

   # helper function
   def call_func_arg_kwarg(func, a, b, x):
       return func(x, a, b=b)

   # Partialise function with fixed parameters
   a = 5.
   b = 4.*np.pi
   partial_rastrigin = partial(call_func_arg_kwarg, rastrigin, a, b)

   # search minimum
   res = opt.differential_evolution(partial_rastrigin, [(xmin,xmax),]*ndim)

Figuratively speaking, :func:`~functools.partial` passes :math:`a` and :math:`b` already during definition
to the function `call_func_arg_kwarg`. :func:`~scipy.optimize.minimize` can then simply call it as
`partial_rastrigin(x)`, which finalizes the call to `rastrigin(x, a, b=b)`.

``partialwrap`` provides a convenience function :func:`~partialwrap.function_wrapper` that
generalises the above helper function `call_func_arg_kwarg` by passing all arguments, given as a
:any:`list`, and keyword arguments, given as a :any:`dict`, to arbitrary functions by the usual
`*args`, `**kwargs` mechanism:

.. code-block:: python

   from partialwrap import function_wrapper

   args   = [20.]
   kwargs = {'b': 1.*np.pi}
   rastra = partial(function_wrapper, rastrigin, args, kwargs)
   res    = opt.differential_evolution(rastra, [(xmin,xmax),]*ndim)

Note that you pass `args` and `kwargs` and not `*args` and `**kwargs` to :func:`~functools.partial`.
The wrapper is simply coded as the following and given for convenience:

.. code-block:: python

   def function_wrapper(func, arg, kwarg, x):
       return func(x, *arg, **kwarg)

Another example where partialisation might come in handy is the use of several CPUs to speed up
multiple evaluations of slow function evaluations. Parallelisation always has on overhead but it
will be beneficial in case of more expensive models, which is most often true for real research
problems. One can use the :func:`~multiprocessing.pool.Pool.map` function of Python's
:mod:`multiprocessing` module. I have 4 processors and can hence evaluate the function 4 times
simultaneously:

.. code-block:: python

   from multiprocessing import Pool

   neval = 100
   x = xmin + (xmax - xmin) * np.random.random_sample((neval,ndim))
   with Pool(4) as pool: 
       y = np.array(pool.map(rastra, x))

Note that all programming guidelines of the :mod:`multiprocessing` module apply
(https://docs.python.org/3/library/multiprocessing.html#programming-guidelines), which might
provide some "gotchas" if multiprocessing is throwing, for example, :any:`RuntimeError` or
:any:`AttributeError`.


Masked parameters in Python functions
=====================================

A common case in numerical optimization is the exclusion of some well-known parameters from
optimization, or fixing correlated parameters during optimization. But the numerical model still
needs to get a parameter value for the excluded/fixed parameters during optimization.
``partialwrap`` provides a convenience function :func:`~partialwrap.function_mask_wrapper` to
include only the masked parameters in the function evaluation and take default values where
`mask==False`:

.. code-block:: python

   from partialwrap import function_mask_wrapper

   x0      = np.array([0.5, 0.0001, 0.5])
   # Do not optimize the second parameter but take its initial value 0.0001
   mask    = [True, False, True]
   mrastra = partial(function_mask_wrapper, rastrigin, x0, mask, args, kwargs)

   res        = opt.differential_evolution(mrastra, [(xmin,xmax),]*np.sum(mask))
   xout       = x0.copy()
   xout[mask] = res.x

The values of `x0` will be taken where `mask==False`, i.e. `mask` could be called an include-mask.

The wrapper is very similar to :func:`~partialwrap.function_wrapper` in that it passes `args` and
`kwargs` to the function, but further needs the default values `x0` and the `mask` as inputs. The
wrapper is simply coded as the following and given for convenience:

.. code-block:: python

   def function_mask_wrapper(func, x0, mask, arg, kwarg, x):
       xx       = np.copy(x0)
       xx[mask] = x
       return func(xx, *arg, **kwarg)


External executables
====================

The great power of ``partialwrap`` is its ability to wrap external executables that cannot directly
be called from Python via :mod:`Cython` or :mod:`numpy.f2py` or similar.

``partialwrap`` provides two wrapper functions to work with external executables:
:func:`partialwrap.exe_wrapper` and :func:`partialwrap.exe_mask_wrapper`. The two wrappers
basically launch the external executable `exe` using Python's :mod:`subprocess` module, while
providing functionality to read and write parameter values and model output. The wrappers write a
parameter set into file(s) `parameterfile` that can be read by the external program `exe`. The
external program `exe` should write its result(s) to (a) file(s) `outputfile`, which will then be
read by the wrappers in return. This means that the two wrappers need to know a function
`parameterwriter` that writes the parameters in the file(s) `parameterfile` suitable for the
external model `exe`. The wrappers also need to know a function `outputreader` that reads the model
output(s) from the file(s) `outputfile`, and possibly calculating an objective value or just
passing back the output value(s).

Consider for simplicity an external Python program (e.g. `rastrigin1.py`)
that calculates the Rastrigin function with :math:`a = 10` and :math:`b = 2 \pi`,
reading in an arbitrary number of parameters :math:`x_i` from a
`parameterfile = params.txt` and writing its output into an
`outputfile = out.txt`:

.. code-block:: python

   # File: rastrigin1.py

   # Rastrigin function a=10, b=2*pi
   import numpy as np
   def rastrigin1(x):
       return 10.*len(x) + np.sum(x**2 - 10.*np.cos(2.*np.pi*x))

   # read parameters
   from partialwrap import standard_parameter_reader
   x = standard_parameter_reader('params.txt')

   # calc function
   y = rastrigin1(x)

   # write output file
   with open('out.txt', 'w') as ff:
       print(y, file=ff)

:func:`partialwrap.standard_parameter_reader` is a convenience functions that reads one parameter
per line from a file without a header.

The external program, which is in full `python3 rastrigin1.py`, can be used with the wrapper
function :func:`partialwrap.exe_wrapper` of ``partialwrap``:

.. code-block:: python

   from functools import partial
   from partialwrap import exe_wrapper, standard_parameter_writer, standard_output_reader
	
   rastrigin_exe  = ['python3', 'rastrigin1.py']
   parameterfile  = 'params.txt'
   outputfile     = 'out.txt'
   rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                            parameterfile, standard_parameter_writer,
                            outputfile, standard_output_reader, {})
   x0  = [0.1, 0.2]
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

:func:`partialwrap.standard_parameter_writer` is another convenience function that writes one
parameter per line in a file without a header. The function
:func:`partialwrap.standard_output_reader` simply reads one value from a file without a header. The
empty dictionary at the end of the partial statement is explained below.

Here I changed from the Differential Evolution algorithm to the quasi-Newton method of Broyden,
Fletcher, Goldfarb, and Shanno (BFGS). Differential Evolution finds the minimum of the Rastrigin
function much better than the gradient-based methods of :func:`scipy.optimize.minimize` but needs
much more function evaluations. It needs more that 2000 function evaluations for the
two-dimensional Rastrigin function, while BFGS stops after about 40-50 function evaluations with
the given good initial value. Calling the subprocess `python3 rastrigin1.py` needs about 180 ms on
my machine. The wrapper adds about another 50 ms for writing parameters, reading output, etc.,
so that the above example takes more than 9 minutes to finish with Differential
Evolution and about 7 seconds with BFGS.

One can see that the external Rastrigin program could have been written in C or Fortran or similar,
compiled, and then used with the :mod:`scipy.optimize` algorithms in Python. A Fortran program
could look like this:

.. code-block:: fortran

   program rastrigin1

       implicit none

       integer, parameter :: dp = kind(1.0d0)

       real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
       real(dp), parameter :: a  = 10.0_dp
       real(dp), parameter :: b  = 2.0_dp * pi

       character(len=*), parameter :: pfile = 'params.txt'
       character(len=*), parameter :: ofile = 'out.txt'

       integer, parameter :: punit = 99
       integer, parameter :: ounit = 101

       real(dp), dimension(100) :: x ! parameters, up to 100 dimensions
       real(dp) :: out               ! output value
       integer  :: n                 ! number of dimensions

       integer  :: ios

       ! read parameters
       open(punit, file=pfile, status='old', action='read')
       ios = 0
       n   = 1
       do while (ios==0)
           read(punit, fmt=*, iostat=ios) x(n)
           n = n + 1
       end do
       n = n - 2
       close(punit)

       ! calc function
       out = a * real(n,dp) + sum(x(1:n)**2 - a*cos(b*x(1:n)))

       ! write output file
       open(ounit, file=ofile)
       write(ounit,*) out
       close(ounit)

   end program rastrigin1

This program can be compiled like:

.. code-block:: bash

   gfortran -O3 -o rastrigin1.exe rastrigin1.f90

and called in Python:

.. code-block:: python

   rastrigin_f90exe  = ['./rastrigin1.exe']
   rastrigin_f90wrap = partial(exe_wrapper, rastrigin_f90exe,
                               parameterfile, standard_parameter_writer,
                               outputfile, standard_output_reader, {})
   x0  = [0.1, 0.2]
   res = opt.minimize(rastrigin_f90wrap, x0, method='BFGS')

The compiled Fortran needs about 5 ms when run in a :mod:`subprocess`, which is about one tenth of
the overhead of the wrapper function :func:`partialwrap.exe_wrapper`. Very fast executables can
hence be minimized in several seconds using ``partialwrap`` (about 4 seconds in the case of the
two-dimensional Rastrigin function) and the gradient-based methods of :mod:`scipy`, or in several
minutes with global search algorithms such as Differential Evolution (1.5 minutes in case of the
two-dimensional Rastrigin function).


Masked parameters with external executables
===========================================

Excluding parameters from, for example, optimization works exactly the same as for Python
functions. One passes the same arguments to :func:`partialwrap.exe_mask_wrapper` than to
:func:`partialwrap.exe_wrapper` plus the default values `x0` and the `mask`:

.. code-block:: python

   from partialwrap import exe_mask_wrapper

   rastrigin_f90exe  = ['./rastrigin1.exe']
   x0   = np.array([0.1, 0.0001, 0.2])
   mask = [True, False, True]
   mrastrigin_f90wrap = partial(exe_mask_wrapper, rastrigin_f90exe, x0, mask,
                                parameterfile, standard_parameter_writer,
                                outputfile, standard_output_reader, {})
   res = opt.minimize(mrastrigin_f90wrap, x0[mask], method='BFGS')
   xout       = x0.copy()
   xout[mask] = res.x

:func:`partialwrap.exe_mask_wrapper` basically does the transformation:

.. code-block:: python

   xx       = np.copy(x0)
   xx[mask] = x

and then calls :func:`~partialwrap.exe_wrapper` with `xx` (instead of `x`). So everything written
in the following about :func:`partialwrap.exe_wrapper` is also valid for
:func:`partialwrap.exe_mask_wrapper`.


Additional arguments for exe_wrapper
====================================

The user can pass further arguments to :func:`~partialwrap.exe_wrapper` via a dictionary at the end
of the call, which was empty at the examples above.

If you need to access shell features such as pipes, wildcards, environment variables, etc., the
external executable `exe` can be called in a shell. Setting the key `shell` to `True` passes
`shell=True` to :func:`subprocess.check_output`, executing the external executable `exe` in a
shell. Note that the `exe` name in :any:`subprocess` must be a string if `shell=True` and a
sequence if `shell=False`. Setting the key `debug` to `True` uses :func:`subprocess.check_call`
instead of :func:`subprocess.check_output` so that any output of the external executable will be
written to the screen (precisely :any:`subprocess.STDOUT`). This especially prints out also any
errors that might occur during execution. The above example using the external python program
`rastrigin1.py` can be debugged as:

.. code-block:: python

   from functools import partial
   from partialwrap import exe_wrapper, standard_parameter_writer, standard_output_reader
	
   rastrigin_exe  = 'python3 rastrigin1.py'
   parameterfile  = 'params.txt'
   outputfile     = 'out.txt'
   rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                            parameterfile, standard_parameter_writer,
                            outputfile, standard_output_reader,
                            {'shell':True, 'debug':True})
   x0  = [0.1, 0.2]
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

Note the change of `rastrigin_exe = ['python3', 'rastrigin1.py']` to `rastrigin_exe = 'python3
rastrigin1.py'` due to the use of `shell=True`.

Both, `parameterfile` and `outputfile` can either be single filenames (string) or a list of
filenames, which will be passed to `parameterwriter` and `outputreader`, respectively.
:func:`~partialwrap.exe_wrapper` deletes the parameter and output files after use. If you want to
keep the files, you can set the keys `keepparameterfile` and `keepoutputfile` to `True`. This can
be useful, for example, if your `parameterwriter` just changes a parameterfile in-place. An example
of such a `parameterwriter` is :func:`~partialwrap.sub_params_names`, which
substitutes all lines `name=.*` with `name=parameter` in the input files. The input file might be a
`parameterfile` for the external executable `exe`, where a parameter is given as
`parameter_name=parameter_value`, for example a Fortran namelist or a file in Python's standard
:mod:`configparser` format. If the first iteration of, for example, an optimization removed the
file, the next iteration could not use it again to insert the new parameter set. The
parameterwriter :func:`~partialwrap.sub_params_names` not only needs the filename(s)
`parameterfile` and the parameter values `params` as input as the above
:func:`standard_parameter_writer(parameterfile, params)` but also the `names` of the parameters. One can
pass additionally arguments `pargs` and keyword arguments `pkwargs` to the `parameterwriter` by
passing the dictionary entries `'pargs':parameterwriter_arguments` and
`'pkwargs':parameterwriter_keywords` to :func:`~partialwrap.exe_wrapper`.

Let's change the above external Python program `rastrigin1.py`, calling it `rastrigin2.py`, so that
it reads its parameters from an input file of the form `name = parameter`.

.. code-block::

   # File: params.txt
     param01 = 0.1 ! Fortran comment
   param03   = 0.3 # Python comment
    param02 = 0.2  // C comment

.. code-block:: python

   # File: rastrigin2.py

   # Rastrigin function a=10, b=2*pi
   import numpy as np
   def rastrigin1(x):
       return 10.*len(x) + np.sum(x**2 - 10.*np.cos(2.*np.pi*x))

   # read parameters
   with open('params.txt', 'r') as fi:
       pdict = {}
       for line in fi:
           ll = line.split()
           if (len(ll)==0) or ll[0].startswith('#'): continue
           pdict[ll[0]] = float(ll[2])
   x = np.array([ pdict[kk] for kk in sorted(pdict.keys()) ])

   # calc function
   y = rastrigin1(x)

   # write output file
   with open('out.txt', 'w') as ff:
       print(y, file=ff)

The parameterwriter :func:`~partialwrap.sub_params_names` will take the
`parameterfile='params.txt'`, searches for the lines that have nothing but whitespace before the
`names=['param01','param02','param03'] and replaces the lines with `names[i] = params[i]`.
`params.txt` will be reused for each iteration during optimization so should not be deleted by
:func:`~partialwrap.exe_wrapper`:

.. code-block:: python

   from functools import partial
   from partialwrap import exe_wrapper, sub_params_names, standard_output_reader
	
   rastrigin_exe  = ['python3', 'rastrigin2.py']
   parameterfile  = 'params.txt'
   outputfile     = 'out.txt'
   x0    = [ 0.1,       0.2,       0.5]
   names = ['param01', 'param02', 'param03']
   rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                            parameterfile, sub_params_names,
                            outputfile, standard_output_reader,
                            {'pargs':[names], 'keepparameterfile':True})
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

Note the list in `'pargs':[names]`. If one put `'pargs':names` than the `*args` mechanism would pass
three single arguments to `sub_params_names`, which would hence wrongly receive 5 instead of 3 arguments.

The same `*args/**kwargs` mechanism is implemented for the `outputreader`, where one can set the
keys `oargs` and `okwargs` to be passed to `outputreader`. This can be used, for example, to pass
observational data and uncertainty to calculate an evaluation metric such as a log-likelihood from
model output.


Provided parameterwriter and outputreader
=========================================

`partialwrap` comes with a few predefined `parameterwriter` and `outputreader`. The most basic ones
were used in the examples above. :func:`~partialwrap.standard_parameter_writer` simply writes
parameter values one by one in a file. :func:`~partialwrap.standard_parameter_reader` reads
parameters line by line from a file that does not contain any header lines, comments or similar.
:func:`~partialwrap.standard_output_reader` similarly reads a single value from a file.

:func:`~partialwrap.standard_parameter_writer_bounds_mask` writes another common format, which
includes one header line (# value min max mask) plus one line per parameter with the following
columns: consecutive parameter number, current parameter value, lower bound of parameter, upper
bound of parameter, 0/1 mask. :func:`~partialwrap.standard_parameter_reader_bounds_mask` reads
exactly these kind of files (lines starting with '#' will be ignored). A last example of an
`outputreader` is :func:`~partialwrap.standard_timeseries_reader` (or
:func:`~partialwrap.standard_time_series_reader`), which reads all lines from an output file into a
`numpy.ndarray`.

These `parameterwriter` and `outputreader` are given rather as examples for users to write their
own readers and writers.

A versatile `parameterwriter` ready to use is, however, :func:`~partialwrap.sub_params_names`,
which was used in the example above. It searches the `parameterfile` for lines that have nothing
but whitespace before given `names` and replaces the right hand side of the equal sign with the
parameter value. This can be used with a large variety of parameter files such as Python's
:mod:`configparser` files or Fortran namelists or similar. It exists in two variants: `names` are
case-sensitive in :func:`~partialwrap.sub_params_names_case` and case-insensitive in
:func:`~partialwrap.sub_params_names_ignorecase`. :func:`~partialwrap.sub_params_names` is simply a
wrapper for the latter case-insensitive function.

Another versatile `parameterwriter` that comes with `partialwrap` is
:func:`~partialwrap.sub_params_ja`. It searches for the strings #JA0000#, #JA0001#, ... in the
`parameterfile` and replaces them with the values of the first parameter, the second parameter, and
so on. The file must be well prepared in advance but the parameters can then be anywhere in the
`parameterfile`, appear several times on the same line or on different lines, etc. After
:func:`~partialwrap.sub_params_ja` was called once, for example by an optimization routine, the
tags #JA0000#, #JA0001#, ... would be gone and a second iteration could not fill in new parameter
values. The best way to use :func:`~partialwrap.sub_params_ja` is thus with the key `'pid':True` to
:func:`~partialwrap.exe_wrapper`, which will be explained in the next section about concurrent
execution of the external program.


Parallel evaluation of external executables
===========================================

Most real-life numerical models have longer run times than just a few milliseconds. One might
hence like to take advantage of more processing units such as simple multi-core processors,
multi-processor nodes or computer clusters. Take the simple parallel evaluation of the Rastrigin
function from above:

.. code-block:: python

   from functools import partial
   from multiprocessing import Pool
   from partialwrap import function_wrapper

   args   = [10.]
   kwargs = {'b': 2.*np.pi}
   rastra = partial(function_wrapper, rastrigin, args, kwargs)

   ndim  = 2
   xmin  = -5.12
   xmax  =  5.12
   neval = 100
   x = xmin + (xmax - xmin) * np.random.random_sample((neval,ndim))
   with Pool(4) as pool: 
       y = np.array(pool.map(rastra, x))

If we want to use the external program `rastrigin1.py` instead of the Python function
:func:`rastrigin`, one would naively do:

.. code-block:: python

   from functools import partial
   from multiprocessing import Pool
   from partialwrap import exe_wrapper, standard_parameter_writer, standard_output_reader

   rastrigin_exe  = ['python3', 'rastrigin1.py']
   parameterfile  = 'params.txt'
   outputfile     = 'out.txt'
   rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                            parameterfile, standard_parameter_writer,
                            outputfile, standard_output_reader, {})
   ndim  = 2
   xmin  = -5.12
   xmax  =  5.12
   neval = 100
   x = xmin + (xmax - xmin) * np.random.random_sample((neval,ndim))
   with Pool(4) as pool: 
       y = np.array(pool.map(rastrigin_wrap, x))

Python would fork 4 times and start 4 concurrent runs of :func:`rastrigin_wrap`. All 4 runs
would write the file 'params.txt', overwriting each other. ``partialwrap`` provides hence the key
`'pid':True` to :func:`~partialwrap.exe_wrapper`, which then passes a unique process identifier
(`pid`) as a keyword to the `parameterwriter`, adds the `pid` to the function call, and then
also passes the `pid` as a keyword to `outputreader`. :func:`~partialwrap.exe_wrapper` normally
deletes `parameterfile` and `outputfile` at its end. In case of `'pid':True`, it deletes
`parameterfile` and `outputfile` suffixed with `.pid`, if present.

So all `parameterwriter` provided by ``partialwrap`` take a keyword `pid` and, if present, write
`parameterfile.pid` rather than simply `parameterfile`. Likewise all `outputreader` provided by
``partialwrap`` take a keyword `pid` and if present read `outputfile.pid` rather than `outputfile`.
:func:`~partialwrap.exe_wrapper` then launches `exe+[str(pid)]` (or `exe+' '+str(pid) in case of
`'shell':True`). The external executable has hence to be able to read the `pid` from the command line
and, if present, read `parameterfile.pid` instead of `parameterfile` and write `outputfile.pid`
instead of `outputfile`. This can be handled with shell scripts if you are unable to change the
external model code (see below).

First, let's change `rastrigin1.py` so that it checks for command line input and uses it for the
`parameterfile` and `outputfile`:

.. code-block:: python

   # File: rastrigin3.py

   # get pid
   import sys
   if len(sys.argv) > 1:
       pid = int(sys.argv[1])
   else:
       pid = None

   # Rastrigin function a=10, b=2*pi
   import numpy as np
   def rastrigin1(x):
       return 10.*len(x) + np.sum(x**2 - 10.*np.cos(2.*np.pi*x))

   # read parameters
   from partialwrap import standard_parameter_reader
   x = standard_parameter_reader('params.txt', pid=pid)

   # calc function
   y = rastrigin1(x)

   # write output file
   if pid:
       fname = 'out.txt'+'.'+str(pid)
   else:
       fname = 'out.txt'
   with open(fname, 'w') as ff:
       print(y, file=ff)

Using `rastrigin3.py` with the key `'pid':True` would now evaluate four times the Rastrigin
function in parallel, every function evaluation using its individual parameter file:

.. code-block:: python

   from functools import partial
   from multiprocessing import Pool
   from partialwrap import exe_wrapper, standard_parameter_writer, standard_output_reader

   rastrigin_exe  = ['python3', 'rastrigin3.py']
   parameterfile  = 'params.txt'
   outputfile     = 'out.txt'
   rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                            parameterfile, standard_parameter_writer,
                            outputfile, standard_output_reader,
                            {'pid':True})

   ndim  = 2
   xmin  = -5.12
   xmax  =  5.12
   neval = 100
   x = xmin + (xmax - xmin) * np.random.random_sample((neval,ndim))
   with Pool(4) as pool: 
       y = np.array(pool.map(rastrigin_wrap, x))

One can take advantage of the workers keyword in :func:`scipy.optimize.differential_evolution` to
use available CPUs to find the minimum of the Rastrigin function:

.. code-block:: python

   from functools import partial
   from partialwrap import exe_wrapper, standard_parameter_writer, standard_output_reader

   rastrigin_exe  = ['python3', 'rastrigin3.py']
   parameterfile  = 'params.txt'
   outputfile     = 'out.txt'
   rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                            parameterfile, standard_parameter_writer,
                            outputfile, standard_output_reader,
                            {'pid':True})

   ndim  = 2
   xmin  = -5.12
   xmax  =  5.12
   res = opt.differential_evolution(rastrigin_wrap, [(xmin,xmax),]*ndim, workers=4)

Or one could use the popular :mod:`emcee` library to calculate parameter uncertainties with the
Markov chain Monte Carlo (MCMC) method. We take the example from the section on parallelization of
the :mod:`emcee` documentation (https://emcee.readthedocs.io/en/stable/tutorials/parallel/) but
code the log-likelihood function as an external Python program:

.. code-block:: python

   # File: logli1.py
   import numpy as np

   # get pid
   import sys
   if len(sys.argv) > 1:
       pid = int(sys.argv[1])
   else:
       pid = None

   # log-likelihood
   def log_prob(theta):
       return -0.5 * np.sum(theta**2)

   # read parameters
   from partialwrap import standard_parameter_reader
   x = standard_parameter_reader('params.txt', pid=pid)

   # calc function
   y = log_prob(x)

   # write output file
   if pid:
       fname = 'out.txt'+'.'+str(pid)
   else:
       fname = 'out.txt'
   with open(fname, 'w') as ff:
       print(y, file=ff)

Partialize it and sample the log-likelihood with :mod:`emcee` using a single processor:

.. code-block:: python

   from functools import partial
   from partialwrap import exe_wrapper, standard_parameter_writer, standard_output_reader

   logli_exe     = ['python3', 'logli1.py']
   parameterfile = 'params.txt'
   outputfile    = 'out.txt'
   logli_wrap    = partial(exe_wrapper, logli_exe,
                           parameterfile, standard_parameter_writer,
                           outputfile, standard_output_reader,
                           {'pid':True})

   # MCMC
   import time
   import emcee
   np.random.seed(42)
   initial = np.random.randn(32, 5)
   nwalkers, ndim = initial.shape
   nsteps  = 8

   sampler = emcee.EnsembleSampler(nwalkers, ndim, logli_wrap)
   start   = time.time()
   sampler.run_mcmc(initial, nsteps, progress=True)
   end     = time.time()
   serial_time = end - start
   print("Serial took {0:.1f} seconds".format(serial_time))

This takes about 80 seconds on my machine. The parallel version using Python's
:mod:`multiprocessing` module is:

.. code-block:: python

   import os
   os.environ["OMP_NUM_THREADS"] = "1"

   from multiprocessing import Pool
   with Pool() as pool:
       sampler = emcee.EnsembleSampler(nwalkers, ndim, logli_wrap, pool=pool)
       start   = time.time()
       sampler.run_mcmc(initial, nsteps, progress=True)
       end     = time.time()
       multi_time = end - start
       print("Multiprocessing took {0:.1f} seconds".format(multi_time))

This needs about 26 seconds on my machine; about 3 times faster.

One can see that the `parameterwriter` :func:`~partialwrap.sub_params_ja` works well with the 'pid'
key. The user prepares the `parameterfile` with the tags #JA0000#, #JA0001#, ....
:func:`~partialwrap.sub_params_ja` then takes `parameterfile` and writes the file
`parameterfile.pid` with the tags replaced by parameter values. No file is overwritten and
`parameterfile` can be reused by the next iteration or from a parallel process.


Using a launch script for the external program
----------------------------------------------

If one cannot change the external program to use a process identifier `pid` from the command line,
one can use a launch script that deals with `pid` by creating individual directories for each model
run and moving and renaming `parameterfile` and `outputfile`. The program `rastrigini1.py`, which
has no `pid` ability, could still be used using a bash script on Unix/Linux systems:

.. code-block:: bash

   # File: rastrigin1.sh

   #!/bin/bash

   set -e

   # get pid
   pid=${1}

   exe=rastrigin1.py
   pfile=params.txt
   ofile=out.txt

   # make individual run directory
   rundir=tmp.${pid}
   mkdir ${rundir}

   # copy individual parameter file
   mv ${pfile}.${pid} ${rundir}/${pfile}

   # run in individual directory
   cd ${rundir}
   ln -s ../${exe}
   python3 ${exe}

   # individualize output file
   mv ${ofile} ../${ofile}.${pid}

   # clean up
   cd ..
   rm -r ${rundir}

This would be used with :func:`~partialwrap.exe_wrapper` like:

.. code-block:: python

   from functools import partial
   from partialwrap import exe_wrapper, standard_parameter_writer, standard_output_reader

   rastrigin_exe  = ['./rastrigin1.sh']
   parameterfile  = 'params.txt'
   outputfile     = 'out.txt'
   rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                            parameterfile, standard_parameter_writer,
                            outputfile, standard_output_reader,
                            {'pid':True})

   x0  = [0.1, 0.2]
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

The bash script could, of course, also be a Python script to work on Windows platforms:

.. code-block:: python

   # File: run_rastrigin1.py

   # get pid
   import sys
   if len(sys.argv) > 1:
       pid = sys.argv[1]
   else:
       raise IOError('This scripts needs a process identifier (pid) as command line argument.')

   import os
   import shutil
   import subprocess

   exe   = 'rastrigin1.py'
   pfile = 'params.txt'
   ofile = 'out.txt'

   # make individual run directory
   rundir = 'tmp.'+pid
   os.mkdir(rundir)

   # copy individual parameter file
   os.rename(pfile+'.'+pid, rundir+'/'+pfile)

   # run in individual directory
   shutil.copyfile(exe, rundir+'/'+exe)
   os.chdir(rundir)
   err = subprocess.check_output(['python3', exe], stderr=subprocess.STDOUT)

   # make output available to exe_wrapper
   os.rename(ofile, '../'+ofile+'.'+pid)

   # clean up
   os.chdir('..')
   shutil.rmtree(rundir)

This would be used with :func:`~partialwrap.exe_wrapper` like:

.. code-block:: python

   from functools import partial
   from partialwrap import exe_wrapper, standard_parameter_writer, standard_output_reader

   rastrigin_exe  = ['python3', 'run_rastrigin1.py']
   parameterfile  = 'params.txt'
   outputfile     = 'out.txt'
   rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                            parameterfile, standard_parameter_writer,
                            outputfile, standard_output_reader,
                            {'pid':True})

   x0  = [0.1, 0.2]
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

That's all Folks!
