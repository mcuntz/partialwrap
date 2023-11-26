**********
User Guide
**********

Partial's ingenious mechanism allows to use even very complex
functions with many arguments and keyword arguments with routines that
need functions in the simple form `func(x)`. This includes Python's
:func:`map` function, the minimizers in :mod:`scipy.optimize`, and
plenty of third-party modules such as `emcee`_ or :mod:`pyeee`. It
also allows easy parallelization of code with, for example, the
parallel `map` function of the :mod:`multiprocessing` module or via
the Message Passing Interface (MPI) :mod:`mpi4py`.

``partialwrap`` provides two easy wrappers so that
:func:`functools.partial` can be used with external programs as well
as Python functions. It provides further convenience functions to deal
with input/output of external programs. The user guide gives examples
how to use pretty much every external executable with any Python
program.


External executables
====================

The great power of ``partialwrap`` is its ability to wrap external
executables that cannot directly be called from Python via `Cython`_,
:mod:`numpy.f2py` or similar.

``partialwrap`` provides two wrapper functions to work with external
executables: :func:`~partialwrap.wrappers.exe_wrapper` and
:func:`~partialwrap.wrappers.exe_mask_wrapper`. The two wrappers
basically launch an external executable `exe` using Python's
:mod:`subprocess` module, while providing functionality to write
parameter files and read in model output. The wrappers write a
parameter set into file(s) `parameterfile` that are needed by the
external program `exe`. The external program `exe` should write its
result(s) to (a) file(s) `outputfile`, which will then be read by the
wrappers in return. This means that the two wrappers need to know a
function `parameterwriter` that writes the parameters in the file(s)
`parameterfile` suitable for the external model `exe`. The wrappers
also need to know a function `outputreader` that reads the model
output(s) from the file(s) `outputfile`, and possibly calculating an
objective value or just passing back the output value(s).

Consider the *Rastrigin function*
(https://en.wikipedia.org/wiki/Rastrigin_function), which is a popular
function for performance testing of optimization algorithms:
:math:`y = a \cdot n + \sum_i^n (x_i^2 - a \cos(b \cdot x_i))`:

.. code-block:: python

   import numpy as np
   def rastrigin(x, a, b=2. * np.pi):
       return a * x.size + np.sum(x**2 - a * np.cos(b * x))

The Rastrigin function has a global minimum of :math:`0` at all
:math:`x_i = 0`. :math:`a` influences mainly the depths of the (local
and global) minima, whereas :math:`b` influences mainly the size of
the minima. A common form uses :math:`a = 10` and :math:`b = 2\pi`.
The parameters :math:`x_i` should then be in the interval
:math:`[-5.12, 5.12]`.

Take an external program that calculates the Rastrigin function with
:math:`a = 10` and :math:`b = 2 \pi`, reading in an arbitrary number
of parameters :math:`x_i` from a `parameterfile = params.txt` and
writing its output into an `outputfile = out.txt`. Take for simplicity
a Python program first (e.g. `rastrigin.py`):

.. code-block:: python

   # File: rastrigin.py
   import numpy as np
   from partialwrap import standard_parameter_reader

   # Rastrigin function a=10, b=2*pi
   def rastrigin(x):
       return 10. * len(x) + np.sum(x**2 - 10. * np.cos(2. * np.pi * x))

   # read parameters
   x = standard_parameter_reader('params.txt')

   # calc function
   y = rastrigin(x)

   # write output file
   with open('out.txt', 'w') as ff:
       print(y, file=ff)

:func:`~partialwrap.std_io.standard_parameter_reader` is a convenience
functions that reads parameters from a file, one per line returning a
:class:`numpy.ndarray`. The external program, which is in full
`python3 rastrigin.py`, can be used with the wrapper function
:func:`~partialwrap.wrappers.exe_wrapper` of ``partialwrap``:

.. code-block:: python

   from functools import partial
   import scipy.optimize as opt
   from partialwrap import exe_wrapper
   from partialwrap import standard_parameter_writer, standard_output_reader
        
   rastrigin_exe   = ['python3', 'rastrigin.py']
   parameterfile   = 'params.txt'
   parameterwriter = standard_parameter_writer
   outputfile      = 'out.txt'
   outputreader    = standard_output_reader
   rastrigin_wrap  = partial(exe_wrapper, rastrigin_exe,
                             parameterfile, parameterwriter,
                             outputfile, outputreader, {})
   x0  = [0.1, 0.2, 0.3]
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

:func:`~partialwrap.std_io.standard_parameter_writer` is another
convenience function that writes parameters into a file, one per
line. The function :func:`~partialwrap.std_io.standard_output_reader`
simply reads one value from a file. The empty dictionary at the end of
the partial statement is explained below.

One can see that the external Rastrigin program could have been
written in a compiled language such as C, Fortran or similar, and then
used with the :mod:`scipy.optimize` algorithms in Python. A Fortran
program could look like this:

.. code-block:: fortran

   program rastrigin

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
       out = a * real(n,dp) + sum(x(1:n)**2 - a * cos(b * x(1:n)))

       ! write output file
       open(ounit, file=ofile)
       write(ounit,*) out
       close(ounit)

   end program rastrigin

This program can be compiled like:

.. code-block:: bash

   gfortran -o rastrigin.exe rastrigin.f90

and used in Python:

.. code-block:: python

   from functools import partial
   import scipy.optimize as opt
   from partialwrap import exe_wrapper
   from partialwrap import standard_parameter_writer, standard_output_reader
        
   rastrigin_exe   = ['./rastrigin.exe']
   parameterfile   = 'params.txt'
   parameterwriter = standard_parameter_writer
   outputfile      = 'out.txt'
   outputreader    = standard_output_reader
   rastrigin_wrap  = partial(exe_wrapper, rastrigin_exe,
                             parameterfile, parameterwriter,
                             outputfile, outputreader, {})
   x0  = [0.1, 0.2, 0.3]
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

Where the only difference to the Python version is that
`rastrigin_exe = ['./rastrigin.exe']` instead of
`rastrigin_exe = ['python3', 'rastrigin.py']`.

The wrapper function :func:`~partialwrap.wrappers.exe_wrapper` adds
about 50 ms of overhead per function evaluation. Compare this to 5 ms
of the compiled Fortran program when run in a :mod:`subprocess` and
180 ms for the Python version. Real-life executables normally have
much longer run times so that the overhead is negligible.


Masked parameters
=================

A common case in numerical optimization is the exclusion of some
well-known or screened parameters from optimization, or fixing
correlated parameters during optimization. But the numerical model
still needs to get a parameter value for the excluded/fixed parameters
during optimization. ``partialwrap`` provides the wrapper function
:func:`~partialwrap.wrappers.exe_mask_wrapper` to include only the
masked parameters in the function evaluation and take default values
`x0` where `mask==False`. Fixing the second parameter to a default
value of 0.0001 changes the above program to:

.. code-block:: python

   from functools import partial
   import numpy as np
   import scipy.optimize as opt
   from partialwrap import exe_mask_wrapper
   from partialwrap import standard_parameter_writer, standard_output_reader
        
   rastrigin_exe   = ['./rastrigin.exe']
   parameterfile   = 'params.txt'
   parameterwriter = standard_parameter_writer
   outputfile      = 'out.txt'
   outputreader    = standard_output_reader
   x0              = np.array([0.1, 0.0001, 0.2])
   mask            = [True, False, True]
   rastrigin_wrap  = partial(exe_mask_wrapper, rastrigin_exe, x0, mask,
                             parameterfile, parameterwriter,
                             outputfile, outputreader, {})
   res = opt.minimize(rastrigin_wrap, x0[mask], method='BFGS')
   xout       = x0.copy()
   xout[mask] = res.x

The values of `x0` will be taken where `mask==False`, i.e. `mask`
could be called an *include mask* in contrast to the mask of numpy's
`masked arrays`_, which is rather an *exclude mask*.

Note that the optimizer :func:`~scipy.optimize.minimize` 'sees' only
the masked parameters so that the initial values `x0` and possible
`bounds` given to the optimizer should also be only the masked values,
here `x0[mask]`.

:func:`~partialwrap.wrappers.exe_mask_wrapper` basically does the
transformation:

.. code-block:: python

   xx       = np.copy(x0)
   xx[mask] = x

and then calls :func:`~partialwrap.wrappers.exe_wrapper` with `xx`
(instead of `x`). So everything written in the following about
:func:`~partialwrap.wrappers.exe_wrapper` is also valid for
:func:`~partialwrap.wrappers.exe_mask_wrapper`.


Additional arguments
====================

The user can pass further arguments to
:func:`~partialwrap.wrappers.exe_wrapper` via a dictionary at the end
of the call, which was empty at the examples above.

**shell** If one needs to access shell features such as pipes,
wildcards, environment variables, etc., in the call to the external
executable `exe`, the latter can be called in a separate shell of the
operating system. Setting the key `shell` to `True` passes
`shell=True` to :func:`subprocess.run`. Note that the interpretation
of `exe` in :func:`subprocess.run` is dependent on the operating
system if `shell=True`. It must hence be a string if `shell=True` and
it should be a sequence (list or tuple) if `shell=False`.

**debug** Setting the key `debug` to `True` writes any output of the
external executable to the screen (precisely
:any:`subprocess.STDOUT`). This especially prints out any errors
that might occur during execution of the executable `exe`. The above
example using the external Python program `rastrigin.py` can be
debugged as:

.. code-block:: python

   from functools import partial
   import scipy.optimize as opt
   from partialwrap import exe_wrapper
   from partialwrap import standard_parameter_writer, standard_output_reader
        
   rastrigin_exe   = 'python3 rastrigin.py'
   parameterfile   = 'params.txt'
   parameterwriter = standard_parameter_writer
   outputfile      = 'out.txt'
   outputreader    = standard_output_reader
   rastrigin_wrap  = partial(exe_wrapper, rastrigin_exe,
                             parameterfile, parameterwriter,
                             outputfile, outputreader,
                             {'shell': True, 'debug': True})
   x0  = [0.1, 0.2, 0.3]
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

Note the change from `rastrigin_exe = ['python3', 'rastrigin.py']` at
the earlier example to `rastrigin_exe = 'python3 rastrigin.py'` here
due to the use of `shell=True`.

**keepparameterfile, keepoutputfile** Both, `parameterfile` and
`outputfile` can either be single filenames (string) or a list of
filenames, which will be passed to `parameterwriter` and
`outputreader`, respectively.
:func:`~partialwrap.wrappers.exe_wrapper` deletes the parameter and
output files on disk after use. If one wants to keep the files, one can
set the keys `keepparameterfile` and `keepoutputfile` to `True`. This
can be useful, for example, if the `parameterwriter` just changes a
parameterfile in-place. An example of such a `parameterwriter` is
:func:`~partialwrap.std_io.sub_params_names`, which substitutes all
lines `name = .*` with `name = parameter` in the parameter file(s). An
input to :func:`~partialwrap.std_io.sub_params_names` could be a
`parameterfile` for the external executable `exe` in which a parameter
is given as `parameter_name = parameter_value`. This is, for example,
the case in Fortran namelists or files for Python's standard
:mod:`configparser`. Imagine an optimization of the external
executable `exe` with such a `parameterwriter` and
:func:`~partialwrap.wrappers.exe_wrapper` deleting the parameter file
after use. In the first iteration, the `parameterwriter` would take
the `parameterfile` and change its content. The `exe` would be run,
the output read by `outputreader`, and then
:func:`~partialwrap.wrappers.exe_wrapper` would delete
`parameterfile`. So there would be no `parameterfile` file anymore in
the second iteration that `parameterwriter` could change. In this
case, one can set `keepparameterfile = True` and
:func:`~partialwrap.wrappers.exe_wrapper` would not delete the
`parameterfile` after use.

**pargs, pkwargs** The parameterwriter
:func:`~partialwrap.std_io.sub_params_names` not only needs the
filename(s) `parameterfile` and the parameter values `params` but also
the `names` of the parameters. One can pass additionally arguments
`pargs` and keyword arguments `pkwargs` to the `parameterwriter` by
passing the dictionary entries `'pargs': parameterwriter_arguments`
and `'pkwargs': parameterwriter_keywords` to
:func:`~partialwrap.wrappers.exe_wrapper`.

Let's change the above external Python program `rastrigin.py`,
calling it `rastrigin_config.py`, so that it reads its parameters from an
input file of the form `name = parameter`:

.. code-block::

   # File: params.txt
     param01 = 0.1 ! Fortran comment
   param03   = 0.3 # Python comment
    param02 = 0.2  // C comment

The parameter names are not aligned in the above example to
demonstrate that whitespace before the names on the right-hand-side
will be ignored but will be preserved in the changed parameter file.

.. code-block:: python

   # File: rastrigin_config.py
   import numpy as np

   # Rastrigin function a=10, b=2*pi
   def rastrigin(x):
       return 10. * len(x) + np.sum(x**2 - 10. * np.cos(2. * np.pi * x))

   # read parameters
   with open('params.txt', 'r') as fi:
       pdict = {}
       for line in fi:
           ll = line.split()
           if (len(ll) == 0) or ll[0].startswith('#'):
               continue
           pdict[ll[0]] = float(ll[2])
   x = np.array([ pdict[kk] for kk in sorted(pdict.keys()) ])

   # calc function
   y = rastrigin(x)

   # write output file
   with open('out.txt', 'w') as ff:
       print(y, file=ff)

The parameterwriter :func:`~partialwrap.std_io.sub_params_names` will
take the `parameterfile = 'params.txt'`, searches for the lines that
have nothing but whitespace before the `names = ['param01', 'param02',
'param03']` and replaces the lines with `names[i] = params[i]` for i
in range(len(params)). `params.txt` will be reused for each iteration
during optimization so should not be deleted by
:func:`~partialwrap.wrappers.exe_wrapper`:

.. code-block:: python

   from functools import partial
   import scipy.optimize as opt
   from partialwrap import exe_wrapper
   from partialwrap import sub_params_names, standard_output_reader
        
   rastrigin_exe   = ['python3', 'rastrigin_config.py']
   parameterfile   = 'params.txt'
   parameterwriter = sub_params_names
   outputfile      = 'out.txt'
   outputreader    = standard_output_reader
   x0              = [0.1, 0.2, 0.3]
   names           = ['param01', 'param02', 'param03']
   rastrigin_wrap  = partial(exe_wrapper, rastrigin_exe,
                             parameterfile, parameterwriter,
                             outputfile, outputreader,
                             {'pargs': [names], 'keepparameterfile': True})
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

Note `pargs` is given a list `'pargs': [names]`, which gives `'pargs':
[['param01', 'param02', 'param03']]`. If one would simply put
`'pargs': names`, than the `*args` mechanism would expand the list
`names` to three individual arguments for `sub_params_names`, so that
the latter would receive wrongly 5 instead of 3 arguments.

**oargs, okwargs** The same `*args/**kwargs` mechanism is implemented
for the `outputreader`, where one can set the keys `oargs` and
`okwargs` to be passed to `outputreader`. This can be used, for
example, to pass observational data and uncertainty to calculate an
evaluation metric such as a log-likelihood from model output.

**pid** This keyword offers another possibility not to delete
`parameterfile` and `outputfile` after use. It is especially useful
with parallel evaluations of external executables `exe` (see
`Parallel evaluation`_ below). If `pid = True` is given,
:func:`~partialwrap.wrappers.exe_wrapper` will write the file
`parameterfile.pid` rather than `parameterfile`, with `pid` being a
random number. It will launch `exe + [str(pid)]` instead of `exe`, and
it will read `outputfile.pid` instead of `outputfile`. The external
executable `exe` has hence to be able to read the `pid` from the
command line, read `parameterfile.pid` instead of `parameterfile` and
write `outputfile.pid` instead of `outputfile`. This can be handled
with shell scripts if one is unable to change the external model code
(also see below).


Parallel evaluation
===================

Most real-life numerical models have longer run times than just a few
milliseconds. One might hence like to take advantage of more processing units
such as simple multi-core processors, multi-processor nodes or computer
clusters.

One may want to find the minimum of the Rastrigin function in
`rastrigin.py` taking advantage of multiple processors. Scipy's
:func:`scipy.optimize.differential_evolution` has the *workers*
keyword to use more than one CPU to find the minimum of a
function. One would naively do:

.. code-block:: python

   from functools import partial
   import scipy.optimize as opt
   from partialwrap import exe_wrapper
   from partialwrap import standard_parameter_writer, standard_output_reader
        
   rastrigin_exe   = ['python3', 'rastrigin.py']
   parameterfile   = 'params.txt'
   parameterwriter = standard_parameter_writer
   outputfile      = 'out.txt'
   outputreader    = standard_output_reader
   rastrigin_wrap  = partial(exe_wrapper, rastrigin_exe,
                             parameterfile, parameterwriter,
                             outputfile, outputreader, {})
   ndim = 3
   bounds = [(-5.12, 5.12),] * ndim
   res = opt.differential_evolution(rastrigin_wrap, bounds, workers=4)

This subdivides the population into 4 worker sections and evaluates
them in parallel (using :func:`multiprocessing.Pool`). This means
`rastrigin.py` will be launched 4 times in parallel, each writing a
`parameterfile` *params.txt*, hence overwriting each other all the
time. Here the `pid` keyword comes in handy. Each invocation of `rastrigin.py` would have its own random number `pid` associated, writing `parameterfile.pid` and reading `outfile.pid`. The rastrigin program would need to be changed to:

.. code-block:: python

   # File: rastrigin_pid.py
   import sys
   import numpy as np
   from partialwrap import standard_parameter_reader

   # Rastrigin function a=10, b=2*pi
   def rastrigin(x):
       return 10. * len(x) + np.sum(x**2 - 10. * np.cos(2. * np.pi * x))

   # get pid
   if len(sys.argv) > 1:
       pid = int(sys.argv[1])
   else:
       pid = None

   # read parameters
   x = standard_parameter_reader('params.txt', pid=pid)

   # calc function
   y = rastrigin(x)

   # write output file
   if pid is None:
       fname = 'out.txt'
   else:
       fname = f'out.txt.{pid}'
   with open(fname, 'w') as ff:
       print(y, file=ff)

All input/output routines provided by ``partialwrap`` take a keyword
`pid` and, if present, read or write `file.pid` rather than `file`
with *pid* being a unique random number. So here
:func:`~partialwrap.std_io.standard_parameter_reader` reads from files
such as `params.txt.158398716` rather than from `params.txt`. Then simply the `pid` keyword has to be set to *True* in the call of `partial`:

.. code-block:: python

   from functools import partial
   import scipy.optimize as opt
   from partialwrap import exe_wrapper
   from partialwrap import standard_parameter_writer, standard_output_reader
        
   rastrigin_exe   = ['python3', 'rastrigin_pid.py']
   parameterfile   = 'params.txt'
   parameterwriter = standard_parameter_writer
   outputfile      = 'out.txt'
   outputreader    = standard_output_reader
   rastrigin_wrap  = partial(exe_wrapper, rastrigin_exe,
                             parameterfile, parameterwriter,
                             outputfile, outputreader, {'pid': True})
   ndim = 3
   bounds = [(-5.12, 5.12),] * ndim
   res = opt.differential_evolution(rastrigin_wrap, bounds, workers=4)
   
Another example could use the popular `emcee`_ library to calculate
parameter uncertainties with the Markov chain Monte Carlo (MCMC)
method. We take the example from the section on parallelization of the
`emcee`_ documentation
(https://emcee.readthedocs.io/en/stable/tutorials/parallel/) but code
the log-likelihood function as an external Python program:

.. code-block:: python

   # File: logli.py
   import sys
   import numpy as np
   from partialwrap import standard_parameter_reader

   # log-likelihood
   def log_prob(theta):
       return -0.5 * np.sum(theta**2)

   # get pid
   if len(sys.argv) > 1:
       pid = int(sys.argv[1])
   else:
       pid = None

   # read parameters
   x = standard_parameter_reader('params.txt', pid=pid)

   # calc function
   y = log_prob(x)

   # write output file
   if pid is None:
       fname = 'out.txt'
   else:
       fname = f'out.txt.{pid}'
   with open(fname, 'w') as ff:
       print(y, file=ff)

Partialize it and sample the log-likelihood with `emcee`_ using a single processor:

.. code-block:: python

   from functools import partial
   from partialwrap import exe_wrapper
   from partialwrap import standard_parameter_writer, standard_output_reader
   import emcee

   logli_exe     = ['python3', 'logli.py']
   parameterfile   = 'params.txt'
   parameterwriter = standard_parameter_writer
   outputfile      = 'out.txt'
   outputreader    = standard_output_reader
   logli_wrap      = partial(exe_wrapper, logli_exe,
                             parameterfile, parameterwriter,
                             outputfile, outputreader, {'pid': True})
   # MCMC
   np.random.seed(42)
   initial = np.random.randn(32, 5)
   nwalkers, ndim = initial.shape
   nsteps  = 8

   # Single processor version
   sampler = emcee.EnsembleSampler(nwalkers, ndim, logli_wrap)
   sampler.run_mcmc(initial, nsteps, progress=True)

or use multiple processors by changing the last two lines to:

.. code-block:: python

   from multiprocessing import Pool
   with Pool() as pool:
       sampler = emcee.EnsembleSampler(nwalkers, ndim, logli_wrap, pool=pool)
       sampler.run_mcmc(initial, nsteps, progress=True)

The serial version takes about 40 seconds on my machine, while the parallel version takes about 6 seconds with 16 cores on my machine.


Using launch scripts
====================

If one cannot change the external program to use a process identifier
`pid` from the command line (i.e. using `rastrigin_pid.py`), one can use
a launch script that deals with `pid` by creating individual
directories for each model run and moving and renaming `parameterfile`
and `outputfile`. The program `rastrigini1.py`, which has no `pid`
ability, could still be used using a bash script on Unix/Linux and macOS:

.. code-block:: bash

   # File: rastrigin.sh
   #!/usr/bin/env bash
   set -e

   # get pid
   pid=${1}

   exe=rastrigin.py
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

This would be used with :func:`~partialwrap.wrappers.exe_wrapper`
in exactly the same way as above:

.. code-block:: python

   from functools import partial
   import scipy.optimize as opt
   from partialwrap import exe_wrapper
   from partialwrap import standard_parameter_writer, standard_output_reader
        
   rastrigin_exe   = ['rastrigin.sh']
   parameterfile   = 'params.txt'
   parameterwriter = standard_parameter_writer
   outputfile      = 'out.txt'
   outputreader    = standard_output_reader
   rastrigin_wrap  = partial(exe_wrapper, rastrigin_exe,
                             parameterfile, parameterwriter,
                             outputfile, outputreader, {'pid': True})
   ndim = 3
   bounds = [(-5.12, 5.12),] * ndim
   res = opt.differential_evolution(rastrigin_wrap, bounds, workers=4)

The bash script could, of course, also be a Python script to work
on Windows platforms as well:

.. code-block:: python

   # File: run_rastrigin.py
   import os
   import shutil
   import subprocess
   import sys

   # get pid
   if len(sys.argv) > 1:
       pid = sys.argv[1]
   else:
       pid = None

   exe   = 'rastrigin.py'
   pfile = 'params.txt'
   ofile = 'out.txt'

   # make individual run directory
   if pid is None:
       rundir = 'tmp'
   else:
       rundir = f'tmp.{pid}'
   os.mkdir(rundir)

   # copy individual parameter file
   if pid is None:
       os.rename(f'{pfile}', f'{rundir}/{pfile}')
   else:
       os.rename(f'{pfile}.{pid}', f'{rundir}/{pfile}')

   # run in individual directory
   shutil.copyfile(exe, f'{rundir}/{exe}')
   os.chdir(rundir)
   err = subprocess.check_output(['python3', exe],
                                 stderr=subprocess.STDOUT)

   # make output available to exe_wrapper
   if pid is None:
       os.rename(ofile, f'../{ofile}')
   else:
       os.rename(ofile, f'../{ofile}.{pid}')

   # clean up
   os.chdir('..')
   shutil.rmtree(rundir)

This Python script could be used exactly as the shell script above:

.. code-block:: python

   from functools import partial
   import scipy.optimize as opt
   from partialwrap import exe_wrapper
   from partialwrap import standard_parameter_writer, standard_output_reader
        
   rastrigin_exe   = ['python3', 'run_rastrigin.py']
   parameterfile   = 'params.txt'
   parameterwriter = standard_parameter_writer
   outputfile      = 'out.txt'
   outputreader    = standard_output_reader
   rastrigin_wrap  = partial(exe_wrapper, rastrigin_exe,
                             parameterfile, parameterwriter,
                             outputfile, outputreader, {'pid': True})
   ndim = 3
   bounds = [(-5.12, 5.12),] * ndim
   res = opt.differential_evolution(rastrigin_wrap, bounds, workers=4)


Input/Output functions
======================

``partialwrap`` comes with a few read and write routines for
parameters and output. All routines support the `pid` keyword,
i.e. take the `pid=something` keyword argument and use `something`, if
present, to suffix the input file(s) with it on output. The provided
input/output routines are seen as examples on which one can build its
own tailored input/output functions.

**standard_output_reader**
:func:`~partialwrap.std_io.standard_output_reader` reads a single
value from a file. For example:

    | 0.123456789

**standard_parameter_reader**
:func:`~partialwrap.std_io.standard_parameter_reader` reads parameters
from a file, one parameter per line. For example:

    | 0.1
    | 0.3
    | 0.2
    | ...

**standard_parameter_writer**
:func:`~partialwrap.std_io.standard_parameter_writer` writes one
parameter per line into a file. For example:

    | 1.000000000000000e-01
    | 3.000000000000000e-01
    | 2.000000000000000e-01
    | ...

**standard_time_series_reader, standard_timeseries_reader**
:func:`~partialwrap.std_io.standard_timeseries_reader` (or
:func:`~partialwrap.std_io.standard_time_series_reader`) reads all
lines from an output file into a :class:`numpy.ndarray`. For example:

    | 0.123456789
    | 0.234567890
    | 0.345678901
    | ...

**standard_parameter_reader_bounds_mask**
:func:`~partialwrap.std_io.standard_parameter_reader_bounds_mask`
reads a common parameter file format. It has one header line (# value
min max mask) plus one line per parameter with the following columns:
consecutive parameter number, current parameter value, lower bound of
parameter, upper bound of parameter, 0/1 mask. For example:

    | # value min max mask
    | 1 3.0e-01 0.0 1. 1
    | 2 2.3e-01 -1.0 1. 1
    | 3 1.44e+01 9.0 20. 1
    | 4 3.0e-01 0.0 1. 0
    | ...

**standard_parameter_writer_bounds_mask**
:func:`~partialwrap.std_io.standard_parameter_writer_bounds_mask` writes a
parameter file with values, bound and mask. It has one header line (#
value min max mask) plus one line per parameter with the following
columns: consecutive parameter number, current parameter value, lower
bound of parameter, upper bound of parameter, 0/1 mask. For example:

    | # value min max mask
    | 1 3.000000000000000e-01 0.000000000000000e+00 1.000000000000000e+00 1
    | 2 2.300000000000000e-01 -1.000000000000000e+00 1.000000000000000e+00 1
    | 3 1.440000000000000e+01 9.000000000000000e+00 2.000000000000000e+01 1
    | 4 3.000000000000000e-01 0.000000000000000e+00 1.000000000000000e+00 0
    | ...

**sub_params_names, sub_params_names_ignorecase**
:func:`~partialwrap.std_io.sub_params_names` (or
:func:`~partialwrap.std_io.sub_params_names_ignorecase`) substitutes
`name = .*` with `name = parameter` in input file(s); names are case
insensitive. It searches the input file(s) for lines that have nothing
but whitespace before given `names` and replaces the right hand side
of the equal sign with the parameter value. An input to
:func:`~partialwrap.std_io.sub_params_names` could be a
`parameterfile` in which a parameters are given as `parameter_name =
parameter_value`. This is, for example, the case in Fortran
namelists. For example:

    | &parameters
    |  param01 = 0.1 ! Fortran comment
    | Param03   = 0.3
    |  PaRaM02 = 0.2
    | /

**sub_params_names_case**
:func:`~partialwrap.std_io.sub_params_names_case` substitutes `name =
.*` with `name = parameter` in input file(s); names are case
sensitive. It is the same as
:func:`~partialwrap.std_io.sub_params_names_ignorecase` but `names`
are case sensitive. An example input file are files for Python's
standard :mod:`configparser`. For example:

    | [Parameters]
    |  param01 = 0.1
    | Param03   = 0.3 # Python comment
    |  PaRaM02 = 0.2

**sub_params_ja** :func:`~partialwrap.std_io.sub_params_ja`
substitutes the strings *#JA0000#*, *#JA0001#*, ... in the input
file(s) with the parameter values. It searches for the tags
*#JA0000#*, *#JA0001#*, ... in the input file(s) and replaces them
with the values of the first parameter, the second parameter, and so
on. The file must be prepared in advance and the parameters can then
be anywhere in the file(s), appear several times on the same line or
on different lines, etc. For example:

    | &parameters
    |  param01 = #JA0000#
    |  param02 = 0.3
    |  Param03 = #JA0001#, #JA0001#, #JA0001#
    | /

The substitution functions
:func:`~partialwrap.std_io.sub_params_names`,
:func:`~partialwrap.std_io.sub_params_names_ignorecase`,
:func:`~partialwrap.std_io.sub_params_names_case`, and
:func:`~partialwrap.std_io.sub_params_ja` change the input
file(s). :func:`~partialwrap.wrappers.exe_wrapper` deletes files after
use by default. The best way to use the substitution functions is
hence with `'keepparameterfile': True` or `'pid': True`. The last
function :func:`~partialwrap.std_io.sub_params_ja` substitutes the
tags *#JA0000#*, *#JA0001#*, ... in the input file(s). Once
substituted, the tags would not be in the input files anymore. Hence
this does not work with `'keepparameterfile': True` but rather works
with only with `'pid': True`.


Python functions
================

``partialwrap`` has also two convenience functions to do the same
mechanisms with Python functions:
:func:`~partialwrap.wrappers.function_wrapper` and
:func:`~partialwrap.wrappers.function_mask_wrapper`. They can, for
example, be used during development of a complex code as substitution
for the exe-equivalents.

Take the Python function *rastrigin*:

.. code-block:: python

   import numpy as np
   def rastrigin(x, a, b=2. ** np.pi):
       return a * x.size + np.sum(x**2 - a * np.cos(b * x))

It's minimum can be found with

.. code-block:: python

   import scipy.optimize as opt

   x0  = [0.1, 0.2, 0.3]
   res = opt.minimize(rastrigin, x0, method='BFGS')

The minimizer :func:`~scipy.optimize.minimize` can pass arguments to
the function *rastrigin* but does not has the functionality to pass
keyword parameters. It is hence *a priori* not possible to search the
minimum of the Rastrigin function with another :math:`b`
parameter. (This is an illustrative example while it is, of course,
possible to just pass `b` as a second argument in this case.) One
could use Python's :func:`functools.partial` in this case:

.. code-block:: python

   from functools import partial
   import scipy.optimize as opt

   # helper function
   def call_func_arg_kwarg(func, a, b, x):
       return func(x, a, b=b)

   # Partialise function with fixed parameters
   a = 5.
   b = 4. * np.pi
   partial_rastrigin = partial(call_func_arg_kwarg, rastrigin, a, b)

   x0  = [0.1, 0.2, 0.3]
   res = opt.minimize(partial_rastrigin, x0, method='BFGS')

``partialwrap`` provides a convenience function
:func:`~partialwrap.wrappers.function_wrapper` that generalises the
above helper function `call_func_arg_kwarg` by passing all arguments,
given as a :any:`list`, and keyword arguments, given as a :any:`dict`,
to arbitrary functions by the usual `*args`, `**kwargs` mechanism:

.. code-block:: python

   from functools import partial
   import scipy.optimize as opt
   from partialwrap import function_wrapper

   # Partialise function with fixed parameters
   args   = [5.]
   kwargs = {'b': 4. * np.pi}
   rastrigin_wrap = partial(function_wrapper, rastrigin, args, kwargs)

   x0  = [0.1, 0.2, 0.3]
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

Note that one passes `args` and `kwargs` and not `*args` and `**kwargs`
to :func:`~functools.partial`.

There is also the masked version
:func:`~partialwrap.wrappers.function_mask_wrapper` to exclude some
parameters from optimization. For example fixing the second parameter
to a default value of 0.0001 changes the above to:

.. code-block:: python

   from functools import partial
   import scipy.optimize as opt
   from partialwrap import function_mask_wrapper

   # Partialise function with fixed parameters
   args   = [5.]
   kwargs = {'b': 4. * np.pi}
   x0     = np.array([0.5, 0.0001, 0.5])
   # Do not optimize the second parameter but take its initial value 0.0001
   mask   = [True, False, True]
   rastrigin_wrap = partial(function_mask_wrapper, rastrigin, x0, mask,
                            args, kwargs)

   res = opt.minimize(rastrigin_wrap, x0[mask], method='BFGS')
   xout       = x0.copy()
   xout[mask] = res.x


A real life example
===================

We use the ecosystem model `MuSICA`_ in our research. The model is
written in Fortran and uses so-called namelists for configuration,
which look like this:

    | &namelist_01
    |  parameter01 = 3.141592653589793
    |  parameter02 = 1, 1, 2, 3, 5, 8, 13
    | /
    |
    | &namelist_02
    |  another_parameter01 = 2.718281828459045
    |  another_parameter02 = 'equilibrium'
    | /

The namelists are organized in three files: *musica.nml*,
*musica_soil.nml*, and one file for each tree species modelled, for
example, *fagus_sylvatica.nml*. There are about 50 parameters in the
namelists that influence the model output.

`MuSICA`_ calculates energy, water, and carbon fluxes in an ecosystem
such as a forest. One primary output is so-called Net Ecosystem
Exchange (NEE). One might want to adapt a few parameters, which were
not measured, for a specific forest stand by optimizing model output
against observed NEE. Below is a documented Python program that
optimizes three parameters for the description of stomatal conductance
in MuSICA: *gs_slope*, *gs_hx_half*, and *gs_hx_shape*. The first one
is the empirical slope in the description of stomatal conductance *gs*
(Ball-Berry *m* or Medlyn *g1*, for example) and the other two
variables are the parameters of a logistic function that links
stomatal conductance *gs* with the water potential in the xylem *hx*
(Tuzet model). Here I prepared the namelist file fagus_sylvatica.nml
and put tags *#JA0000#*, *#JA0001#*, and *#JA0002#* instead of values
on the right-hand-side of the equal sign and will then use the
substitution function :func:`~partialwrap.std_io.sub_params_ja`:

    | &leaflevelctl
    | ...
    | ! slope of stomata model
    | GS_SLOPE = #JA0000#
    | ! intercept of stomata model
    | GS_INTERCEPT = 1.e-3
    | ! minimum stomatal conductance
    | GS_MIN = 1.e-3  ! check influence on very dry sites
    | ! Tuzet parameters
    | GS_HX_HALF = #JA0001#
    | GS_HX_SHAPE = #JA0002#
    | ...
    | /
    
The Python program for optimizing *gs_slope*, *gs_hx_half*, and
*gs_hx_shape* uses the Root Mean Square Error (RMSE) between modelled
and observed NEE:

.. code-block:: python

   # File: optimize_musica.py
   import sys
   from functools import partial
   import numpy as np
   import scipy.optimize as opt
   import netCDF4 as nc
   from partialwrap import exe_wrapper, sub_params_ja

   #
   # Functions
   #

   # Read NEE from MuSICA's netCDF output
   def read_model_nee(ofile):
       with nc.Dataset(ofile, 'r') as fi:
           nee = fi.variables['nee'][:].squeeze()
       return nee

   # The objective: RMSE(obs, mod)
   def rmse(obs, mod):
       return np.sqrt(((obs - mod)**2).mean())

   # RMSE given the output file and the observations
   def calc_rmse(ofile, obs):
       mod = read_model_nee(ofile)
       return rmse(obs, mod)

   # Read NEE observations
   # Assume a csv file such as submitted to europe-fluxdata.org
   #   TIMESTAMP_END,H_1_1_1,LE_1_1_1,FC_1_1_1,FC_SSITC_TEST_1_1_1,TA_1_1_1,...
   #   201601010030,-20.4583,-1.8627,1.9019,0,5.5533,...
   #   ...
   # Assume that timestamps are the same as MuSICA output file
   def read_obs_nee(ofile):
       with open(ofile, 'r') as fi:
           head = fi.readline().strip().split(',')
       iivar = head.index('FC_1_1_1')
       iiflag = head.index('FC_SSITC_TEST_1_1_1')
       dat = np.loadtxt(ofile, delimiter=',', skiprows=1)
       nee = np.ma.array(dat[:, iivar], mask=(dat[:, iiflag] > 0))
       return nee

   # RMSE is around 1-10 (umol m-2 s-1).
   # Use a large random number in case of model error because of
   # odd parameter combinations.
   def err(x):
       return (1. + np.random.random()) * 1000.

   #
   # Setup
   #
   
   # namelist files
   nfiles = ['musica.nml', 'musica_soil.nml', 'fagus_sylvatica.nml']
   # parameter names (not used, only for info)
   names = ['GS_SLOPE', 'GS_HX_HALF', 'GS_HX_SHAPE']
   # lower bounds
   lb = [1.0, -4.5, 2.0]
   # upper bounds
   ub = [13.0, -1.0, 20.]

   # observations
   obsfile = 'FR-Hes_europe-fluxdata_2017.txt'
   obs = read_obs_nee(obsfile)

   # Use a Python wrapper to run musica.exe, which deals with pid
   # and works on all operating systems
   exe             = ['python3', 'run_musica.py']
   parameterfile   = nfiles
   parameterwriter = sub_params_ja
   outputfile   = 'musica_out.nc'
   outputreader = calc_rmse
   oargs        = [obs]  # additional outputreader arguments
   # use exe_wrapper for run_musica.py
   wrap = partial(exe_wrapper, exe,
                  parameterfile, parameterwriter,
                  outputfile, outputreader,
                  {'oargs': oargs,
                   'pid': True, 'error': err})

   #
   # Optimize
   #

   print('Start optimization')
   ncpu = 4
   bounds = list(zip(lb, ub))
   res = opt.differential_evolution(wrap, bounds, workers=ncpu)
   print('Best parameters:', res.x, ' with objective:', res.fun)

   # write parameter files with optimized parameters and suffix .opt
   print('Write parameter files with optimized parameters')
   sub_params_ja('fagus_sylvatica.nml', res.x, pid='opt')

Here we could have also used
:func:`~partialwrap.std_io.sub_params_names` instead of
:func:`~partialwrap.std_io.sub_params_ja` to replace the parameters in
the namelists. This would have changed the setup of *partial*
slightly to:

.. code-block:: python
   
   exe             = ['python3', 'run_musica.py']
   parameterfile   = nfiles
   parameterwriter = sub_params_names
   pargs        = [names]  # additional parameterwriter arguments
   outputfile   = 'musica_out.nc'
   outputreader = calc_rmse
   oargs        = [obs]  # additional outputreader arguments
   # use exe_wrapper for run_musica.py
   wrap = partial(exe_wrapper, exe,
                  parameterfile, parameterwriter,
                  outputfile, outputreader,
                  {'pargs': pargs, 'oargs': oargs,
                   'pid': True, 'error': err})

Only the `parameterwriter` changed and the *names* go into `pargs`,
additional arguments passed to the
`parameterwriter = sub_params_names`. We did not use it here because
`MuSICA`_ can simulate several tree species in the same forest having
one namelist file per species with the same parameters,
e.g. *fagus_sylvatica.nml* and *picea_abeis.nml*.
:func:`~partialwrap.std_io.sub_params_names` would
replace the same value into both species files, while I can prepare
separate tags such as *#JA0000#* for *gs_slope* in
*fagus_sylvatica.nml* and *#JA0003#* for *gs_slope* in
*picea_abies.nml*.

We use almost the same Python wrapper as `run_rastrigin.py` in section
`Using launch scripts`_ to deal with `pid`. We only have to deal with
several parameter files instead of just one file here:

.. code-block:: python

   # File: run_musica.py
   import os
   import shutil
   import subprocess
   import sys

   # get pid
   if len(sys.argv) > 1:
       pid = sys.argv[1]
   else:
       pid = None

   exe    = './musica.exe'
   pfiles = ['musica.nml', 'musica_soil.nml', 'fagus_sylvatica.nml']
   ofile  = 'musica_out.nc'

   # make individual run directory
   if pid is None:
       rundir = 'tmp'
   else:
       rundir = f'tmp.{pid}'
   os.mkdir(rundir)

   # copy individual parameter files
   for pfile in pfiles:
       if pid is None:
           os.rename(f'{pfile}', f'{rundir}/{pfile}')
       else:
           os.rename(f'{pfile}.{pid}', f'{rundir}/{pfile}')

   # run in individual directory
   shutil.copyfile(exe, f'{rundir}/{exe}')
   os.chdir(rundir)
   err = subprocess.check_output([exe], stderr=subprocess.STDOUT)

   # make output available to exe_wrapper
   if pid is None:
       os.rename(ofile, f'../{ofile}')
   else:
       os.rename(ofile, f'../{ofile}.{pid}')

   # clean up
   os.chdir('..')
   shutil.rmtree(rundir)

That's all Folks!

.. _emcee: https://emcee.readthedocs.io/en/latest/
.. _Cython: https://cython.readthedocs.io/
.. _masked arrays: https://numpy.org/doc/stable/reference/maskedarray.html
.. _MuSICA: https://ecofun.ispa.bordeaux.inrae.fr/index.php/musica-model/
