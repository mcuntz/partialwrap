partialwrap – Wrappers of external executables and Python functions for functools.partial
=========================================================================================

.. raw:: html

   <!-- pandoc -f gfm -o README.html -t html README.md -->

A small Python library providing wrappers for external executables and
Python functions so that they can easily be partialized with Python’s
functools.partial.

|DOI| |PyPI version| |Conda version| |License| |Build Status| |Coverage Status|


About partialwrap
-----------------

**partialwrap** is a Python library providing easy wrapper functions to use with
Python’s `functools.partial`_. It allows to use any external executable as well
as any Python function with arbitrary arguments and keywords to be used with
libraries that call functions simply in the form ``func(x)``. This includes the
functions of the `scipy.optimize`_ package or external packages such as `emcee`_
or `pyeee`_, and allows the use of distributed function evaluations with
Python’s `multiprocessing`_ or via `MPI`_.


Documentation
-------------

The complete documentation for **partialwrap** is available from Read The
Docs.

   https://mcuntz.github.io/partialwrap/


Quick usage guide
-----------------

Context
~~~~~~~

The easiest wrapper in **partialwrap** is for Python functions.

Consider the `Rastrigin function`_, which is a popular function for performance
testing of optimization algorithms: ``y = a*n + sum_i^n(x_i^2 - a*cos(b*x_i))``

.. code:: python

   import numpy as np
   def rastrigin(x, a, b=2. * np.pi):
       return a * len(x) + np.sum(x**2 - a * np.cos(b * x))

The Rastrigin function has a global minimum of ``0`` at all ``x_i = 0``.
``a`` influences mainly the depth of the (local and global) minima,
whereas ``b`` influences mainly the size of the minima.

A common form uses ``a = 10`` and ``b = 2*pi``. The parameters ``x_i``
should then be in the interval [-5.12, 5.12].

.. code:: python

   import numpy as np
   def rastrigin1(x):
       return 10. * len(x) + np.sum(x**2 - 10. * np.cos(2. * np.pi * x))

One could use the `scipy.optimize`_ package to try to find the minimum of the
Rastrigin function:

.. code:: python

   import scipy.optimize as opt

   x0 = np.array([0.5, 0.5])
   res = opt.minimize(rastrigin1, x0, method='Nelder-Mead')
   print(res.x)
   # returns: [0.99497463 0.99498333]

   res = opt.minimize(rastrigin1, x0, method='BFGS')
   print(res.x)
   # returns: [-7.82960597e-09 -7.82960597e-09]

`scipy.optimize`_ allows passing arguments to the function to minimize.
One could hence use the general Rastrigin function ``rastrigin``
(instead of ``rastrigin1``) to get the same result:

.. code:: python

   res = opt.minimize(rastrigin, x0, args=(10.), method='BFGS')


Simple Python functions
~~~~~~~~~~~~~~~~~~~~~~~

Not all optimizers allow the passing of arguments. And notably `scipy.optimize`_
does not allow the passing of keyword arguments. One can use `partial`_ of
Python’s `functools`_ in this case:

.. code:: python

   from functools import partial

   def call_func_arg_kwarg(func, a, b, x):
      return func(x, a, b=b)

   # Partialize function with fixed parameters
   a = 5.
   b = 4. * np.pi
   partial_rastrigin = partial(call_func_arg_kwarg, rastrigin, a, b)

   res = opt.minimize(partial_rastrigin, x0, method='BFGS')

Figuratively speaking, `partial`_ passes ``a`` and ``b`` to the
function ``call_func_arg_kwarg`` already during definition. ``minimize``
can then simply call it as ``partial_rastrigin(x)``, which finalizes
the call to ``rastrigin(x, a, b=b)``.

**partialwrap** provides a convenience function ``partialwrap.function_wrapper``
passing all arguments, given as a *list*, and keyword arguments, given as a
*dictionary*, to arbitrary functions:

.. code:: python

   from partialwrap import function_wrapper

   args   = [20.]
   kwargs = {'b': 1. * np.pi}
   rastra = partial(function_wrapper, rastrigin, args, kwargs)

   res = opt.minimize(rastra, x0, method='BFGS')

Masking parameters
~~~~~~~~~~~~~~~~~~

A common case in numerical optimization are bound parameters and specifically
the exclusion of some well-known or correlated parameters from optimization.
**partialwrap** provides a convenience function
``partialwrap.function_mask_wrapper`` to include only the masked parameters in
the function evaluation:

.. code:: python

       from partialwrap import function_mask_wrapper

       x0      = np.array([0.5,  0.0001, 0.5])
       mask    = np.array([True, False,  True])
       mrastra = partial(function_mask_wrapper, rastrigin, x0, mask, args, kwargs)

       res        = opt.minimize(mrastra, x0[mask], method='BFGS')
       xout       = x0.copy()
       xout[mask] = res.x

The values of ``x0`` will be taken where ``mask==False``, i.e. ``mask`` could be
called an include-mask.

External executables
~~~~~~~~~~~~~~~~~~~~

**partialwrap** provides wrapper functions to work with external
executables: ``partialwrap.exe_wrapper`` and
``partialwrap.exe_mask_wrapper``.

**partialwrap** writes the sampled parameter sets into files that can be
read by the external program. The program writes its result to a file
that will then be read by **partialwrap** in return. This means
**partialwrap** needs to have a function ``parameterwriter`` that writes
the parameter file ``parameterfile`` needed by the executable ``exe``.
It then needs to have a function ``outputreader`` for reading the output
file ``outputfile`` of the external executable ``exe``.

Consider for simplicity an external Python program
(e.g. ``rastrigin1.py``) that calculates the Rastrigin function with
``a = 10`` and ``b = 2*pi``, reading in an arbitrary number of
parameters ``x_i`` from a ``parameterfile = params.txt`` and writing its
output into an ``outputfile = out.txt``:

.. code:: python

   # File: rastrigin1.py

   # Rastrigin function a=10, b=2*pi
   import numpy as np
   def rastrigin1(x):
       return 10. * len(x) + np.sum(x**2 - 10. * np.cos(2. * np.pi * x))

   # read parameters
   from partialwrap import standard_parameter_reader
   x = standard_parameter_reader('params.txt')

   # calc function
   y = rastrigin1(x)

   # write output file
   with open('out.txt', 'w') as ff:
       print(y, file=ff)

This program can be called on the command line (if `params.txt` is present) with:

.. code:: bash

   python rastrigin1.py

The external program can be used with Python’s ``functools.partial`` and
the wrapper function ``partialwrap.exe_wrapper``:

.. code:: python

   from functools import partial
   from partialwrap import exe_wrapper, standard_parameter_writer, standard_output_reader

   rastrigin_exe  = ['python3', 'rastrigin1.py']
   parameterfile  = 'params.txt'
   outputfile     = 'out.txt'
   rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                            parameterfile, standard_parameter_writer,
                            outputfile, standard_output_reader, {})

   x0  = [0.1, 0.2, 0.3]
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

``partialwrap.standard_parameter_reader`` and
``partialwrap.standard_parameter_writer`` are convenience functions that read
and write one parameter per line in a file without a header. The function
``partialwrap.standard_output_reader`` simply reads one value from a file
without header. The empty dictionary at the end is explained in the
`userguide`_.

One can easily imagine to replace the Python program ``rastrigin1.py`` by any
compiled executable from C, Fortran or alike. See the `userguide`_ for details.


Installation
------------

The easiest way to install is via `pip`:

.. code-block:: bash

   pip install partialwrap

or via `conda`:

.. code-block:: bash

   conda install -c conda-forge partialwrap


Requirements:
-------------

-  `NumPy <https://www.numpy.org>`__


License
-------

**partialwrap** is distributed under the MIT License. See the `LICENSE`_ file
for details.

Copyright (c) 2016-2023 Matthias Cuntz

The project structure is based on a `template`_ provided by `Sebastian Müller`_.


.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3893705.svg
   :target: https://doi.org/10.5281/zenodo.3893705
.. |PyPI version| image:: https://badge.fury.io/py/partialwrap.svg
   :target: https://badge.fury.io/py/partialwrap
.. |Conda version| image:: https://anaconda.org/conda-forge/partialwrap/badges/version.svg
   :target: https://anaconda.org/conda-forge/partialwrap
.. |License| image:: http://img.shields.io/badge/license-MIT-blue.svg?style=flat
   :target: https://github.com/mcuntz/partialwrap/blob/master/LICENSE
.. |Build Status| image:: https://github.com/mcuntz/partialwrap/workflows/Continuous%20Integration/badge.svg?branch=main
   :target: https://github.com/mcuntz/partialwrap/actions
.. |Coverage Status| image:: https://coveralls.io/repos/github/mcuntz/partialwrap/badge.svg?branch=master
   :target: https://coveralls.io/github/mcuntz/partialwrap?branch=master

.. _functools.partial: https://docs.python.org/3/library/functools.html#functools.partial
.. _scipy.optimize: https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html
.. _emcee: https://github.com/dfm/emcee
.. _pyeee: https://github.com/mcuntz/pyeee
.. _multiprocessing: https://docs.python.org/3/library/multiprocessing.html
.. _MPI: https://bitbucket.org/mpi4py/mpi4py
.. _Rastrigin function: https://en.wikipedia.org/wiki/Rastrigin_function
.. _partial: https://docs.python.org/3/library/functools.html#functools.partial
.. _functools: https://docs.python.org/3/library/functools.html
.. _userguide: https://mcuntz.github.io/partialwrap/html/userguide.html
.. _LICENSE: https://github.com/mcuntz/partialwrap/LICENSE
.. _template: https://github.com/MuellerSeb/template
.. _Sebastian Müller: https://github.com/MuellerSeb
