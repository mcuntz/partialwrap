Quickstart
==========

``partialwrap``: a small Python library providing wrappers for
external executables to be used easily with Python’s
:py:func:`functools.partial`.

|DOI| |PyPI version| |Conda version| |License| |Build Status| |Coverage Status|


About partialwrap
-----------------

``partialwrap`` is a Python library providing easy wrapper functions
to use with Python’s :py:func:`functools.partial`. Partial's ingenious
mechanism allows to use even very complex functions with many
arguments and keyword arguments with routines that need functions in
the simple form `func(x)`. This includes Python's :func:`map`
function, the minimizers in :mod:`scipy.optimize`, and plenty of
third-party modules such as `emcee`_ or :mod:`pyeee`. ``partialwrap``
allows to use any external executable as well as any Python function
with arbitrary arguments and keywords to be used with
:py:func:`functools.partial`. It also allows to use the executables or
functions with easy parallelization of code with, for example, the
parallel `map` function of the :mod:`multiprocessing` module or via
the Message Passing Interface (MPI) :mod:`mpi4py`.


Quick usage guide
-----------------

``partialwrap`` provides two wrapper functions to work with external
executables: :func:`~partialwrap.wrappers.exe_wrapper` and
:func:`~partialwrap.wrappers.exe_mask_wrapper`.

``partialwrap`` writes the input arguments to the wrapper functions
into files that can be read by the external executable. The executable
should write its result to a file that will then be read by
``partialwrap`` in return. This means ``partialwrap`` needs to have a
function `parameterwriter` that writes the parameter file
`parameterfile` in the form needed by the executable `exe`. It then
also needs to have a function `outputreader` for reading the output
file `outputfile` of the external executable, perhaps calculating an
objective function value.

Take the *Rastrigin function*
(https://en.wikipedia.org/wiki/Rastrigin_function), which is a popular
function for performance testing of optimization algorithms: :math:`y
= a \cdot n + \sum_i^n (x_i^2 - a \cos(b \cdot x_i))`. It has a global
minimum of :math:`0` at all :math:`x_i = 0`. :math:`a` influences
mainly the depth of the (local and global) minima, whereas :math:`b`
influences mainly the size of the minima. A common form uses :math:`a
= 10` and :math:`b = 2 \pi`. The parameters :math:`x_i` should then be
in the interval :math:`[-5.12, 5.12]`.

Consider for simplicity an external Python program
(e.g. `rastrigin1.py`) that calculates the Rastrigin function with
:math:`a = 10` and :math:`b = 2 \pi`, reading in an arbitrary number
of parameters :math:`x_i` from a `parameterfile = params.txt` and
writing its output into an `outputfile = out.txt`:

.. code-block:: python

   # File: rastrigin1.py
   import numpy as np
   from partialwrap import standard_parameter_reader

   # Rastrigin function a=10, b=2*pi
   def rastrigin1(x):
       return 10. * len(x) + np.sum(x**2 - 10. * np.cos(2. * np.pi * x))

   # read parameters
   x = standard_parameter_reader('params.txt')

   # calc function
   y = rastrigin1(x)

   # write output file
   with open('out.txt', 'w') as ff:
       print(y, file=ff)

This program can be called on the command line (if `params.txt` is
present) with:

.. code-block:: bash

   python3 rastrigin1.py

The external program calculating the Rastrigin function could, of
course, also be written in any compiled language such as C or
Fortran. See the `User Guide <userguide.html>`_ for details. The
external program, here the Python version, can be used with Python's
:py:func:`functools.partial` and the wrapper function
:func:`~partialwrap.wrappers.exe_wrapper`:

.. code-block:: python

   import scipy.optimize as opt
   from functools import partial
   from partialwrap import exe_wrapper
   from partialwrap import standard_parameter_writer, standard_output_reader

   rastrigin_exe   = ['python3', 'rastrigin1.py']
   parameterfile   = 'params.txt'
   parameterwriter = standard_parameter_writer
   outputfile      = 'out.txt'
   outputreader    = standard_output_reader
   rastrigin_wrap  = partial(exe_wrapper, rastrigin_exe,
                             parameterfile, parameterwriter,
                             outputfile, outputreader, {})

   x0  = [0.1, 0.2, 0.3]
   res = opt.minimize(rastrigin_wrap, x0, method='BFGS')

The :mod:`scipy.optimize` function :func:`~scipy.optimize.minimize`
passes its sampled parameters to `exe_wrapper`, which writes it to the
file `parameterfile = 'params.txt'`. It then calls `rastrigin_exe =
'python3 rastrigin1.py'` and reads its `outputfile = 'out.txt'`.
:func:`~partialwrap.std_io.standard_parameter_reader` and
:func:`~partialwrap.std_io.standard_parameter_writer` are convenience
functions that read and write one parameter per line in a file without
a header. The empty dictionary at the end is explained in the
`User Guide <userguide.html>`_.

More elaborate input/output of the external program can simply be done
by replacing :func:`~partialwrap.std_io.standard_parameter_reader` and
:func:`~partialwrap.std_io.standard_parameter_writer` with appropriate
functions, while the rest stays pretty much the same.


Installation
------------

The easiest way to install is via `pip`:

.. code-block:: bash

   pip install partialwrap

or via `conda`:

.. code-block:: bash

   conda install -c conda-forge partialwrap


Requirements
------------

- :mod:`numpy`


License
-------

``partialwrap`` is distributed under the MIT License. See the
`LICENSE`_ file for details.

Copyright (c) 2016-2023 Matthias Cuntz

The project structure is based on a `template`_ provided by `Sebastian Müller`_.


Index and Tables
----------------

* :ref:`genindex`
* :ref:`modindex`


.. |DOI|
   image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3893705.svg
   :target: https://doi.org/10.5281/zenodo.3893705
.. |PyPI version|
   image:: https://badge.fury.io/py/partialwrap.svg
   :target: https://badge.fury.io/py/partialwrap
.. |Conda version|
   image:: https://anaconda.org/conda-forge/partialwrap/badges/version.svg
   :target: https://anaconda.org/conda-forge/partialwrap
.. |License|
   image:: http://img.shields.io/badge/license-MIT-blue.svg?style=flat
   :target: https://github.com/mcuntz/partialwrap/blob/master/LICENSE
.. |Build Status|
   image:: https://github.com/mcuntz/partialwrap/workflows/Continuous%20Integration/badge.svg?branch=main
   :target: https://github.com/mcuntz/partialwrap/actions
.. |Coverage Status|
   image:: https://coveralls.io/repos/github/mcuntz/partialwrap/badge.svg?branch=master
   :target: https://coveralls.io/github/mcuntz/partialwrap?branch=master

.. _emcee: https://emcee.readthedocs.io/en/latest/
.. _MPI: https://bitbucket.org/mpi4py/mpi4py
.. _LICENSE: https://github.com/mcuntz/partialwrap/LICENSE
.. _template: https://github.com/MuellerSeb/template
.. _Sebastian Müller: https://github.com/MuellerSeb
