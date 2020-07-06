# partialwrap -- Wrappers of external executables and Python functions for functools.partial
<!-- pandoc -f gfm -o README.html -t html README.md -->

A small Python library providing wrappers for external executables and Python
functions so that they can easily be partialised with Python's
functools.partial.

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3893706.svg)](https://doi.org/10.5281/zenodo.3893706) -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3893705.svg)](https://doi.org/10.5281/zenodo.3893705)
[![PyPI version](https://badge.fury.io/py/partialwrap.svg)](https://badge.fury.io/py/partialwrap)
[![License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/mcuntz/partialwrap/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/mcuntz/partialwrap.svg?branch=master)](https://travis-ci.org/mcuntz/partialwrap)
[![Coverage Status](https://coveralls.io/repos/github/mcuntz/partialwrap/badge.svg?branch=master)](https://coveralls.io/github/mcuntz/partialwrap?branch=master)
[![Documentation Status](https://readthedocs.org/projects/partialwrap/badge/?version=latest)](https://partialwrap.readthedocs.io/en/latest/?badge=latest)

## About partialwrap

*partialwrap* is a Python library providing easy wrapper functions to use with
Python's
[functools.partial](https://docs.python.org/3/library/functools.html#functools.partial).
It allows to use any external executable as well as any Python function with
arbitrary arguments and keywords to be used with libraries that call functions
simply in the form `func(x)`. This includes the functions of the
[scipy.optimize](https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html)
package or external packages such as [emcee](https://github.com/dfm/emcee) or
[pyeee](https://github.com/mcuntz/pyeee), and allows the use of distributed
function evaluations with Python's
[multiprocessing](https://docs.python.org/3/library/multiprocessing.html) or via
[MPI](https://bitbucket.org/mpi4py/mpi4py).


## Documentation

The complete documentation for *partialwrap* is available from Read The Docs.

   http://partialwrap.readthedocs.org/en/latest/


## Quick usage guide


### Context

The easiest wrapper in ``partialwrap`` is for Python functions.

Consider the [Rastrigin
function](https://en.wikipedia.org/wiki/Rastrigin_function), which is a popular
function for performance testing of optimization algorithms: `y = a*n +
sum_i^n(x_i^2 - a*cos(b*x_i))`

```python
    import numpy as np
    def rastrigin(x, a, b=2.*np.pi):
        return a*len(x) + np.sum(x**2 - a*np.cos(b*x))
```

The Rastrigin function has a global minimum of 0 at all `x_i = 0`. `a`
influences mainly the depth of the (local and global) minima, whereas `b`
influences mainly the size of the minima.

A common form uses `a = 10` and `b = 2*pi`. The parameters `x_i` should then be
in the interval [-5.12,5.12].

```python
    import numpy as np
    def rastrigin1(x):
        return 10.*len(x) + np.sum(x**2 - 10.*np.cos(2.*np.pi*x))
```

One could use the
[scipy.optimize](https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html)
package to try to find the minimum of the Rastrigin function:

```python
    import scipy.optimize as opt
    x0 = np.array([0.5, 0.5]) 
    res = opt.minimize(rastrigin1, x0, method='Nelder-Mead')
    print(res.x)
    # retruns: [0.99497463 0.99498333]
    res = opt.minimize(rastrigin1, x0, method='BFGS') 
    print(res.x)
    # returns: [-7.82960597e-09 -7.82960597e-09]
```

`scipy.optimize` allows passing arguments to the function to minimize. One could
hence use the general Rastrigin function `rastrigin` (instead of `rastrigin1`)
to get the same result:

```python
    res = opt.minimize(rastrigin, x0, args=(10.), method='BFGS') 
```


### Simple Python functions

Not all optimizers allow the passing of arguments. And notably `scipy.optimize`
does not allow the passing of keyword arguments. One can use
[partial](https://docs.python.org/3/library/functools.html#functools.partial) of
Python's [functools](https://docs.python.org/3/library/functools.html) in this
case:

```python
    from functools import partial

    def call_func_arg_kwarg(func, a, b, x):
       return func(x, a, b=b)

    # Partialise function with fixed parameters
    a = 5.
    b = 4.*np.pi
    partial_rastrigin = partial(call_func_arg_kwarg, rastrigin, a, b)

    res = opt.minimize(partial_rastrigin, x0, method='BFGS')
```

Figuratively speaking, `partial` passes `a` and `b` to the function `call_func_arg_kwarg`
already during definition. `minimize` can then simply calls it as
`partial_rastrigin(x)`, which finalizes the call to `rastrigin(x, a, b=b)`.

*partialwrap* provides a convenience function `function_wrapper` passing
all arguments, given as a *list*, and keyword arguments, given as a
*dictionary*, to arbitrary functions:

```python
    from partialwrap import function_wrapper

    args   = [20.]
    kwargs = {'b': 1.*np.pi}
    rastra = partial(function_wrapper, rastrigin, args, kwargs)

    res = opt.minimize(rastra, x0, method='BFGS')
```


### Masking parameters

A common case in numerical optimization are bound parameters and specifically
the exclusion of some well-known or correlated parameters from optimization.
*partialwrap* provides a convenience function `function_mask_wrapper` to include
only the masked parameters in the function evaluation:

```python
    from partialwrap import function_mask_wrapper

    x0      = np.array([0.5,  0.0001, 0.5])
    mask    = np.array([True, False,  True])
    mrastra = partial(function_mask_wrapper, rastrigin, x0, mask, args, kwargs)

    res        = opt.minimize(mrastra, x0[mask], method='BFGS')
    xout       = x0.copy()
    xout[mask] = res.x
```

The values of `x0` will be taken where `mask==False`.


External executables
--------------------

*partialwrap* provides wrapper functions to work with external executables:
`partialwrap.exe_wrapper` and `partialwrap.exe_mask_wrapper`.

*partialwrap* writes the sampled parameter sets into files that can be read by
the external program. The program writes its result to a file that will then be
read by *partialwrap* in return. This means *partialwrap* needs to have a
function `parameterwriter` that writes the parameter file `parameterfile` needed
by the executable `exe`. It then needs to have a function `outputreader` for
reading the output file `outputfile` of the external executable `exe`.

Consider for simplicity an external Python program (e.g. `rastrigin1.py`) that
calculates the Rastrigin function with `a = 10` and `b = 2*pi`, reading in an
arbitrary number of parameters `x_i` from a `parameterfile = params.txt` and
writing its output into an `outputfile = out.txt`:

```python
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
```

This program can be called on the command line with:

```bash
    python rastrigin1.py
```

The external program can be used with Python's `functools.partial` and the
wrapper function `partialwrap.exe_wrapper`:

```python
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
```

`partialwrap.standard_parameter_reader` and
`partialwrap.standard_parameter_writer` are convenience functions that read and
write one parameter per line in a file without a header. The function
`partialwrap.standard_output_reader` simply reads one value from a file without
header. The empty dictionary at the end is explained in the
[userguide](https://partialwrap.readthedocs.io/en/latest/userguide.html).

One can easily imagine to replace the Python program `rastrigin1.py` by any
compiled executable from C, Fortran or alike. See the
[userguide](https://partialwrap.readthedocs.io/en/latest/userguide.html) for
details.


## Installation

The easiest way to install is via `pip`:

```bash
    pip install partialwrap
```
	
See the [installation instructions](http://partialwrap.readthedocs.io/en/latest/install.html)
in the [documentation](http://partialwrap.readthedocs.io) for more information.


## Requirements:

- [NumPy](https://www.numpy.org)


## License

*partialwrap* is distributed under the MIT License.  
See the [LICENSE](https://github.com/mcuntz/partialwrap/LICENSE) file for details.

Copyright (c) 2016-2020 Matthias Cuntz

The project structure is based on a
[template](https://github.com/MuellerSeb/template) provided by [Sebastian
Müller](https://github.com/MuellerSeb).

## Contributing to partialwrap

Users are welcome to submit bug reports, feature requests, and code
contributions to this project through GitHub.

More information is available in the
[Contributing](http://partialwrap.readthedocs.org/en/latest/contributing.html)
guidelines.
