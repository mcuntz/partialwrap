#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    This is the unittest for the wrappers module.

    python -m unittest -v tests/test_wrappers.py
    python -m pytest --cov-report term-missing -v tests/test_wrappers.py
"""
import unittest
import numpy as np


# --------------------------------------------------------------------
# Rastrigin function
def rastrigin(x, a, b=2.*np.pi):
    return a*x.size + np.sum(x**2 - a*np.cos(b*x))


# Rastrigin function with a=10 and b=2*pi
def rastrigin1(x):
    return 10.*x.size + np.sum(x**2 - 10.*np.cos(2.*np.pi*x))


# --------------------------------------------------------------------
# wrappers.py
# Missing coverage:
#     114-117: TypeError in try: of internal _tolist helper function,
#              to catch input=None.
class TestWrappers(unittest.TestCase):

    def setUp(self):
        import numpy as np
        # seed for reproducible results
        seed = 1234
        np.random.seed(seed=seed)

    # function_wrapper
    def test_function_wrapper(self):
        from functools import partial
        import numpy as np
        from partialwrap import function_wrapper

        ndim   = 5
        xmin   = -5.12
        xmax   = 5.12
        args   = [20.]
        kwargs = {'b': 1.*np.pi}
        rastra = partial(function_wrapper, rastrigin, args, kwargs)

        res = rastra(np.ones(ndim)*(xmin+0.5*(xmax-xmin)))
        self.assertEqual('{:.4g}'.format(res), '0')

    # function_mask_wrapper
    def test_function_mask_wrapper(self):
        from functools import partial
        import numpy as np
        from partialwrap import function_mask_wrapper

        xmin   = -5.12
        xmax   = 5.12
        x0     = np.array([0.5, 0.0001, 0.5])
        # Do not optimize the second parameter but take its initial
        # value 0.0001
        mask   = [True, False, True]
        args   = [10.]
        kwargs = {'b': 2.*np.pi}
        rastra = partial(function_mask_wrapper, rastrigin, x0, mask, args,
                         kwargs)

        res = rastra(np.ones(np.sum(mask))*(xmin+0.5*(xmax-xmin)))
        self.assertEqual('{:.4g}'.format(res), '1.984e-06')

    # exe_wrapper, w/o kwarg
    def test_exe_wrapper(self):
        from functools import partial
        from partialwrap import exe_wrapper, standard_parameter_writer
        from partialwrap import standard_output_reader

        rastrigin_exe  = ['python3', 'tests/rastrigin1.py']
        ndim           = 2
        xmin           = -5.12
        xmax           = 5.12
        parameterfile  = 'params.txt'
        outputfile     = 'out.txt'
        rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                                 parameterfile, standard_parameter_writer,
                                 outputfile, standard_output_reader, {})

        x0  = [0.1, 0.2]
        res = rastrigin_wrap(x0)
        self.assertEqual('{:.4g}'.format(res), '8.87')

        rastrigin_wrap = partial(exe_wrapper, rastrigin1,
                                 parameterfile, standard_parameter_writer,
                                 outputfile, standard_output_reader, {})
        self.assertRaises(TypeError, rastrigin_wrap, x0)

    # exe_wrapper, w/o kwarg, pid
    def test_exe_wrapper_pid(self):
        from functools import partial
        from partialwrap import exe_wrapper, standard_parameter_writer
        from partialwrap import standard_output_reader

        rastrigin_exe  = ['python3', 'tests/rastrigin1.py']
        ndim           = 2
        xmin           = -5.12
        xmax           = 5.12
        parameterfile  = 'params.txt'
        outputfile     = 'out.txt'
        rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                                 parameterfile, standard_parameter_writer,
                                 outputfile, standard_output_reader,
                                 {'pid': True})

        x0  = [0.1, 0.2]
        res = rastrigin_wrap(x0)
        self.assertEqual('{:.4g}'.format(res), '8.87')

    # exe_wrapper, shell, debug
    def test_exe_wrapper_shell_debug(self):
        from functools import partial
        from partialwrap import exe_wrapper, standard_parameter_writer
        from partialwrap import standard_output_reader

        rastrigin_exe  = 'python3 tests/rastrigin1.py'
        ndim           = 2
        xmin           = -5.12
        xmax           = 5.12
        parameterfile  = 'params.txt'
        outputfile     = 'out.txt'
        rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                                 parameterfile, standard_parameter_writer,
                                 outputfile, standard_output_reader,
                                 {'shell': True, 'debug': True})

        x0  = [0.1, 0.2]
        res = rastrigin_wrap(x0)
        self.assertEqual('{:.4g}'.format(res), '8.87')

    # exe_wrapper, shell, debug, pid
    def test_exe_wrapper_shell_debug_pid(self):
        from functools import partial
        from partialwrap import exe_wrapper, standard_parameter_writer
        from partialwrap import standard_output_reader

        rastrigin_exe  = 'python3 tests/rastrigin1.py'
        ndim           = 2
        xmin           = -5.12
        xmax           = 5.12
        parameterfile  = 'params.txt'
        outputfile     = 'out.txt'
        rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                                 parameterfile, standard_parameter_writer,
                                 outputfile, standard_output_reader,
                                 {'shell': True, 'debug': True, 'pid': True})

        x0  = [0.1, 0.2]
        res = rastrigin_wrap(x0)
        self.assertEqual('{:.4g}'.format(res), '8.87')

    # exe_wrapper, pargs, keepparameterfile
    def test_exe_wrapper_pargs_keepparameterfile(self):
        import os
        from functools import partial
        from partialwrap import exe_wrapper, sub_params_names
        from partialwrap import standard_output_reader

        parameterfile  = 'params.txt'
        ff = open(parameterfile, 'w')
        print('# File: params.txt', file=ff)
        print(' param02 = 0.2 // C comment', file=ff)
        print('param01 = 0.1   ! Fortran comment', file=ff)
        ff.close()

        rastrigin_exe  = ['python3', 'tests/rastrigin2.py']
        ndim           = 2
        xmin           = -5.12
        xmax           = 5.12
        outputfile     = 'out.txt'
        x0             = [ 0.1,       0.2]
        names          = ['param01', 'param02']
        rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                                 parameterfile, sub_params_names,
                                 outputfile, standard_output_reader,
                                 {'pargs': [names], 'keepparameterfile': True})

        res = rastrigin_wrap(x0)
        self.assertEqual('{:.4g}'.format(res), '8.87')
        self.assertTrue(os.path.exists(parameterfile))
        self.assertFalse(os.path.exists(outputfile))

        # clean up
        if os.path.exists(parameterfile):
            os.remove(parameterfile)

    # exe_wrapper, pargs, keepparameterfile, keepoutputfile
    def test_exe_wrapper_pargs_keepparameterfile_keepoutputfile(self):
        import os
        from functools import partial
        from partialwrap import exe_wrapper, sub_params_names
        from partialwrap import standard_output_reader

        parameterfile  = 'params.txt'
        ff = open(parameterfile, 'w')
        print('# File: params.txt', file=ff)
        print(' param02 = 0.2 // C comment', file=ff)
        print('param01 = 0.1   ! Fortran comment', file=ff)
        ff.close()

        rastrigin_exe  = ['python3', 'tests/rastrigin2.py']
        ndim           = 2
        xmin           = -5.12
        xmax           = 5.12
        outputfile     = 'out.txt'
        x0             = [ 0.1,       0.2]
        names          = ['param01', 'param02']
        rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                                 parameterfile, sub_params_names,
                                 outputfile, standard_output_reader,
                                 {'pargs': [names], 'keepparameterfile': True,
                                  'keepoutputfile': True})

        res = rastrigin_wrap(x0)
        self.assertEqual('{:.4g}'.format(res), '8.87')
        self.assertTrue(os.path.exists(parameterfile))
        self.assertTrue(os.path.exists(outputfile))

        # clean up
        if os.path.exists(parameterfile):
            os.remove(parameterfile)
        if os.path.exists(outputfile):
            os.remove(outputfile)

    # exe_wrapper, parameterfiles
    def test_exe_wrapper_parameterfiles(self):
        import os
        from functools import partial
        from partialwrap import exe_wrapper, sub_params_names
        from partialwrap import standard_output_reader

        parameterfile1  = 'params1.txt'
        ff = open(parameterfile1, 'w')
        print('# File: params1.txt', file=ff)
        print('param01 = 0.1   ! Fortran comment', file=ff)
        ff.close()
        parameterfile2  = 'params2.txt'
        ff = open(parameterfile2, 'w')
        print('# File: params2.txt', file=ff)
        print(' param02 = 0.2 // C comment', file=ff)
        ff.close()

        rastrigin_exe  = ['python3', 'tests/rastrigin3.py']
        ndim           = 2
        xmin           = -5.12
        xmax           = 5.12
        parameterfile  = (parameterfile1, parameterfile2)
        outputfile     = 'out.txt'
        x0             = [ 0.1,       0.2]
        names          = ['param01', 'param02']
        rastrigin_wrap = partial(exe_wrapper, rastrigin_exe,
                                 parameterfile, sub_params_names,
                                 outputfile, standard_output_reader,
                                 {'pargs': [names], 'keepparameterfile': True})

        res = rastrigin_wrap(x0)
        self.assertEqual('{:.4g}'.format(res), '8.87')
        self.assertTrue(os.path.exists(parameterfile1))
        self.assertTrue(os.path.exists(parameterfile2))
        self.assertFalse(os.path.exists(outputfile))

        # clean up
        if os.path.exists(parameterfile1):
            os.remove(parameterfile1)
        if os.path.exists(parameterfile2):
            os.remove(parameterfile2)

    # exe_mask_wrapper, w/o kwarg
    def test_exe_mask_wrapper(self):
        from functools import partial
        import numpy as np
        from partialwrap import exe_mask_wrapper, standard_parameter_writer
        from partialwrap import standard_output_reader

        rastrigin_exe  = ['python3', 'tests/rastrigin1.py']
        xmin           = -5.12
        xmax           = 5.12
        x0             = np.array([0.1, 0.0001, 0.3])
        # Do not optimize the second parameter but take its initial
        # value 0.0001
        mask           = [True, False, True]
        parameterfile  = 'params.txt'
        outputfile     = 'out.txt'
        rastrigin_wrap = partial(exe_mask_wrapper, rastrigin_exe, x0, mask,
                                 parameterfile, standard_parameter_writer,
                                 outputfile, standard_output_reader, {})

        res = rastrigin_wrap(x0[mask])
        self.assertEqual('{:.4g}'.format(res), '15.1')

    # exe_mask_wrapper, pargs, keepparameterfile, keepoutputfile
    def test_exe_mask_wrapper_pargs_keepparameterfile_keepoutputfile(self):
        import os
        from functools import partial
        import numpy as np
        from partialwrap import exe_mask_wrapper, sub_params_names
        from partialwrap import standard_output_reader

        parameterfile  = 'params.txt'
        ff = open(parameterfile, 'w')
        print('# File: params.txt', file=ff)
        print(' param02 = 0.0001 // C comment', file=ff)
        print('param01 = 0.1      ! Fortran comment', file=ff)
        print('  param03 = 0.3    # Python comment', file=ff)
        ff.close()

        rastrigin_exe  = ['python3', 'tests/rastrigin2.py']
        ndim           = 3
        xmin           = -5.12
        xmax           = 5.12
        outputfile     = 'out.txt'
        x0             = np.array([ 0.1,       0.0001,    0.3])
        names          = ['param01', 'param02', 'param03']
        # Do not optimize the second parameter but take its initial
        # value 0.0001
        mask           = [ True, False, True]
        rastrigin_wrap = partial(exe_mask_wrapper, rastrigin_exe, x0, mask,
                                 parameterfile, sub_params_names,
                                 outputfile, standard_output_reader,
                                 {'pargs': [names], 'keepparameterfile': True,
                                  'keepoutputfile': True})

        res = rastrigin_wrap(x0[mask])
        self.assertEqual('{:.4g}'.format(res), '15.1')
        self.assertTrue(os.path.exists(parameterfile))
        self.assertTrue(os.path.exists(outputfile))

        # clean up
        if os.path.exists(parameterfile):
            os.remove(parameterfile)
        if os.path.exists(outputfile):
            os.remove(outputfile)


if __name__ == "__main__":
    unittest.main()
