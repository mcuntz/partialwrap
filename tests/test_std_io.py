#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    This is the unittest for Standard I/O module.

    python -m unittest -v test_std_io.py
    python -m pytest --cov-report term-missing -v test_std_io.py
"""
import unittest


# --------------------------------------------------------------------
# std_io.py
# Missing coverage:
#     518: skip line.startswith('#') in standard_parameter_reader
class TestStd_io(unittest.TestCase):

    def test_std_io_sub_ja(self):
        import os
        import numpy as np
        from partialwrap import sub_params_ja

        filename1 = 'params1.txt'
        filename2 = 'params2.txt'
        pid       = 1234
        params    = np.arange(10, dtype=np.float)

        ff = open(filename1, 'w')
        print('param0 = #JA0000#', file=ff)
        print('param1 = #JA0001#', file=ff)
        print('param2 = #JA0002#', file=ff)
        print('param3 = #JA0003#', file=ff)
        print('param4 = #JA0004#', file=ff)
        ff.close()

        ff = open(filename2, 'w')
        print('param4 = #JA0004#', file=ff)
        print('param5 = #JA0005#', file=ff)
        print('param6 = #JA0006#', file=ff)
        print('param7 = #JA0007#', file=ff)
        ff.close()

        # single file with pid
        sub_params_ja(filename1, params, pid)

        f = open(filename1+'.'+str(pid), 'r')
        lines1 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines1 ],
                         ['param0 = 0.00000000000000e+00',
                          'param1 = 1.00000000000000e+00',
                          'param2 = 2.00000000000000e+00',
                          'param3 = 3.00000000000000e+00',
                          'param4 = 4.00000000000000e+00'])

        # files with pid
        sub_params_ja([filename1, filename2], params, pid)

        f = open(filename1+'.'+str(pid), 'r')
        lines1 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines1 ],
                         ['param0 = 0.00000000000000e+00',
                          'param1 = 1.00000000000000e+00',
                          'param2 = 2.00000000000000e+00',
                          'param3 = 3.00000000000000e+00',
                          'param4 = 4.00000000000000e+00'])

        f = open(filename2+'.'+str(pid), 'r')
        lines2 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines2 ],
                         ['param4 = 4.00000000000000e+00',
                          'param5 = 5.00000000000000e+00',
                          'param6 = 6.00000000000000e+00',
                          'param7 = 7.00000000000000e+00'])

        # no pid
        sub_params_ja(filename2, params)

        f = open(filename2, 'r')
        lines2 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines2 ],
                         ['param4 = 4.00000000000000e+00',
                          'param5 = 5.00000000000000e+00',
                          'param6 = 6.00000000000000e+00',
                          'param7 = 7.00000000000000e+00'])

        # input does not exist
        filename3 = 'params3.txt'
        self.assertRaises(IOError, sub_params_ja, filename3, params)

        # clean up
        if os.path.exists(filename1):
            os.remove(filename1)
        if os.path.exists(filename2):
            os.remove(filename2)
        if os.path.exists(filename1+'.'+str(pid)):
            os.remove(filename1+'.'+str(pid))
        if os.path.exists(filename2+'.'+str(pid)):
            os.remove(filename2+'.'+str(pid))

    def test_std_io_sub_params_names(self):
        import os
        import numpy as np
        from partialwrap import sub_params_names_ignorecase
        from partialwrap import sub_params_names_case, sub_params_names

        # ignore case
        filename1 = 'params11.txt'
        filename2 = 'params21.txt'
        pid       = 1234
        params    = np.arange(10, dtype=np.float)
        names     = ['param0', 'param1', 'param2', 'param3', 'param4',
                     'param5', 'param6', 'param7', 'param8', 'param9']

        ff = open(filename1, 'w')
        print('param0 = -6', file=ff)
        print('Param1 = -7', file=ff)
        print('param2 = -8', file=ff)
        print('Param3 = -9', file=ff)
        print('param4 = -10', file=ff)
        ff.close()

        ff = open(filename2, 'w')
        print('param4 = -10', file=ff)
        print('param5 = 3', file=ff)
        print('PARAM6 = 4', file=ff)
        print('param7 = 5', file=ff)
        ff.close()

        sub_params_names_ignorecase([filename1, filename2], params, names, pid)

        f = open(filename1+'.'+str(pid), 'r')
        lines1 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines1 ],
                         ['param0 = 0.00000000000000e+00',
                          'Param1 = 1.00000000000000e+00',
                          'param2 = 2.00000000000000e+00',
                          'Param3 = 3.00000000000000e+00',
                          'param4 = 4.00000000000000e+00'])

        f = open(filename2+'.'+str(pid), 'r')
        lines2 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines2 ],
                         ['param4 = 4.00000000000000e+00',
                          'param5 = 5.00000000000000e+00',
                          'PARAM6 = 6.00000000000000e+00',
                          'param7 = 7.00000000000000e+00'])

        if os.path.exists(filename1):
            os.remove(filename1)
        if os.path.exists(filename2):
            os.remove(filename2)
        if os.path.exists(filename1+'.'+str(pid)):
            os.remove(filename1+'.'+str(pid))
        if os.path.exists(filename2+'.'+str(pid)):
            os.remove(filename2+'.'+str(pid))

        # wrapper for ignore case
        filename1 = 'params12.txt'
        filename2 = 'params22.txt'
        pid       = 1234
        params    = np.arange(10, dtype=np.float)
        names     = ['param0', 'param1', 'param2', 'param3', 'param4',
                     'param5', 'param6', 'param7', 'param8', 'param9']

        ff = open(filename1, 'w')
        print('param0 = -6', file=ff)
        print('Param1 = -7', file=ff)
        print('param2 = -8', file=ff)
        print('Param3 = -9', file=ff)
        print('param4 = -10', file=ff)
        ff.close()

        ff = open(filename2, 'w')
        print('param4 = -10', file=ff)
        print('param5 = 3', file=ff)
        print('PARAM6 = 4', file=ff)
        print('param7 = 5', file=ff)
        ff.close()

        sub_params_names([filename1, filename2], params, names, pid)

        f = open(filename1+'.'+str(pid), 'r')
        lines1 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines1 ],
                         ['param0 = 0.00000000000000e+00',
                          'Param1 = 1.00000000000000e+00',
                          'param2 = 2.00000000000000e+00',
                          'Param3 = 3.00000000000000e+00',
                          'param4 = 4.00000000000000e+00'])

        f = open(filename2+'.'+str(pid), 'r')
        lines2 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines2 ],
                         ['param4 = 4.00000000000000e+00',
                          'param5 = 5.00000000000000e+00',
                          'PARAM6 = 6.00000000000000e+00',
                          'param7 = 7.00000000000000e+00'])

        if os.path.exists(filename1):
            os.remove(filename1)
        if os.path.exists(filename2):
            os.remove(filename2)
        if os.path.exists(filename1+'.'+str(pid)):
            os.remove(filename1+'.'+str(pid))
        if os.path.exists(filename2+'.'+str(pid)):
            os.remove(filename2+'.'+str(pid))

        # case sensitive
        filename1 = 'params13.txt'
        filename2 = 'params23.txt'
        pid       = 1234
        params    = np.arange(10, dtype=np.float)
        names     = ['param0', 'param1', 'param2', 'param3', 'param4',
                     'param5', 'param6', 'param7', 'param8', 'param9']

        ff = open(filename1, 'w')
        print('    param0= -6', file=ff)
        print('   Param1 = -7', file=ff)
        print('  param2  = -8', file=ff)
        print(' Param3   = -9', file=ff)
        print('param4    = -10', file=ff)
        ff.close()

        ff = open(filename2, 'w')
        print('param4    = -10', file=ff)
        print('param5   = 3', file=ff)
        print('PARAM6  = 4', file=ff)
        print('param7 = 5', file=ff)
        ff.close()

        sub_params_names_case(filename1, params, names, pid)
        sub_params_names_case(filename2, params, names, pid)

        f = open(filename1+'.'+str(pid), 'r')
        lines1 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines1 ],
                         ['    param0= 0.00000000000000e+00',
                          '   Param1 = -7',
                          '  param2  = 2.00000000000000e+00',
                          ' Param3   = -9',
                          'param4    = 4.00000000000000e+00'])

        f = open(filename2+'.'+str(pid), 'r')
        lines2 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines2 ],
                         ['param4    = 4.00000000000000e+00',
                          'param5   = 5.00000000000000e+00',
                          'PARAM6  = 4',
                          'param7 = 7.00000000000000e+00'])

        # no pid
        sub_params_names_case([filename1, filename2], params, names)

        f = open(filename1, 'r')
        lines1 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines1 ],
                         ['    param0= 0.00000000000000e+00',
                          '   Param1 = -7',
                          '  param2  = 2.00000000000000e+00',
                          ' Param3   = -9',
                          'param4    = 4.00000000000000e+00'])

        f = open(filename2, 'r')
        lines2 = f.readlines()
        f.close()

        self.assertEqual([ i.rstrip() for i in lines2 ],
                         ['param4    = 4.00000000000000e+00',
                          'param5   = 5.00000000000000e+00',
                          'PARAM6  = 4',
                          'param7 = 7.00000000000000e+00'])

        if os.path.exists(filename1):
            os.remove(filename1)
        if os.path.exists(filename2):
            os.remove(filename2)
        if os.path.exists(filename1+'.'+str(pid)):
            os.remove(filename1+'.'+str(pid))
        if os.path.exists(filename2+'.'+str(pid)):
            os.remove(filename2+'.'+str(pid))

    def test_std_io_standard(self):
        import os
        import numpy as np
        from partialwrap import standard_parameter_writer
        from partialwrap import standard_parameter_reader
        from partialwrap import standard_parameter_writer_bounds_mask
        from partialwrap import standard_parameter_reader_bounds_mask
        from partialwrap import standard_output_reader
        from partialwrap import standard_timeseries_reader
        from partialwrap import standard_time_series_reader

        # standard_parameter_reader/writer without pid
        filename = 'params.txt'
        params   = np.arange(10, dtype=np.float)
        standard_parameter_writer(filename, params)

        iparams = standard_parameter_reader(filename)

        self.assertEqual(list(iparams), list(params))

        if os.path.exists(filename):
            os.remove(filename)

        # standard_parameter_writer with pid
        filename = 'params.txt'
        pid      = 1234
        params   = np.arange(10, dtype=np.float)
        standard_parameter_writer(filename, params, pid)

        iparams = standard_parameter_reader(filename, pid)

        self.assertEqual(list(iparams), list(params))

        if os.path.exists(filename+'.'+str(pid)):
            os.remove(filename+'.'+str(pid))

        # standard_parameter_reader/writer_bounds_mask
        filename = 'params.txt'
        pid      = 1234
        params   = np.arange(10, dtype=np.float)
        pmin     = params - 1.
        pmax     = params + 1.
        mask     = np.ones(10, dtype=np.bool)
        standard_parameter_writer_bounds_mask(filename, params, pmin, pmax,
                                              mask, pid)

        ids, iparams, ipmin, ipmax, imask = (
            standard_parameter_reader_bounds_mask(filename, pid))

        self.assertEqual(list(ids),     list([ str(i)
                                               for i in np.arange(10)+1 ]))
        self.assertEqual(list(iparams), list(params))
        self.assertEqual(list(ipmin),   list(pmin))
        self.assertEqual(list(ipmax),   list(pmax))
        self.assertEqual(list(imask),   list(mask))

        if os.path.exists(filename+'.'+str(pid)):
            os.remove(filename+'.'+str(pid))

        # standard_parameter_reader_bounds_mask - IOError
        filename = 'params.txt'
        params   = np.arange(10, dtype=np.float)
        pmin     = params - 1.
        pmax     = params + 1.
        mask     = np.ones(10, dtype=np.bool)
        ff = open(filename, 'w')
        for i in range(10):
            dstr = '{:d} {:.14e} {:.14e} {:.14e}'.format(
                i+1, params[i], pmin[i], pmax[i])
            print(dstr, file=ff)
        ff.close()

        self.assertRaises(IOError, standard_parameter_reader_bounds_mask,
                          filename)

        if os.path.exists(filename):
            os.remove(filename)

        # standard_parameter_writer_bounds_mask - no pid
        filename = 'params.txt'
        params   = np.arange(10, dtype=np.float)
        pmin     = params - 1.
        pmax     = params + 1.
        mask     = np.ones(10, dtype=np.bool)
        standard_parameter_writer_bounds_mask(
            filename, params, pmin, pmax, mask)

        ids, iparams, ipmin, ipmax, imask = (
            standard_parameter_reader_bounds_mask(filename))

        self.assertEqual(list(ids),     list([ str(i)
                                               for i in np.arange(10)+1 ]))
        self.assertEqual(list(iparams), list(params))
        self.assertEqual(list(ipmin),   list(pmin))
        self.assertEqual(list(ipmax),   list(pmax))
        self.assertEqual(list(imask),   list(mask))

        if os.path.exists(filename):
            os.remove(filename)

        # standard_parameter_writer_bounds_mask - pid=None
        filename = 'params.txt'
        pid      = None
        params   = np.arange(10, dtype=np.float)
        pmin     = params - 1.
        pmax     = params + 1.
        mask     = np.ones(10, dtype=np.bool)
        standard_parameter_writer_bounds_mask(filename, params, pmin, pmax,
                                              mask, pid)

        ids, iparams, ipmin, ipmax, imask = (
            standard_parameter_reader_bounds_mask(filename, pid))

        self.assertEqual(list(ids),     list([ str(i)
                                               for i in np.arange(10)+1 ]))
        self.assertEqual(list(iparams), list(params))
        self.assertEqual(list(ipmin),   list(pmin))
        self.assertEqual(list(ipmax),   list(pmax))
        self.assertEqual(list(imask),   list(mask))

        if os.path.exists(filename):
            os.remove(filename)

        # standard_output_reader
        filename = 'obj.txt'

        ff = open(filename, 'w')
        print('{:.14e}'.format(1234.), file=ff)
        ff.close()

        obj = standard_output_reader(filename)
        self.assertEqual(obj, 1234.)

        if os.path.exists(filename):
            os.remove(filename)

        # standard_output_reader, pid
        pid = 1234
        filename = 'obj.txt'

        ff = open(filename+'.'+str(pid), 'w')
        print('{:.14e}'.format(1234.), file=ff)
        ff.close()

        obj = standard_output_reader(filename, pid)
        self.assertEqual(obj, 1234.)

        if os.path.exists(filename+'.'+str(pid)):
            os.remove(filename+'.'+str(pid))

        # standard_time_series_reader
        filename = 'ts.txt'
        params   = np.arange(10, dtype=np.float)

        ff = open(filename, 'w')
        for i in params:
            print('{:.14e}'.format(i), file=ff)
        ff.close()

        ts = standard_time_series_reader(filename)
        self.assertEqual(list(ts), list(params))

        if os.path.exists(filename):
            os.remove(filename)

        # standard_timeseries_reader, pid
        pid = 1234
        filename = 'ts.txt'
        params   = np.arange(10, dtype=np.float)

        ff = open(filename+'.'+str(pid), 'w')
        for i in params:
            print('{:.14e}'.format(i), file=ff)
        ff.close()

        ts = standard_timeseries_reader(filename, pid)
        self.assertEqual(list(ts), list(params))

        if os.path.exists(filename+'.'+str(pid)):
            os.remove(filename+'.'+str(pid))


if __name__ == "__main__":
    unittest.main()
