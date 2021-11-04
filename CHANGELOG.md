# Changelog

All notable changes after its initial development up to June 2020 (v1.0) are documented in this file.

### v1.5 (Nov 2021)
    - Allow empty right-hand sides in files, such as `name1 =`.

### v1.4.x (Oct 2021)
    - Corrected error in docs still using directory structure without src.
    - Updated Zenodo, MANIFEST.in, and documentation.
    - Removed stray print statement.

### v1.4 (Oct 2021)
    - Move to new pip structure using pyproject.toml.
    - Move to Github actions.
    - Remove History sections from docstrings of wrappers and IO routines.

### v1.3.5 (May 2021)
    - Escape backslash also in params in sub_params_ja routines.

### v1.3.4 (May 2021)
    - Escape backslash in params in sub_params_names routines.

### v1.3.3 (May 2021)
    - Deleted trailing print statement. Deleting and replacing v1.3.2 did
      not replace whl on PyPI.

### v1.3.2 (May 2021)
    - Protect saved groups in replacement pattern in sub_params_names routines.

### v1.3.1 (May 2021)
    - Respect space after = sign in sub_params_names routines.
    - Replace everything after space after = sign up to next space.

### v1.3 (May 2021)
    - Correct bug in exe_wrapper if pid=True.
    - Allow non-numeric parameters in std_io functions.

### v1.2 (Dec 2020)
    - Removed scipy.optimize from tests of wrappers.
    - Use build instead of cibuildwheel to build pure Python wheels.

### v1.1 (Dec 2020)
    - Check only Linux on TravisCI.
    - Made all flake8 compatible.

### v1.0 (Jun 2020)
    - Initial release on PyPI.
