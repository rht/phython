from distutils.core import setup
from Cython.Build import cythonize

setup(
        name = "star in galaxy",
        ext_modules = cythonize('physics.py'), # accepts a glob pattern
        )
