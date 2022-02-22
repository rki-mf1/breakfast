from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    name='_pairwise_fast',
    ext_modules=cythonize("_pairwise_fast.pyx"),
    zip_safe=False,
    include_dirs=[numpy.get_include()]
)
