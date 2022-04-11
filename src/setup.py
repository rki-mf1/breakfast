#from setuptools import setup, find_packages
#import os
#from distutils.core import setup
#from distutils.extension import Extension

#from Cython.Build import cythonize

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy


#os.environ["CC"] = "g++-5"
'''
setup(
    name='_pairwise_fast',
    ext_modules=cythonize("_pairwise_fast.pyx"),
    zip_safe=False,
    #libraries=["m"],
    cython_compile_time_env=dict(OPENMP=True),
    include_dirs=[numpy.get_include()],
    extra_compile_args=["-O3", "-ffast-math", "-march=native","-fopenmp"],
    extra_link_args=['-fopenmp'],
)
'''

ext_module = Extension(
    "_pairwise_fast",
    ["_pairwise_fast.pyx"],
    extra_compile_args=['-fopenmp'],
    extra_link_args=['-fopenmp'],
)

setup(
    name = "_pairwise_fast",
    cmdclass = {'build_ext': build_ext},
    include_dirs=[numpy.get_include()],
    ext_modules = [ext_module],
)



