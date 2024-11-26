from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import platform

if platform.system() == 'Linux':
    compilation_args = ['-fopenmp']
    link_args = ['-fopenmp']
elif platform.system() == 'Darwin':
    compilation_args = ['-Xpreprocessor', '-fopenmp', '-lomp']
    link_args = ['-lomp']

setup(
    ext_modules=cythonize(["spincalc.pyx"]),
    include_dirs=[numpy.get_include()]
)    
