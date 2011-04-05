from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("compare_sim", ["compare_sim.pyx"]),
                    Extension("fib", ["fib.c"]),
                    Extension("multiprocess_sim", ["multiprocess_sim.pyx"])]
)
