from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ Extension("compare_sim",
                        ["compare_sim.pyx", "c_sim.c"],
                        extra_compile_args=['-O3', '-fopenmp', '-mmmx', '-msse', '-msse2', '-msse3', '-mssse3',
                            '-msse4', '-msse4.1', '-msse4.2', '-maes'],
                        extra_link_args=['-fopenmp']
                    )
                   ]
)
