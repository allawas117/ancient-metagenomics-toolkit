from setuptools import setup, Extension
import numpy as np
import os

def find_htslib():
    for path in ['/usr/local', '/usr']:
        if os.path.exists(os.path.join(path, 'include', 'htslib')):
            return path
    return None

htslib_path = find_htslib()
if not htslib_path:
    raise RuntimeError("HTSlib not found. Please install HTSlib first.")

damage_patterns_module = Extension(
    'amt.cpp_extensions.damage_patterns',
    sources=['amt/cpp_extensions/damage_patterns/damage_patterns.cpp'],
    include_dirs=[
        np.get_include(),
        os.path.join(htslib_path, 'include')
    ],
    library_dirs=[os.path.join(htslib_path, 'lib')],
    libraries=['hts'],
    extra_compile_args=['-std=c++11', '-O3'],
)

setup(
    name='amt_cpp_extensions',
    version='0.1',
    description='C++ extensions for AncientMetagenomicsToolkit',
    ext_modules=[damage_patterns_module],
    python_requires='>=3.6',
)