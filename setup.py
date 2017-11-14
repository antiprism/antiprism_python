"""A setuptools based setup module.
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import listdir
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

prog_dir = 'anti_lib_progs'
src_path = path.join(here, prog_dir)
script_progs = []
for f in listdir(src_path):
    if path.isfile(path.join(src_path, f)) and f not in [
            '__init__.py', 'anti_lib.py']:
        script_progs.append(f+'='+prog_dir+'.'+f[:-3]+':main')

setup(
    name='antiprism_python',
    version='0.1.2',
    description='Scripts to generate various types of polyhedra',
    long_description=long_description,
    url='https://github.com/antiprism/antiprism_python',

    author='Adrian Rossiter',
    author_email='adrian@antiprism.com',

    license='MIT',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Visualization',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3 :: Only',
    ],

    keywords='antiprism polyhdron polyhedra',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=[]),
    py_modules=['anti_lib'],

    entry_points={
        'console_scripts': script_progs,
    },
)
