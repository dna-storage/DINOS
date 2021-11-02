# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.core import setup, Extension

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='overhang',
    version='1.0.0',
    description='DINOS experiments code',
    long_description=readme,
    author='Kevin Volkel',
    author_email='kvolkel@ncsu.edu',
    url='https://github.ncsu.edu/kvolkel/',
    license=license,
    packages=find_packages(exclude=('tests','docs', 'tools'))
)
