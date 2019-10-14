#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
        'numpy >= 1.13.3',
        'scipy >= 1.1.0',
        'matplotlib >= 2.2.2',
        'shapely >= 1.6.4']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Exneyder A. Montoya-Araque & Ludger O. Suarez-Burgoa",
    author_email='eamontoyaa@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    description="Application software to evaluate the stability of slopes made of Blocks-In-Matrix materials",
    entry_points={
        'console_scripts': [
            'pybimstab=pybimstab.cli:main',
        ],
    },
    install_requires=requirements,
    license="BSD 2-Clause License",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords=['slope stability', 'GLE', 'bimsoil', 'bimrock', 'tortuosity', 'A-star', 'Python', 'application software'],
    name='pybimstab',
    packages=find_packages(exclude=['docs', 'tests']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/eamontoyaa/pybimstab',
    version='0.1.4',
    zip_safe=False,
)
