#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import glob


requirements = [
    'numpy',
    'pandas',
    'xlrd'
]

test_requirements = [
    # TODO: put package test requirements here
]


long_desc = """
Python interface to TECAN Fluent liquid handling robot.
"""


setup(
    name='pyTecanFluent',
    version='0.4.1',
    description='Python interface to TECAN Fluent liquid handling robot',
    long_description=long_desc,
    author="Nick Youngblut",
    author_email='leylabmpi@gmail.com',
    entry_points={
        'console_scripts': [
            'pyTecanFluent = pyTecanFluent.__main__:main'
        ]
    },
    url='https://github.com/leylabmpi/pyTecanFluent',
    packages=find_packages(),
    package_dir={'pyTecanFluent':
                 'pyTecanFluent'},
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='pyTecanFluent',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
    ],
    test_suite='tests',
    tests_require=test_requirements
)
