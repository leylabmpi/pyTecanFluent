#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import glob

with open('HISTORY.rst') as history_file:
    history = history_file.read()

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
    version='0.0.1',
    description='Python interface to TECAN Fluent liquid handling robot',
    long_description=long_desc + '\n\n' + history,
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
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
