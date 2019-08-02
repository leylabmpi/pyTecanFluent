[![Build Status](https://travis-ci.org/leylabmpi/pyTecanFluent.svg?branch=master)](https://travis-ci.org/leylabmpi/pyTecanFluent)
[![PyPI version](https://badge.fury.io/py/pyTecanFluent.svg)](http://badge.fury.io/py/pyTecanFluent)


pyTecanFluent
=============

Python interface to TECAN Fluent liquid handling robot

> WARNING: this package is under heavy development and is subject to change at any time.


#### Sections

- [Contents](#contents)
- [Examples](#examples)
- [Installation](#installation)
- [Usage](#usage)
- [Changelog](#changelog)
- [License](#license)


## Contents

[[top](#sections)]



## Examples

[[top](#sections)]

* See the bash scripts in:
  * `./examples/`


## Installation

[[top](#sections)]

### Altering the FluentControl database

The database includes:

* labware
* liquid classes
* tip types
* target positions

To edit the database, just alter the JSON files in `./pyTecanFluent/database/` prior to installing.


### With pip

`pip install pyTecanFluent`

### From source

Optional testing:

`python setup.py test`

Install:

`python setup.py install`


## Changelog

[[top](#sections)]


# License

[[top](#sections)]

* Free software: MIT license

