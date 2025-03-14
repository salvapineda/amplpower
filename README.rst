========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - |github-actions| |codecov|
    * - package
      - |version| |wheel| |supported-versions| |supported-implementations| |commits-since|
.. |docs| image:: https://readthedocs.org/projects/amplpower/badge/?style=flat
    :target: https://readthedocs.org/projects/amplpower/
    :alt: Documentation Status

.. |github-actions| image:: https://github.com/salvapineda/amplpower/actions/workflows/github-actions.yml/badge.svg
    :alt: GitHub Actions Build Status
    :target: https://github.com/salvapineda/amplpower/actions

.. |codecov| image:: https://codecov.io/gh/salvapineda/amplpower/branch/main/graphs/badge.svg?branch=main
    :alt: Coverage Status
    :target: https://app.codecov.io/github/salvapineda/amplpower

.. |version| image:: https://img.shields.io/pypi/v/amplpower.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/amplpower

.. |wheel| image:: https://img.shields.io/pypi/wheel/amplpower.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/amplpower

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/amplpower.svg
    :alt: Supported versions
    :target: https://pypi.org/project/amplpower

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/amplpower.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/amplpower

.. |commits-since| image:: https://img.shields.io/github/commits-since/salvapineda/amplpower/v0.0.0.svg
    :alt: Commits since latest release
    :target: https://github.com/salvapineda/amplpower/compare/v0.0.0...main



.. end-badges

AMPL package for power systems

* Free software: MIT license

Installation
============

::

    pip install amplpower

You can also install the in-development version with::

    pip install https://github.com/salvapineda/amplpower/archive/main.zip


Documentation
=============


https://amplpower.readthedocs.io/


Development
===========

To run all the tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
