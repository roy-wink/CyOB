=========================
Quantum chemical packages
=========================

As briefly discussed in the introduction, CyOB uses two packes created by Dr. Ivo Filot for the quantum chemical calculations.
On this page, a bit more information about these packages can be found. How to install these packages can be found on the 
:ref:`installation` page.

CyOB needs at least PyQInt version 0.7.2.2 and PyTessel version 0.1.2. In order to check the version of the installed packages, 
type the following command in your Spyder console::

	conda list

PyQInt
------
PyQInt is a Python package for calculating one- and two-electron integrals as encountered in electronic structure calculations. 
Since integral evaluation can be quite computationally intensive, they are programmed in C++ and connected to Python using Cython.

More information can be found on the `GitHub page of PyQInt`_.


PyTessel
--------
PyTessel is a python package for building isosurfaces from 3D scalar fields.

More information can be found on the `GitHub page of PyTessel`_.


.. _GitHub page of PyQInt: https://github.com/ifilot/pyqint
.. _GitHub page of PyTessel: https://github.com/ifilot/pytessel 