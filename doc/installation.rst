============
Installation
============

Cython Orbital Builder (CyOB) is a python module for easy calculation and visualization of molecular orbitals. 
CyOB works best in Spyder and needs a few other packages to be installed in order to work properly. 


Spyder
------

Spyder, one of the best Python development environments, can be installed via the `Anaconda website`_.
This link should take you directly to the installation page for Spyder stand-alone. It is also possible to install the entire
Anaconda navigator, but this is not necessairy, since you will only need Spyder.


Extra packages
--------------

When Spyder is installed, head over to its console (in the bottom right corner) and run the following commands separately::
	
	conda install -c ifilot pyqint
	conda install -c ifilot pytessel

This should install the PyQInt and PyTessel packages. More information about these packages can be found in the :ref:`Packages<Quantum chemical packages>` page.

To be sure you have the newest version of some other packages on which CyOB depends (``Cython``, ``NumPy`` and ``Matplotlib``), you could optionally enter the 
following command to update those packages::

	conda install cython numpy matplotlib


CyOB
----

CyOB is installed in a similar fashion as described above. Run the following command in the Spyder console::

	conda install -c rwink cyob

If CyOB is already installed, but you want to upgrade to the newest version (at least version 0.1.2 is required),
use the following command::

	conda update -c rwink cyob



.. _Anaconda website: https://docs.spyder-ide.org/current/installation.html#standalone-installers-ref