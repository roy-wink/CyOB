=====
Input
=====

As input, CyOB has to be given a dictionary containing a certain amount of input. This dictionary needs to contain the correct information,
formatted in a certain way. On the `GitHub page of CyOB`_, an example is found for |H2| (called ``H2.py``), which can also be found below. It is highly recommended to change this file with
the settings you need, instead of creating your own file. On this page, extra information can be found about the blocks of information in the file.

.. tip::
	Download or copy-paste the input file from GitHub, instead of copying it from down below.

.. code-block:: python
	:caption: Example of a CyOB execution script for |H2|
	
	def main():
		# CyOB needs to be imported **after** main call on Windows
		from cyob import CyOB as cyob
		
		# Set functionalities - do not delete, but set to False!
		do_single_plyfiles = False
		do_double_plyfiles = False
		do_coefficients = False
		do_energies = True
		do_abofile = True
		do_plots = True

		# Set some settings - do not leave empty
		name = 'H2'                 
		isovalue = 0.05
		basis_set = 'sto3g'
		plot_show_atoms = 'in_plane'

		# Set atoms
		# Add atoms in this format:
		#       ('element', x-pos, y-pos, z-pos, 'unit')
		# Place molecule on the xy-plane if possible.
		atom_list = [('H', 0, 0.36628549, 0, 'angstrom'),
			     ('H', 0, -0.36628549, 0, 'angstrom')
			    ]

		# Create dictionary - do not change!
		input_dict = {'do_single_plyfiles': do_single_plyfiles,
			      'do_double_plyfiles': do_double_plyfiles,
			      'do_coefficients': do_coefficients,
			      'do_abofile': do_abofile,
			      'do_plots': do_plots,
			      'atom_list': atom_list,
				  'plot_show_atoms': plot_show_atoms,
			      'isovalue': isovalue,
			      'basis_set': basis_set,
			      'name': name}

		# Call CyOB to start the calculations
		cyob.handler(input_dict)

	if __name__ == '__main__':
		main()
	

Import
------
On the first line, CyOB is imported as ``cyob``. It is important that the import is formulated in this exact way. 
Also, the import needs to be placed after the main call when running on Windows, otherwise Spyder will crash.

Input information
-----------------
After the import, the ``main()`` function starts, which first contains all the important information, creates the dictionary and finally calls CyOB.


Functionalities
~~~~~~~~~~~~~~~
The first few lines contain the `functionalities` (also called the *doables*). These functionalities determine what CyOB will calculate and create. 
All these settings have to be ``Booleans`` (True/False, with capitals) and have to be present. This means that the user should not remove the 
functionality settings, and should also not use strings or numbers as input.

The following functionalities exist:

do_single_plyfiles
	Determines whether CyOB creates polygon (.ply) files of single lobes of the molecular orbitals.
do_double_plyfiles
	Determines whether CyOB creates polygon (.ply) files of whole molecular orbitals (both lobes).
do_coefficients
	Determines whether CyOB exports the molecular orbital coefficients as .txt file.
do_energies
	Determines whether CyOB exports the energy of the molecule and energies of the molecular orbitals into a .txt file.
do_abofile
	Determines whether CyOB creates an .abo file, which can be opened by Managlyph.
do_plots
	Determines whether CyOB creates .png images of the contour plots of the molecular orbitals.

Examples of these outputs can be found on the :ref:`outputs` page.

.. tip::
	
	Think wisely about what you need. Everything costs computing power, and thus time. The bigger the molecule, the longer everything takes.
	Therefore, the more time you waste on useless calculations.


Settings
~~~~~~~~

Then, some settings have to be set. The following settings are important:

name
	This is the name of the molecule. It is used to export the desired files in a comprehensive way, as described on the :ref:`outputs` page. 
	The name has to	be a ``string``.
isovalue
	The isovalue determines how large the orbitals will be after the tessellation. The isovalue has to be a ``float or integer``.
	The lower the isovalue, the larger the shown orbitals.
basis_set
	The basis set determines the set of functions that are used to represent the electronic wave function in the Hartree–Fock calculations. The basis set
	has to be a ``string``, and there are four options: ``p321``, ``p631``, ``sto3g``, ``sto6g``.
plot_show_atoms
	This argument decides whether the plots will have labels for atoms at their respective positions. There are three options: ``none``, ``all`` and ``in_plane``.
	The first two options will place no labels, or labels for all atoms, respectively. The final option only places labels for atoms that are placed on the xy-plane,
	since this is the plane of the plot. If do_plots is set to ``False``, this argument is disregarded.


Atoms
~~~~~
Next up, the atoms are given to CyOB. Adding the atoms is quite a tedious process, and a small formatting error is easily made. Luckily, CyOB was 
constructed in such a way that it prints very clearly what goes wrong in the case something is not fully correct.

Atoms are written down as ``(tuples)`` inside a ``[list]``. Hence, for every atom, a separate tuple has to be made, which are grouped together in the list.
All tuples are separated by a comma, and for ease, you are allowed to place the tuples on separate lines. Be sure to only insert a new line after a comma.

Tuples for atoms are formatted as follows::

	('element', <x-pos>, <y-pos>, <z-pos>, 'unit')
	
Herein, the element and the unit have to be ``'strings'``. One must give the abbreviation of the element as known from the periodic table, not the atomic number.
The x-, y-, and z-positions have to be ``integers or floats``. The only allowed units are ``bohr`` (Bohr atomic units of length), ``angstrom``
(Ångström; |Å| meter) and ``pm`` (picometer; |pm| meter).

Note that the density plots show a cross section of the molecule through the xy-plane. Therefore, it is wise to place the molecule in the xy-plane 
as much as possible, since this leads to the best results.

.. important::

	Do not forget to use periods (``.``) as decimal separators, and not commas.


Dictionary
----------
Following, the dictionary of the input has to be made, in order to give this in a quick and comprehensive way to CyOB. 
Please do not change anything to this part of the code.


Start calculations
------------------
Finally, the calculations can be started. This is done by calling the main handler of CyOB via::

	cyob.handler(input_dict)
	
The input dictionary (``input_dict``) is given to the handler as an input parameter. This is validated for correctness, and then the calculations are started.

In Spyder, press ``F5`` to run the script.
In the console, CyOB will keep the user updated on the current state of the calculations.


.. |H2| replace:: H\ :sub:`2`\
.. |Å| replace::  10\ :sup:`-10`\
.. |pm| replace:: 10\ :sup:`-12`\
.. _GitHub page of CyOB: https://github.com/roy-wink/CyOB


Debug options
-------------
.. warning::

	The debug options should only be used for testing, since this could alter the output.
	Only use these options if you know what you are doing!
	
Two other flags could be added to the input dictionary:

'parent_folder': ``str``
	Changes the location of where the output files are saved. Inside the specified folder, a parent folder with the name of the molecule is still made.
	Be sure to use relative folder paths.
'always_overwrite': ``Boolean``
	If set to true, CyOB does not check whether the parent folder already exists, and hence overwrites as a standard. This is useful for a testing environment 
	where no keyboard inputs can be given via a console e.g.
