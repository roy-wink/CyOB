=======
Outputs
=======

In the same folder as the file from which you call the CyOB handler, a folder will be created with the name you have given for the molecule. 
In this folder, a file named ``input.txt`` is saved. This file contains all the input information the user gave to CyOB. 

Furthermore, for every functionality executed, another folder is created inside the aforementioned folder, containing the output files.
Below, these files are discussed, as well as how to open them.


.abo
----

The ``.abo`` file is saved in the ``abofile`` folder. These kind of files can be opened by `Managlyph`_. Managlyph is a neat little program,
designed in order to easily visualize orbitals. 

.. figure:: files/AboH2_1.png
   :width: 400
   :alt: The bonding orbital of |H2| as shown using Managlyph
   
   The bonding orbital of |H2| as shown using Managlyph

.. figure:: files/AboH2_2.png
   :width: 400
   :alt: The antibonding orbital of |H2| as shown using Managlyph
   
   The antibonding orbital of |H2| as shown using Managlyph

.. _Managlyph: https://www.managlyph.nl/


.ply
----

The polygon (``.ply``) files are saved in the folders ``single_plyfiles`` and ``double_plyfiles``. On most Windows systems, Print 3D is already
installed, which can easily open these types of files. Otherwise, `Blender`_ or `OpenCTM`_ could be used.

.. figure:: files/dplyH2_1.png
   :width: 400
   :alt: Double ``.ply`` file of the bonding orbital of |H2| shown using OpenCTM
   
   Double ``.ply`` file of the bonding orbital of |H2| shown using OpenCTM

.. figure:: files/dplyH2_2.png
   :width: 400
   :alt: Double ``.ply`` file of the antibonding orbital of |H2| shown using OpenCTM
   
   Double ``.ply`` file of the antibonding orbital of |H2| shown using OpenCTM

.. _Blender: https://www.blender.org/
.. _OpenCTM: http://openctm.sourceforge.net/


Density plots
-------------

For the density plots, ``.png`` files are stored in the ``density_plots`` folder. 
CyOB tries to auto-detect nodal planes. In this case, the plot will be blank and a prompt about this will be given in de console.
In the current version of CyOB, there is no option for choosing the colors used in these plots.

.. figure:: files/H2_figure1.png
   :width: 400
   :alt: Density plot of the bonding orbital of |H2|
   
   Density plot of the bonding orbital of |H2|

.. figure:: files/H2_figure2.png
   :width: 400
   :alt: Density plot of the antibonding orbital of |H2|
   
   Density plot of the antibonding orbital of |H2|

.. note::
   Due to the way PyQInt works, along with the fact that there are an (almost) infinity numer of solutions to the Schrödinger equation,
   there exists a chance that CyOB fails to flag a nodal plane. Always take another look at the plot and it's scale yourself.
   Luckily, false positives are very rare.

Orbital coefficients
--------------------

For the orbital coefficients, a single ``.txt`` file is stored in the ``orbital_coeff`` folder. Herein, the rows depict the atomic orbitals and the 
columns represent the molecular orbitals.





.. |H2| replace:: H\ :sub:`2`\