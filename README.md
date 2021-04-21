<img src="https://badgen.imc-tue.nl/badge/CyOB/v.0.2?icon=gitlab" />

# Purpose
Cython orbital builder (CyOB) is a python package (module) for easily calculating and visualizing molecular orbitals.

# Prerequisites
CyOB works in a Spyder environment. Spyder can be downloaded via the Anaconda website.

One needs to install the PyQInt and PyTessel packages via:

```
conda install -c ifilot pyqint
conda install -c ifilot pytessel
```

# Installation
CyOB can be installed via the Spyder console:

```
conda install -c rwink cyob
```

# Excecution and more information
Check the CyOB documentation via https://cyob.ddoc.dev/ for all important information.

# Changelog version 0.2
In version 0.2, the possibility was added to export energies of the molecular orbitals and the entire system.
More user-friendlyness was added in the event there is a mistake in the input file. Now, CyOB check whether all functionalities and settings are present. Also, redundant input exits execution, instead of simply being discarded.
The density plots got a big overhaul. The size of the plots and the scale on the axis are now dynamic and depend on the size of the molecule. Labels for the atoms can be placed in the plots using the plot_show_atoms argument. When CyOB thinks a nodal plane is present, the plot will be blank and a prompt will be given.
Finally, the documentation was addapted accordingly and spelling errors were fixed. Also, the beta-testers got their name on the appropriate place.