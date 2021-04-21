from pyqint import PyQInt, pyqint, Molecule, HF
from pytessel import PyTessel
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import json
import time
import sys
import os
import re


def _version():
    """
    Function to return the version of CyOB in ``string`` format.
    
    """
    cyob_version = '0.2'
    return cyob_version


def _ang2bohr():
    """
    Function to return the conversion factor from Angstrom to Bohr a.u. in
    ``float`` format.
    
    """
    ang2bohr = 1.88973
    return ang2bohr


def handler(input_dict):
    """
    The main handler for Cython Orbital Builder. Takes the information from the
    input dictionary to determine what functions to call.

    """

    # start time
    start = time.time()

    # Check for valid input
    if not validity_check(input_dict):
        sys.exit('Invalid input. Unknown error.')

    print('CyOB execution for %s started' % input_dict['name'])

    # Read functionalities
    try:
        do_single_plyfiles = input_dict['do_single_plyfiles']
        do_double_plyfiles = input_dict['do_double_plyfiles']
        do_coefficients = input_dict['do_coefficients']
        do_energies = input_dict['do_energies']
        do_abofile = input_dict['do_abofile']
        do_plots = input_dict['do_plots']
    except Exception as e:
        sys.exit('Missing input functionality: %s' % e)

# perform HF calculations
    print('Hartree-Fock calculations for %s started' % input_dict['name'])
    cgfs, coeff, energy, energies = calculate_molecule(input_dict['name'], 
                                                       input_dict['atom_list'], 
                                                       input_dict['basis_set'])

# make a dict of tesselated things for efficiency
    print('Making tesselation dictionary')
    tesselation_needed = do_single_plyfiles or do_double_plyfiles or do_abofile
    if tesselation_needed:
        orbital_info = tesselate_all(cgfs, coeff, input_dict['isovalue'])

# start doing things
    # create parent folder
    if 'parent_folder' in input_dict:
        parent_folder = os.path.join(input_dict['parent_folder'], input_dict['name'])
    else:
        parent_folder = os.path.join(os.getcwd(), input_dict['name'])
    print('Working in parent folder: %s' % parent_folder)

    if not os.path.exists(parent_folder):
        os.mkdir(parent_folder)
    else:
        if 'always_overwrite' not in input_dict or not input_dict['always_overwrite']:
            print('\nWARNING: PREVIOUS RESULTS MIGHT BE OVERWRITTEN')
            cont = input('Continue? (y/n)\n')
            if cont not in ['y', 'Y', 'yes', 'Yes', 'j', 'J', 'ja', 'Ja']:
                sys.exit('Manually stopped to prevent overwriting')

    # write input info to file
    with open(parent_folder + '\\input.txt', 'w') as info_file:
        for line in input_dict:
            print((line + ': ' + str(input_dict[line])), file=info_file)

    # make .ply files
    print('Creating .ply files for single lobes of %s' % input_dict['name'])
    if do_single_plyfiles:
        # generate folder
        sp_folder = os.path.join(parent_folder, 'single_plyfiles')
        if not os.path.exists(sp_folder):
            os.mkdir(sp_folder)

        for orb in orbital_info:
            filename = '%s\\%s_orb%s.ply' % (sp_folder, input_dict['name'], str(orb['orbital']) + '_pos')
            single_polygon(filename, orb, 'positive')
            filename = '%s\\%s_orb%s.ply' % (sp_folder, input_dict['name'], str(orb['orbital']) + '_neg')
            single_polygon(filename, orb, 'negative')
        print('\tDone!')

    if do_double_plyfiles:
        print('Creating .ply files for full orbitals of %s' % input_dict['name'])
        # generate folder
        dp_folder = os.path.join(parent_folder, 'double_plyfiles')
        if not os.path.exists(dp_folder):
            os.mkdir(dp_folder)

        for orb in orbital_info:
            filename = '%s\\%s_orb%s.ply' % (dp_folder, input_dict['name'], str(orb['orbital']))
            double_polygon(filename, orb)
        print('\tDone!')

    if do_abofile:
        print('Creating .abo file for %s' % input_dict['name'])
        abo_folder = os.path.join(parent_folder, 'abofile')
        if not os.path.exists(abo_folder):
            os.mkdir(abo_folder)

        build_abo(orbital_info, abo_folder + '\\' + input_dict['name'], input_dict['atom_list'], input_dict['basis_set'])
        print('\tDone!')

    if do_plots:
        print('Creating density plots for %s' % input_dict['name'])
        p_folder = os.path.join(parent_folder, 'density_plots')
        if not os.path.exists(p_folder):
            os.mkdir(p_folder)

        filename = os.path.join(p_folder, input_dict['name'])
        make_plots(input_dict['atom_list'], cgfs, coeff, filename, input_dict['plot_show_atoms'])
        print('\tDone!')

    if do_coefficients:
        print('Exporting orbital coefficients for %s' % input_dict['name'])
        c_folder = os.path.join(parent_folder, 'orbital_coeff')
        if not os.path.exists(c_folder):
            os.mkdir(c_folder)

        filename = '%s\\%s_orbital_coefficients.txt' % (c_folder, input_dict['name'])
        np.savetxt(filename, coeff)
        print('\tDone!')
    
    if do_energies:
        print('Exporting energies for %s' % input_dict['name'])
        e_folder = os.path.join(parent_folder, 'energies')
        if not os.path.exists(e_folder):
            os.mkdir(e_folder)

        filename = '%s\\%s_energies.txt' % (e_folder, input_dict['name'])
        export_energies(energy, energies, filename)
        print('\tDone!')


    runtime = time.time() - start
    print('Finished CyOB execution for %s in %s seconds' % (input_dict['name'], '{:.4f}'.format(runtime)))


def validity_check(input_dict):
    """
    Function in order to check whether the given input dictionary is valid.

    Exits the program if a formatting error is found. Returns ``True`` if the
    dictionary is valid.

    :param input_dict: The input ``dictionary``.
    
    """

    # check doables
    all_dos = ['do_single_plyfiles', 'do_double_plyfiles', 'do_coefficients',
               'do_energies', 'do_abofile', 'do_plots']
    do_list = [do for do in input_dict if re.match('do_*', do)]
    for do in do_list:
        if not isinstance(input_dict[do], bool):
            sys.exit('Formatting error found in input functionality %s (not a boolean)' % do)
        if do not in all_dos:
            sys.exit('Unknown functionality found in input: %s' % do)

    # check if all settings are present
    all_settings = ['name', 'isovalue', 'basis_set']
    for setting in all_settings:
        if setting not in all_settings:
            sys.exit('Missing setting argument: %s' % setting)

    # check isovalue
    try:
        isovalue = float(input_dict['isovalue'])
    except ValueError:
        sys.exit('Incorrect isovalue found (not a number): %s' % isovalue)

    if input_dict['isovalue'] <= 0:
        sys.exit('Please use a positive isovalue')

    # check name
    if not isinstance(input_dict['name'], str):
        sys.exit('Please provide name as string.')

    # check basis set
    if input_dict['basis_set'] not in ('p321', 'p631', 'sto3g', 'sto6g'):
        sys.exit('Unknown basis set: %s' % input_dict['basis_set'])

    # check atoms plot flag (does not have to be present if not plotting)
    if input_dict['do_plots']:
        if 'plot_show_atoms' not in input_dict:
            sys.exit('Missing setting argument: plot_show_atoms')
        if input_dict['plot_show_atoms'] not in ('all', 'none', 'in_plane'):
            sys.exit('Unknown argument for plot_show_atoms: %s' % input_dict['plot_show_atoms'])

    # check atoms
    if not isinstance(input_dict['atom_list'], list):
        sys.exit('Formatting error found in atoms (not a list)')
    elif len(input_dict['atom_list']) == 0:
        sys.exit('No atoms in input')

    # first locate basis set to check for valid atoms
    basis_loc = '%s\\basis\\%s.json' % (os.path.dirname(pyqint.__file__), input_dict['basis_set'])
    basis_data = json.load(open(basis_loc))

    for atom in input_dict['atom_list']:
        if not isinstance(atom, tuple):
            sys.exit('Formatting error found in atoms (not a tuple): \n%s' % atom)

        if len(atom) not in [4, 5]:
            sys.exit('Unknown formatting error in atoms: %s' % atom)

        element, x, y, z = atom[:4]
        if element not in basis_data:
            sys.exit('Unknown element type: %s' % element)

        for value in [x, y, z]:
            try:
                value = float(value)
            except ValueError:
                sys.exit('Incorrect position found for %s (not a number): %s' % (element, value))

        if len(atom) == 5:
            unit = atom[4]
            if unit not in ['bohr', 'angstrom', 'pm']:
                sys.exit('Unknown unit found in atoms: %s' % unit)

    return True


def calculate_molecule(name, atoms, basis):
    """
    Function for performing the Hartree-Fock calculations on the molecule
    using the PyQInt package.

    Outputs the contracted Gaussian basis functions and the orbital
    coefficients.

    :param name: The name of the molecule (``str``).
    :param atoms: The ``list`` of atom ``tuples``.
    :param basis: The basis set (``str``).
    
    """

    mol = Molecule(_name=name)
    
    for atom in atoms:
        element = str(atom[0])
        x, y, z = float(atom[1]), float(atom[2]), float(atom[3])
        if len(atom) == 5:
            unit = atom[4]
            if unit == 'pm':
                x, y, z = x / 100, y / 100, z / 100
                unit = 'angstrom'
            mol.add_atom(element, x, y, z, unit)
        else:
            mol.add_atom(element, x, y, z)

    try:
        result = HF().rhf(mol, basis)
    except Exception as e:
        print('An error occured: %s' % e)

    return result['cgfs'], result['orbc'], result['energy'], result['energies']


def tesselate(cgfs, coeff, isovalue=0.01):
    """
    Function for building isosurfaces (tessellation) of the given orbital,
    using the PyTessel and PyQInt packages.

    Outputs arrays for the vertices, normals and indices (faces).

    :param cgfs: Contracted Gaussion basis functions of the molecule.
    :param coeff: Coefficients of the desired molecular orbital (``list``).
    :param isovalue: The desired isovalue (``float`` or ``int``, standard is 0.01).
    
    """

    # initiate
    integrator = PyQInt()
    pytessel = PyTessel()

    # build grid
    div = 100
    length = 5
    ang2bohr = _ang2bohr()
    grid = integrator.build_rectgrid3d(-length, length, div)
    unitcell = np.multiply(np.diag(np.ones(3) * 2 * length), 1 / ang2bohr)
    # To go from bohr to A, since this is the unit for managlyph
    scalarfield = np.reshape(integrator.plot_wavefunction(grid, coeff, cgfs), (div, div, div))

    vertices, normals, indices = \
        pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, \
                                unitcell.flatten(), isovalue)

    return vertices, normals, indices


def tesselate_all(cgfs, coeff, isovalue=0.01):
    """
    Function for creating a dictionaries for tessellated surfaces for all
    molecular orbitals. This increases the efficiency of CyOB.
    
    Outputs a list with a dict per MO. The dicts contain information for
    the positive and negative lobes.
    
    :param cgfs: Contracted Gaussion basis functions of the molecule.
    :param coeff: Coefficients of the desired molecular orbital (``list``).
    :param isovalue: The desired isovalue (``float`` or ``int``, standard is 0.01).
    
    """

    orbital_info = []
    for i in range(len(coeff)):
        orb = coeff[:,i]
        p_vertices, p_normals, p_indices = tesselate(cgfs, orb, isovalue)
        n_vertices, n_normals, n_indices = tesselate(cgfs, orb, -isovalue)
        orb_dict = {'orbital': (i + 1),
                    'p_vertices': p_vertices,
                    'p_normals': p_normals,
                    'p_indices': p_indices,
                    'n_vertices': n_vertices,
                    'n_normals': n_normals,
                    'n_indices': n_indices,}
        orbital_info.append(orb_dict)
    return orbital_info
    

def single_polygon(filename, orbital_dict, pos_neg='positive'):
    """
    Function for creating polygon files of a separate lobe of an orbital.
    
    :param filename: Filename under which the created file will be saved \
        (``str``, including suffix)
    :param orbital_dict: ``Dictionary`` of the information of the desired \
        molecular orbital as created by :py:func:`CyOB.tesselate_all`.
    :param pos_neg: Flag to determine whether the positive- or negative lobe \
        has to be made into a polygon file (``'positive'`` or ``'negative'``)
    
    """

    if pos_neg not in ['positive', 'negative']:
        sys.exit('Unknown handler in def single_polygon(): ' % pos_neg)

    pytessel = PyTessel()
    if pos_neg == 'positive':
        vertices = orbital_dict['p_vertices']
        normals = orbital_dict['p_normals']
        indices = orbital_dict['p_indices']
    if pos_neg == 'negative':
        vertices = orbital_dict['n_vertices']
        normals = orbital_dict['n_normals']
        indices = orbital_dict['n_indices']
    pytessel.write_ply(filename, vertices, normals, indices)


def double_polygon(filename, orbital_dict):
    """
    Function for creating polygon files of a both lobe of an orbital.
    Vertices, normals and indices are combined before creating the file.
    
    :param filename: Filename under which the created file will be saved\
        (``str``, including suffix)
    :param orbital_dict: ``Dictionary`` of the information of the desired\
        molecular orbital as created by :py:func:`CyOB.tesselate_all`.
    
    """

    pytessel = PyTessel()

    p_vertices = orbital_dict['p_vertices']
    p_normals = orbital_dict['p_normals']
    p_indices = orbital_dict['p_indices']

    n_vertices = orbital_dict['n_vertices']
    n_normals = orbital_dict['n_normals']
    n_indices = orbital_dict['n_indices']

    vertices = np.matrix(np.vstack([p_vertices, n_vertices]), dtype=np.float32)
    normals = np.matrix(np.vstack([p_normals, n_normals]), dtype=np.float32)
    indices = np.array(np.concatenate((p_indices, [x+len(p_vertices) for x in n_indices])), dtype=np.uint32)
    pytessel.write_ply(filename, vertices, normals, indices)



def build_abo(orbital_info, name, atoms, basis):
    """
    Function for building the file to be opened in Managlyph (``.abo``).
    
    :param orbital_list: ``List`` of ``dictionaries`` containing the information\
        about the molecular orbitals, as created by :py:func:`CyOB.tesselate_all`.
    :param name: The name of the molecule, used to determine the name of the\
        file (``str``).
    :param atoms: The ``list`` of atom ``tuples``.
    :param basis: The basis set (``str``).
    
    """

    outfile = ('%s.abo' % name)
    f = open(outfile, 'wb')

    # change element names to atomic numbers and count electrons
    basis_loc = '%s\\basis\\%s.json' % (os.path.dirname(pyqint.__file__), basis)
    basis_data = json.load(open(basis_loc))

    atoms_num = []
    nelec = 0
    for atom in atoms:
        atom = list(atom)
        atom[0] = basis_data[atom[0]]['atomic_number']
        atoms_num.append(tuple(atom))
        nelec += atom[0]

    # set electron counter
    el_counter = 0

    # write number of frames
    f.write(len(orbital_info).to_bytes(2, byteorder='little'))

    # loop over orbitals
    for k in range(len(orbital_info)):
        # write frame index (k)
        f.write(k.to_bytes(2, byteorder='little'))

        # write description
        description = 'orbital ' + str(k + 1)
        f.write(len(description).to_bytes(2, byteorder='little'))
        f.write(bytearray(description, encoding='utf8'))

        # write atoms
        f.write(len(atoms).to_bytes(2, byteorder='little'))
        # loop over atoms
        for atom in atoms_num:
            atom = list(atom)
            # pm has to be converted to angstrom
            if len(atom) > 4 and atom[4] == 'pm':
                for pos in range(1,4):
                    atom[pos] = atom[pos] / 100

            # bohr has to be converted to angstrom
            if len(atom) > 4 and atom[4] == 'bohr':
                ang2bohr = _ang2bohr()
                for pos in range(1,4):
                    atom[pos] = atom[pos] / ang2bohr

            # write atom information
            f.write(atom[0].to_bytes(1, byteorder='little'))
            f.write(np.array(atom[1:4], dtype=np.float32).tobytes())

        # determine colors
        el_counter += 2
        if el_counter <= nelec:
            # occupied
            color_pos = (0, 0, 1, 0.8)
            color_neg = (1, 0, 0, 0.8)
        else:
            color_pos = (0, 1, 0.8, 0.7)
            color_neg = (1, 0, 1, 0.8)

        # determine models (pos and neg)
        # R, G, B, a, vertices, normals, indices
        orbital = orbital_info[k]
        positive = color_pos + (orbital['p_vertices'], orbital['p_normals'], orbital['p_indices'])
        negative = color_neg + (orbital['n_vertices'], orbital['n_normals'], orbital['n_indices'])
        models = [positive, negative]

        # write number of models
        f.write(len(models).to_bytes(2, byteorder='little'))

        # loop over models
        for m in range(len(models)):
            # write model index
            f.write(m.to_bytes(2, byteorder='little'))

            # write model color
            f.write(np.array(models[m][0:4], dtype=np.float32).tobytes())

            # write number of vertices
            vertices = models[m][4]
            normals = models[m][5]
            f.write(len(vertices).to_bytes(4, byteorder='little'))

            # loop over vertices and normals
            for v in range(len(vertices)):
                f.write(np.array(vertices[v], dtype=np.float32).tobytes())
                f.write(np.array(normals[v], dtype=np.float32).tobytes())

            # write number of faces
            indices = models[m][6]
            f.write((int(len(indices) / 3)).to_bytes(4, byteorder='little'))
                # geen idee waarom, maar die /3 moet...
                # aantal faces = indices / 3

            # write faces
            f.write(np.array(indices, dtype=np.uint32).tobytes())

    f.close()


def make_plots(atom_list, cgfs, coeff, name, plot_show_atoms):
    """
    Function for creating the contour plots of the orbitals using the
    PyQInt package.
    
    :param atoms: The ``list`` of atom ``tuples``.
    :param cgfs: Contracted Gaussion gasis functions of the molecule.
    :param coeff: Coefficients of the desired molecular orbital (``list``).
    :param name: The name of the molecule, used to determine the name of the\
        files (str).
    :param plot_show_atoms: Argument that tells whether or not to plot atom\
        letters in the plots.
    
    """

    integrator = PyQInt()
    
    # First set all dimensions to Bohr a.u.
    for n in range(len(atom_list)):
        atom = atom_list[n]
        x, y, z = atom[1:4]
        unit = atom[4] if len(atom) == 5 else 'bohr'
        if unit == 'pm':
            x, y = x / 100, y / 100
            unit = 'angstrom'
        
        if unit == 'angstrom':
            ang2bohr = _ang2bohr()
            x, y = x * ang2bohr, y * ang2bohr
            unit = 'bohr'
        
        atom_list[n] = (atom[0], x, y, z, unit)

    # get x- and y-locations to make the best grid
    max_x_loc = np.max([np.abs(atom[1]) for atom in atom_list])
    max_y_loc = np.max([np.abs(atom[2]) for atom in atom_list])

    x_range = np.max([3 * np.sqrt(max_x_loc), 2.5])
    y_range = np.max([3 * np.sqrt(max_y_loc), 2.5])

    #build grid
    x = np.linspace(-x_range, x_range, 100)
    y = np.linspace(-y_range, y_range, 100)
    xx, yy = np.meshgrid(x,y)
    zz = np.zeros(len(x) * len(y))
    grid = np.vstack([xx.flatten(), yy.flatten(), zz]).reshape(3,-1).T

    # Turn interactive plotting off and save current state
    old_plt_query = plt.isinteractive()
    if old_plt_query:
        plt.ioff()

    for i in range(len(coeff)):
        dens = integrator.plot_wavefunction(grid, coeff[:,i], cgfs).reshape((len(y), len(x)))

        fig = plt.figure(i + 1)
        axes = fig.add_subplot(1, 1, 1)
        
        # define benchmark to spot nodal planes
        benchmark = 1e-2
        limit = max(max(abs(np.min(dens)), abs(np.max(dens))), benchmark)
        
        im = axes.imshow(dens, origin='lower', interpolation='bilinear',
            extent=[-x_range,x_range,-y_range,y_range], cmap='RdBu', vmin=-limit, vmax=limit)
            # cmap is kleur. PiYG voor paars-groen, RdBu voor blauw-rood, RdGy voor rood zwart
        axes.set_xlabel('Distance a.u. (x-axis)')
        axes.set_ylabel('Distance a.u. (y-axis)')
        divider = make_axes_locatable(axes)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        
        if limit == benchmark:
            print('\tOrbital number %i probably has a nodal plane on the xy-plane.' % (i + 1))
            print('\tRotate the molecule over the z-axis and run again to be certain.\n')
            im.remove()
        
        if plot_show_atoms in ['all', 'in_plane']:
            for atom in atom_list:
                element = atom[0]
                x_atom_loc, y_atom_loc, z_atom_loc = atom[1:4]
                if z_atom_loc == 0 or plot_show_atoms == 'all':
                    axes.text(x_atom_loc + 0.01, y_atom_loc - 0.01, element, horizontalalignment='center', 
                              verticalalignment='center', fontsize=12, color='w')
                    axes.text(x_atom_loc + 0.00, y_atom_loc - 0.00, element, horizontalalignment='center', 
                              verticalalignment='center', fontsize=12, color='k')

        filename = ('%s_figure%s.png' % (name, str(i + 1)))
        fig.savefig(filename, dpi=300)
        plt.close(fig)

    # Turn interactive plotting back on if it was like that
    if old_plt_query:
        plt.ion()


def export_energies(energy, energies, filename):
    """
    Function for exporting the system energy and orbital energies.
    
    :param energy: Energy of the system (``float``).
    :param energies: Energies of the molecular orbitals (``list``).
    :param filename: filename of the .txt file where the function exports to \
        (``str``).

    """
    
    f = open(filename, 'w')
    
    f.write('Molecule energy:\n')
    f.write('%f\n\n' % energy)

    f.write('Orbital energies:\n')
    for e in energies:
        f.write('%f\n' % e)

    f.close()
