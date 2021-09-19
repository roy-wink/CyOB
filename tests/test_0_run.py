import unittest

class TestTest(unittest.TestCase):

    def test_test(self):
        print('Tests for CyOB running propperly')
        assert True

    def test_run(self):
        run_H2cyob()
        assert True
        

def run_H2cyob():
    from cyob import CyOB as cyob

    # to save time. Also, if double works, single will work as well.
    do_single_plyfiles = False
    do_double_plyfiles = True
    do_coefficients = True
    do_energies = True
    do_abofile = True
    do_plots = True

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
                  'do_energies': do_energies,
                  'do_abofile': do_abofile,
                  'do_plots': do_plots,
                  'atom_list': atom_list,
                  'plot_show_atoms': plot_show_atoms,
                  'isovalue': isovalue,
                  'basis_set': basis_set,
                  'name': name,
                  # to always skip the overwriting question
                  'always_overwrite': True,
                  # to force the parent folder
                  'parent_folder': 'tests'}

    cyob.handler(input_dict)
