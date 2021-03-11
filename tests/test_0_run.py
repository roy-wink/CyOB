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
    do_abofile = True
    do_plots = True

    name = 'H2'
    isovalue = 0.05
    basis_set = 'sto3g'

    atom_list = [('H', 0, 0.36628549, 0, 'angstrom'),
                 ('H', 0, -0.36628549, 0, 'angstrom')
                ]

    input_dict = {'do_single_plyfiles': do_single_plyfiles,
                  'do_double_plyfiles': do_double_plyfiles,
                  'do_coefficients': do_coefficients,
                  'do_abofile': do_abofile,
                  'do_plots': do_plots,
                  'atom_list': atom_list,
                  'isovalue': isovalue,
                  'basis_set': basis_set,
                  'name': name,
                  # to skip the overwriting question
                  'always_overwrite': True,
                  # to force the parent folder
                  'parent_folder': 'tests'}

    cyob.handler(input_dict)
