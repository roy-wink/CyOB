def main():
    # CyOB needs to be imported **after** main call on Windows
    import sys
    sys.path.append('.')
    from cyob import CyOB as cyob
    
    # Set functionalities - do not delete, but set to False!
    do_single_plyfiles = True
    do_double_plyfiles = True
    do_coefficients = True
    do_abofile = True
    do_plots = True

    # Set some settings - do not leave empty
    name = 'H2'                 
    isovalue = 0.05
    basis_set = 'sto3g'

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
                  'isovalue': isovalue,
                  'basis_set': basis_set,
                  'name': name}

    # Call CyOB to start the calculations
    cyob.handler(input_dict)

if __name__ == '__main__':
    main()
