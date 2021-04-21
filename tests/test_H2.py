import numpy as np
import unittest

class TestH2(unittest.TestCase):
    def test_H2abo(self):
        # Test the abofile
        abo_res = open('tests\\results\\abofile\\H2.abo', 'rb').read()
        abo_tes = open('tests\\H2\\abofile\\H2.abo', 'rb').read()

        np.testing.assert_equal(abo_res, abo_tes)

    # Testing the plots is a pain due to inaccuracies. A better way has to be found.
    # def test_H2plot(self):
    #     # Test the first plot
    #     plot_1_res = open('tests\\results\\density_plots\\H2_figure1.png', 'rb').read()
    #     plot_1_tes = open('tests\\H2\\density_plots\\H2_figure1.png', 'rb').read()

    #     np.testing.assert_equal(plot_1_res, plot_1_tes)

    #     # Test the second plot
    #     plot_2_res = open('tests\\results\\density_plots\\H2_figure2.png', 'rb').read()
    #     plot_2_tes = open('tests\\H2\\density_plots\\H2_figure2.png', 'rb').read()

    #     np.testing.assert_equal(plot_2_res, plot_2_tes)

    def test_H2ply(self):
        # Test the first plyfile
        ply_1_res = open('tests\\results\\double_plyfiles\\H2_orb1.ply', 'rb').read()
        ply_1_tes = open('tests\\H2\\double_plyfiles\\H2_orb1.ply', 'rb').read()

        np.testing.assert_equal(ply_1_res, ply_1_tes)

        ply_2_res = open('tests\\results\\double_plyfiles\\H2_orb2.ply', 'rb').read()
        ply_2_tes = open('tests\\H2\\double_plyfiles\\H2_orb2.ply', 'rb').read()

        np.testing.assert_equal(ply_2_res, ply_2_tes)

    def test_H2coeff(self):
        # Test the MO coefficients
        coeff_res = open('tests\\results\\orbital_coeff\\H2_orbital_coefficients.txt', 'r').read()
        coeff_tes = open('tests\\H2\\orbital_coeff\\H2_orbital_coefficients.txt', 'r').read()

        np.testing.assert_equal(coeff_res, coeff_tes)

    def test_H2input(self):
        #Test the input output
        input_res = open('tests\\results\\input.txt', 'r').read()
        input_tes = open('tests\\H2\\input.txt', 'r').read()

        np.testing.assert_equal(input_res, input_tes)
    
    def test_H2energy(self):
        # Test the energy output:
        energy_res = open('tests\\results\\energies\\H2_energies.txt', 'r').read()
        energy_tes = open('tests\\H2\\energies\\H2_energies.txt', 'r').read()
        
        np.testing.assert_equal(energy_res, energy_tes)
