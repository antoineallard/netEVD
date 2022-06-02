import os
import sys
sys.path.insert(1, '../probabilistic_epidemic_forecasting/')
from probabilistic_epidemic_forecasting import Simulation
import numpy as np


config_dict = {'param_A': {'distribution': 'uniform',              # R0 parameter (average number of secondary infections)
                              'uniform': {'low': 0.01,
                                          'high': 5.}},

               'param_B': {'distribution': 'loguniform',           # k parameter (dispersion coefficient)
                              'loguniform': {'low': 1e-4,
                                             'high': 10.}},

               'param_C': {'distribution': 'uniform',              # Unused parameter.
                              'uniform': {'low': 0.,
                                          'high': 1.}},

               'sigma': {'distribution': 'uniform',                # sigma parameter (average serial interval between epidemic generations)
                              'uniform': {'low': 0.,
                                          'high': 75}},

               'gamma_shape': {'distribution': 'uniform',          # alpha parameter (shape parameter for the gamma distribution)
                              'uniform': {'low': 0.,
                                          'high': 25.}},

               'seed': 1,                                          # number of initial infected individuals
               'max_cases': 1e4,                                   # maximum cases before the single simulation is stopped
               'max_time': 100,                                    # length of the simulations in days (day 0 included)
               'N_simulations': 2000000,                          # total number of ABC iteration
               'excess_outward_degree_distributions': 'NBinom',    # distribution of the number of secondary infections
               'custom_p0': None,                                  # p0 will be fixed automatically
               'tolerances': [[90],[30]],                          # rough window of tolerance [[day,day,day...],[tol,tol,tol...]]
               'data_path': "../../SLED_data/lab-confirmed_database.txt",
               'folder_path': '../rough_results/',
               'folder_name': 'simulations_NBinom_default_p0_46',
               'notes': ''}


for i in range(100):

    folder_name = 'simulations_NBinom_default_p0_{:03d}'.format(i)

    if not os.path.isdir(config_dict['folder_path'] + folder_name):

        print('\n\n\n' + folder_name + '\n')

        config_dict['folder_name'] = folder_name

        simulation = Simulation(config_dict)
        simulation.run_simulations()
        simulation.run_undirected_contagion_network_projection()
        simulation.save_simulation_results(keep_rejected_simulations=False)  # Recommended
        simulation.refine_tolerance()
        simulation.generate_prior_posterior_parameter_distributions_figure()
        simulation.generate_ABC_figure()
        simulation.compute_and_draw_small_size_components_distributions(color='orange', s_max=100, normalized=True)

        os.remove(config_dict['folder_path'] + folder_name + '/ABC_figure.png')
        os.remove(config_dict['folder_path'] + folder_name + '/explored_parameters.txt')
        os.remove(config_dict['folder_path'] + folder_name + '/prior_&_posterior_param_distr.png')
        os.remove(config_dict['folder_path'] + folder_name + '/simulation_results.txt')
        os.remove(config_dict['folder_path'] + folder_name + '/small_comp_distr.png')
        os.remove(config_dict['folder_path'] + folder_name + '/small_comp_distr.txt')
