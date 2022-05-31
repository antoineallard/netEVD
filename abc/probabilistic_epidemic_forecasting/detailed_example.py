from probabilistic_epidemic_forecasting import Simulation
import numpy as np

if __name__ == "__main__":
    config_dict = {'param_A': {'distribution': 'uniform',  # mean R0 (NBinom) or kapppa (PowExpCut)
                                                           # high R0 increases dramatically simulation time
                          # Change the parameters corresponding to the chosen distribution
                          'uniform': {'low': 0.1,
                                      'high': 4},
                          'gamma': {'a': 1.5,              # a is corresponding to shape and b is scale, see
                                    'b': 1.},              # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gamma.html for mode details
                          'normal': {'mean': 1.5,
                                     'std': 0.5}},

                   'param_B': {'distribution': 'loguniform',              # dispersion k (NBinom) or tau (PowExpCut)
                         'loguniform': {'low': 1e-4,
                                     'high': 0.5},
                         'normal': {'mean': 0.1,
                                     'std': 0.2}},

                   'param_C': {'distribution': 'uniform',                # Probability of having zero secondary infections 0 =< alpha =< 1
                             'uniform': {'low': 0.30,                  # only used for PowExpcut distr, ignored otherwise
                                         'high': 0.80}},

                   'sigma': {'distribution': 'uniform',                # Mean serial time between two infected generations [days]
                             'uniform': {'low': 0,
                                         'high': 10}},

                   'gamma_shape': {'distribution': 'uniform',          # definition of gamma_shape
                                   'uniform': {'low': 1,
                                               'high': 40.}},

                   'seed': 1,                                               # Number of initial infected persons
                   'max_cases': 1e4,                                        # Maximum cases before the simulation is stopped
                   'max_time': 100,                                         # Maximum simulation days (day 0 included)
                   'N_simulations': 5000,                                   # Total number of iteration in ABC 
                   'excess_outward_degree_distributions': 'NBinom',  # Secondary infections
                   'custom_p0': None,
                   'tolerances': [[100],[100]],                  # Refined windows of tolerance [[day,day,day...],[tol,tol,tol...]]
                   'data_path': "..\Epidemic_forecasting_data\combined_database.txt",
                   'folder_path': None,
                   'folder_name': None,
                   'notes': ''}

    simulation_example = Simulation(config_dict)
    simulation_example.run_simulations()
    simulation_example.run_undirected_contagion_network_projection()
    simulation_example.save_simulation_results(keep_rejected_simulations=False)  # Recommended
    simulation_example.refine_tolerance() 
    simulation_example.generate_prior_posterior_parameter_distributions_figure()
    simulation_example.generate_ABC_figure()
    simulation_example.compute_and_draw_small_size_components_distributions(color='orange', s_max=100, normalized=False)
    