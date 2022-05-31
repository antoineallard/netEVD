'''
File that can be ran by anyone to directly launch the ABC and create figures containing results 
with minimal customization parameters
'''

from probabilistic_epidemic_forecasting import App
from probabilistic_epidemic_forecasting import Simulation
import matplotlib.pyplot as plt

interface = App()
interface.mainloop()
config_dict = interface.config_dict
simulation_example = Simulation(config_dict)
simulation_example.run_simulations()
simulation_example.run_undirected_contagion_network_projection()
simulation_example.save_simulation_results(keep_rejected_simulations=False)
simulation_example.refine_tolerance() 
simulation_example.generate_prior_posterior_parameter_distributions_figure()
simulation_example.generate_ABC_figure()
simulation_example.compute_and_draw_small_size_components_distributions(color='orange', s_max=40)
plt.show()