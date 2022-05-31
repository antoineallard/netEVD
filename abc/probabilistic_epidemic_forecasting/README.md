## Probabilistic Epidemic Forecasting
Python package for early epidemic parameters inference using Approximate Bayesian Computation (ABC). 

It can also be used to create contagion networks and evaluate the giant component, and the secondary components sizes.

### Installation

Download zip file and use the package manager [pip](https://pip.pypa.io/en/stable/) to install package.
```bash
 pip3 install [path/to/extracted/zip]
```
### Usage

#### Configure, run ABC simulations, computes undirected final size projection and visualize results:

The *config* dictionary object must have precise structure, a template can be found in **detailed_example.py**. If config_dict['folder_name'] == None, the results
and the configuration of the simulation are saved in the automatically created folder **~../Epidemic_forecasting_data/simulation_[datetime]** 
```python
from probabilistic_epidemic_forecasting import Simulation
from detailed_example import config_dict

simulation_example = Simulation(config_dict)
simulation_example.run_simulations()
simulation_example.run_undirected_contagion_network_projection()
simulation_example.save_simulation_results(keep_rejected_simulations=False)  
simulation_example.refine_tolerance() 
simulation_example.generate_prior_posterior_parameter_distributions_figure()
simulation_example.generate_ABC_figure()
simulation_example.compute_and_draw_small_size_components_distributions(color='orange', s_max=100, normalized=False)
```

#### Refine results once the whole process is done

The **PostProd** class can be used to rewrite results and redraw figures.
```python
from probabilistic_epidemic_forecasting import PostProd, refine_tolerance, load_config_dict

folder_path = 'folder_in_which_results_are_stored'
raw_path = '../Epidemic_forecasting_data/combined_database.txt' #file path to raw data
simulation_path = folder_path+'/refined_simulation_results.txt'	#file path to simulated data
new_path = folder_path+'/newly_refined_simulation_results.txt'	#new file path where newly refined results will be stored
tol = {50:(45,30),70:(45,30)}					#new tolerances

refine_tolerance(raw_path, new_path, simulation_path, tol)
config_dict = load_config_dict(folder_path)
test = PostProd(config_dict)
test.generate_prior_posterior_parameter_distributions_figure()
test.generate_ABC_figure(asym_tol=tol)
test.compute_and_draw_small_size_components_distributions(color='blue', s_max=40, normalized=True, fig_name="refined_small_comp_distr")
```

##### Run simulation for new users

The **run_me.py** file can be used to launch a graphic interface helping new users to run their own simulations, with almost all customization options available.



### Reference
The network epidemiology of an Ebola epidemic: Validating contact tracing and disease surveillance<br>
L. Hébert-Dufresne, J. Bedson, L. Skrip, D. Pedi, M. F. Jalloh, B. Raulier, O. Lapointe-Gagné, J.-G. Young, A. Allard and B. M. Althouse<br>
