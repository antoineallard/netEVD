import sys
sys.path.insert(1, '../probabilistic_epidemic_forecasting/')

from probabilistic_epidemic_forecasting.refine_results import PostProd, refine_tolerance, load_config_dict

folder_path = 'results'
raw_path = 'validation_data/time_series.txt'
simulation_path = folder_path + '/refined_simulation_results.txt'
new_path = folder_path + '/newly_refined_simulation_results.txt'
tol = {15:(30,30), 30:(30,30), 45:(30,30), 60:(30,30), 75:(30,30), 90:(30,30)}

refine_tolerance(raw_path, new_path, simulation_path, tol)
config_dict = load_config_dict(folder_path)
postprod = PostProd(config_dict, new_path)
postprod.compute_and_draw_small_size_components_distributions(color='blue', s_max=100, normalized=True)
