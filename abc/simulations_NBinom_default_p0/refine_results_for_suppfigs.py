import os
import sys
sys.path.insert(1, '../probabilistic_epidemic_forecasting/')

from probabilistic_epidemic_forecasting.refine_results import PostProd, refine_tolerance, load_config_dict

folder_path = 'results'
raw_path = '../../SLED_data/lab-confirmed_database.txt'
simulation_path = folder_path + '/refined_simulation_results.txt'
new_path = folder_path + '/newly_refined_simulation_results'



suffix = '_si_every_45_days_error_30'
tol = {45:(30,30), 90:(30,30)}

refine_tolerance(raw_path, new_path + suffix + '.txt', simulation_path, tol)
config_dict = load_config_dict(folder_path)
postprod = PostProd(config_dict, new_path + suffix + '.txt')
# postprod.refined_dataframe = pd.read_csv(new_path + suffix + '.txt').iloc[: , 1:]
postprod.compute_and_draw_small_size_components_distributions(color='blue', s_max=100, normalized=True, fig_name="small_comp_distr" + suffix)
os.remove(folder_path + '/' + "small_comp_distr" + suffix + '.png')
print("")


suffix = '_si_every_45_days_error_15'
tol = {45:(15,15), 90:(15,15)}

refine_tolerance(raw_path, new_path + suffix + '.txt', simulation_path, tol)
config_dict = load_config_dict(folder_path)
postprod = PostProd(config_dict, new_path + suffix + '.txt')
# postprod.refined_dataframe = pd.read_csv(new_path + suffix + '.txt').iloc[: , 1:]
postprod.compute_and_draw_small_size_components_distributions(color='blue', s_max=100, normalized=True, fig_name="small_comp_distr" + suffix)
os.remove(folder_path + '/' + "small_comp_distr" + suffix + '.png')
print("")



suffix = '_si_every_30_days_error_30'
tol = {30:(30,30), 60:(30,30), 90:(30,30)}

refine_tolerance(raw_path, new_path + suffix + '.txt', simulation_path, tol)
config_dict = load_config_dict(folder_path)
postprod = PostProd(config_dict, new_path + suffix + '.txt')
# postprod.refined_dataframe = pd.read_csv(new_path + suffix + '.txt').iloc[: , 1:]
postprod.compute_and_draw_small_size_components_distributions(color='blue', s_max=100, normalized=True, fig_name="small_comp_distr" + suffix)
os.remove(folder_path + '/' + "small_comp_distr" + suffix + '.png')
print("")


suffix = '_si_every_30_days_error_15'
tol = {30:(15,15), 60:(15,15), 90:(15,15)}

refine_tolerance(raw_path, new_path + suffix + '.txt', simulation_path, tol)
config_dict = load_config_dict(folder_path)
postprod = PostProd(config_dict, new_path + suffix + '.txt')
# postprod.refined_dataframe = pd.read_csv(new_path + suffix + '.txt').iloc[: , 1:]
postprod.compute_and_draw_small_size_components_distributions(color='blue', s_max=100, normalized=True, fig_name="small_comp_distr" + suffix)
os.remove(folder_path + '/' + "small_comp_distr" + suffix + '.png')
print("")


suffix = '_si_every_15_days_error_30'
tol = {15:(30,30), 30:(30,30), 45:(30,30), 60:(30,30), 75:(30,30), 90:(30,30)}

refine_tolerance(raw_path, new_path + suffix + '.txt', simulation_path, tol)
config_dict = load_config_dict(folder_path)
postprod = PostProd(config_dict, new_path + suffix + '.txt')
# postprod.refined_dataframe = pd.read_csv(new_path + suffix + '.txt').iloc[: , 1:]
postprod.compute_and_draw_small_size_components_distributions(color='blue', s_max=100, normalized=True, fig_name="small_comp_distr" + suffix)
os.remove(folder_path + '/' + "small_comp_distr" + suffix + '.png')
print("")


suffix = '_si_every_15_days_error_15'
tol = {15:(15,15), 30:(15,15), 45:(15,15), 60:(15,15), 75:(15,15), 90:(15,15)}

refine_tolerance(raw_path, new_path + suffix + '.txt', simulation_path, tol)
config_dict = load_config_dict(folder_path)
postprod = PostProd(config_dict, new_path + suffix + '.txt')
# postprod.refined_dataframe = pd.read_csv(new_path + suffix + '.txt').iloc[: , 1:]
postprod.compute_and_draw_small_size_components_distributions(color='blue', s_max=100, normalized=True, fig_name="small_comp_distr" + suffix)
os.remove(folder_path + '/' + "small_comp_distr" + suffix + '.png')
print("")


suffix = '_si_every_10_days_error_30'
tol = {10:(30,30), 20:(30,30), 30:(30,30), 40:(30,30), 50:(30,30), 60:(30,30), 70:(30,30), 80:(30,30), 90:(30,30)}

refine_tolerance(raw_path, new_path + suffix + '.txt', simulation_path, tol)
config_dict = load_config_dict(folder_path)
postprod = PostProd(config_dict, new_path + suffix + '.txt')
# postprod.refined_dataframe = pd.read_csv(new_path + suffix + '.txt').iloc[: , 1:]
postprod.compute_and_draw_small_size_components_distributions(color='blue', s_max=100, normalized=True, fig_name="small_comp_distr" + suffix)
os.remove(folder_path + '/' + "small_comp_distr" + suffix + '.png')
print("")


suffix = '_si_every_10_days_error_15'
tol = {10:(15,15), 20:(15,15), 30:(15,15), 40:(15,15), 50:(15,15), 60:(15,15), 70:(15,15), 80:(15,15), 90:(15,15)}

refine_tolerance(raw_path, new_path + suffix + '.txt', simulation_path, tol)
config_dict = load_config_dict(folder_path)
postprod = PostProd(config_dict, new_path + suffix + '.txt')
# postprod.refined_dataframe = pd.read_csv(new_path + suffix + '.txt').iloc[: , 1:]
postprod.compute_and_draw_small_size_components_distributions(color='blue', s_max=100, normalized=True, fig_name="small_comp_distr" + suffix)
os.remove(folder_path + '/' + "small_comp_distr" + suffix + '.png')
print("")