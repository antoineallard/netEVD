'''
File used to refine results once the simulations are done and the files are written.
'''

import pandas as pd
import numpy as np
import sys
from tqdm import tqdm
import pickle
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# local imports
from .main import Simulation
from .other_functions import draw_asymmetric_tolerance_windows, draw_ABC_figure

def get_epicurve(dataframe):
    '''
    Takes a dataframe containing epidemic data and returns the epicurve, and the epicurve_dict
    '''
    key = "Date of symptom onset "
    time_serie_string = list(dataframe[key])
    time_serie = mdates.datestr2num(time_serie_string)

    n_bins = int(time_serie[-1] - time_serie[0])
    histogram = plt.hist(time_serie, n_bins)
    plt.close()
    epicurve = np.cumsum(histogram[0])
    epicurve_dictionnary = {}
    for days_since_outbreak, nb_infected in zip(np.arange(n_bins), epicurve):
        epicurve_dictionnary[days_since_outbreak] = nb_infected
    return epicurve #, epicurve_dictionnary

def refine_tolerance(raw_data_path, new_path, simulated_data_path, tolerances):
    '''
    Redefines windows of tolerance for the simulation and writes a new .txt containing only refined succesful simulations

    Params: raw_data_path(str): path to raw data, new_path(str): path where to save newly refined data,
            simulated_data_path(str): path to refined_data.txt, tolerances(dict): indicates new refined tolerances e.g: {30:(10,50), 45:(10,30)}
    '''
    refined_results_path, refined_dataframe = multiple_asymmetric_selector(raw_data_path, new_path, simulated_data_path, tolerances)
    refined_dataframe.to_csv(refined_results_path)
    return tolerances


    
def multiple_asymmetric_selector(raw_data_path, new_path, simulated_data_path, tolerances):
    '''
    Selects succesful simulations according to new windows of tolerance
    Params: raw_data_path(str): path to raw data, new_path(str): path where to save newly refined data, 
                simulated_data_path(str): path to refined_data.txt, tolerances(dict): indicates new refined tolerances e.g: {30:(10,50), 45:(10,30)}  
    Return: The saving path, the refined dataframe
    '''
    simulation, dataframe = get_simulation_arrays(simulated_data_path)
    selection_vector = selection_vector_generator(raw_data_path, simulation, tolerances)
    new_data_array = dataframe.to_numpy()
    new_dataframe = pd.DataFrame(new_data_array[selection_vector.astype(bool), :], columns=dataframe.columns[:])
    saving_path = new_path
    return saving_path, new_dataframe

def get_simulation_arrays(path): #OK
    '''
    Returns a tuple containing the simulation results array and then the parameters array
    Param: The path to cvs file containing saved simulation results
    '''
    dataframe = pd.read_csv(path).iloc[: , 1:]
    data_array = dataframe.to_numpy()
    sim_lenght = int((dataframe.columns[-2])[1:])
    n_params = len(dataframe.columns)-sim_lenght-1
    simulation_results_array = data_array[:, n_params: n_params+sim_lenght]
    configuration_array = data_array[:, :n_params]
    return simulation_results_array, dataframe

def selection_vector_generator(raw_data_path, simulation_array, tolerances): #OK
    '''
    Creates a boolean vector in which True corresponds to accepted simulations and false corresponds to rejected simulations
    Params: raw_data_path(str): path to raw data, simulation_array: simulation array without parameters, 
            tolerances: Dict containing {day:(down_tol,up_tol), day:(down_tol, up_tol)}
    Return: List of boolean
    '''
    evaluation_days = np.array(list(tolerances.keys()))-1
    selection_vector = np.array([],dtype=bool)
    simulation_array_at_eval_days = simulation_array[:, evaluation_days]

    real_epicurve = get_epicurve(pd.read_csv(raw_data_path))
    real_curve_at_eval_days = real_epicurve[evaluation_days]
    real_array = np.array([real_curve_at_eval_days for i in range(len(simulation_array_at_eval_days))])
    
    uptol_array = np.zeros((len(simulation_array_at_eval_days),len(evaluation_days)))
    downtol_array = np.zeros((len(simulation_array_at_eval_days),len(evaluation_days)))
    pbar = tqdm(total=2*len(tolerances)+len(uptol_array[:,0]),colour='green')

    #down tol
    for day in range(len(tolerances)):
        pbar.update()
        downtol_array[: , day] = np.greater(simulation_array_at_eval_days[: , day], real_array[: , day]*(1-tolerances[evaluation_days[day]+1][0]/100))

    #up tol
    for day in range(len(tolerances)):
        pbar.update()
        uptol_array[: , day] = np.greater( real_array[: , day]*(1+tolerances[evaluation_days[day]+1][1]/100), simulation_array_at_eval_days[: , day])

    tolerance_array = np.logical_and(uptol_array,downtol_array)

    for row in tolerance_array:
        pbar.update()
        if all(row):
            selection_vector = np.append(selection_vector, True)
        else:
            selection_vector = np.append(selection_vector, False)
    if any(selection_vector):
        return selection_vector
    else: 
        tqdm.write('No simulation accepted')
        sys.exit()

def load_config_dict(path):
    """
    Loads config dict (pickle) saved from main.py
    """
    with open(path+'/config_dict.pickle', 'rb') as handle:
        return(pickle.load(handle))


class PostProd(Simulation):
    """
    Child of Simulation, this class is used to refine results, write new results and redraw figures.
    """

    def set_path(self):
        pass

    def create_folder(self):
        pass
        
    def save_dict(self, config_dict):
        pass

    def set_S(self):
        try:
            self.S = pd.read_csv(self.new_path).iloc[: , -1]
        except:
            self.S = pd.read_csv(self.path+'/refined_simulation_results.txt').iloc[: , -1]
        self.S = self.S.to_numpy()

    def __init__(self, config_dict, new_path):
        super(PostProd,self).__init__(config_dict)
        self.path = config_dict['path']
        self.new_path = new_path
        self.set_S()

    def reset_and_load_results_arrays(self):
        try:
            self.refined_dataframe = pd.read_csv(self.new_path).iloc[: , 1:]
        except:
            self.refined_dataframe = pd.read_csv(self.path+'/refined_simulation_results.txt').iloc[: , 1:]
        super(PostProd,self).reset_and_load_results_arrays()


    def generate_ABC_figure(self, fig_name='ABC_figure.png',asym_tol=False):
        self.reset_and_load_results_arrays()
        draw_ABC_figure(self.simulation_array, self.post_config_array,  self.real_epicurve,
                        self.EDD_obj.params,self.max_time, self.tolerances, self.n_simulation, asym_tol=asym_tol)
        plt.savefig(self.path + '/' + fig_name, dpi=600, bbox_inches='tight')
        


if __name__ == "__main__":
    folder_path = 'D:\STAGE\Epidemic_forecasting_data\simulation_2021-08-06-15_43_44'
    raw_path = '../Epidemic_forecasting_data/combined_database.txt'
    simulation_path = folder_path+'/refined_simulation_results.txt'
    new_path = folder_path+'/haha_refined_simulation_results.txt'
    tol = {50:(45,30),70:(45,30)}

    refine_tolerance(raw_path, new_path, simulation_path, tol)
    config_dict = load_config_dict(folder_path)
    test = PostProd(config_dict, new_path)
    test.generate_prior_posterior_parameter_distributions_figure()
    test.generate_ABC_figure(asym_tol=tol)
    test.compute_and_draw_small_size_components_distributions(color='blue', s_max=40, normalized=True, fig_name="refined_small_comp_distr")