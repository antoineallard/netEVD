'''
Main file used to run the ABC and draw figures presenting results.
'''

import os
import pickle
import datetime
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.stats import nbinom, gamma, norm, loguniform
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import sys
import pickle

# local imports
from .other_functions import *
from . import excess_degree_distribution as EDD


class InvalidDistribution(Exception):
    def __init__(self, variable, ):
        self.variable = variable
        self.message = 'This variable has not been assigned a valid distribution: ' + variable
        super().__init__(self.message)


class Simulation:
    '''
    Main class of the project used to compute the ABC, the small comp sized distribution, the standard epidemic size projection
    and produce figures to illustrates the results.
    '''


    # ============================================ #
    #       Initialization of parameters           #
    # ============================================ #


    def __init__(self, config_dict):

        self.seed = config_dict['seed']                                                 # Number of initial cases
        self.max_cases = config_dict['max_cases']                                       # Maximum number of cases
        self.max_time = config_dict['max_time']                                         # Maximum time of the simulation
        self.n_simulation = config_dict['N_simulations']                                # Number of simulation
        self.variable_couples = {}                                                      # Dictionnary containing all the variable couple
        self.degree_dist = config_dict['excess_outward_degree_distributions']           # Excess Degree Distribution of network
        self.tolerances = config_dict['tolerances']                                     # Refined windows of tolerance
        self.configuration = config_dict                                                # Rest of dictionnary for param distributions
        self.data_path = config_dict['data_path']
        self.real_epicurve = get_epicurve(pd.read_csv(self.data_path))
        self.folder = config_dict['folder_path']
        self.folder_name = config_dict['folder_name']
        if self.folder_name is None:
            self.folder_name = "/simulation_{}".format(str(datetime.datetime.now())[0:19].replace(":", "_").replace(' ', '-'))
        else:
            self.folder_name = "/" + self.folder_name

        #EDD_obj : Excess Degree Distribution object
        self.p0 = config_dict['custom_p0']
        if self.degree_dist == 'NBinom':
            self.EDD_obj = EDD.NBinom(custom_p0=self.p0)
        if self.degree_dist == 'Binom':
            self.EDD_obj = EDD.Binom(custom_p0=self.p0)
        if self.degree_dist == 'ShiftedPowExpcut':
            self.EDD_obj = EDD.ShiftedPowExpcut(custom_p0=self.p0)
        if self.degree_dist == 'ShiftedPowerLaw':
            self.EDD_obj = EDD.ShiftedPowerLaw(custom_p0=self.p0)
        if self.degree_dist == 'Poisson':
            self.EDD_obj = EDD.Poisson(custom_p0=self.p0)
        if self.degree_dist == 'Exponential':
            self.EDD_obj = EDD.Exponential(custom_p0=self.p0)

        self.param_names = self.EDD_obj.params
        self.n_params = self.EDD_obj.num_param
        self.draw_variables()

        self.refined_results_path = None
        self.simulation_array = None
        self.post_config_array = None
        self.projection_array = None
        self.set_path()
        self.temporary_data_array = None
        self.dataframe = None



        self.create_folder()
        self.save_dict(config_dict)

    def set_path(self):
        if self.folder is not None:
            self.path = self.folder + self.folder_name
        else:
            self.path = "../Epidemic_forecasting_data" + self.folder_name
        self.configuration['path']= self.path

    def create_folder(self):
        os.makedirs(self.path)

    def save_dict(self, config_dict):
        with open(self.path+'/config_dict.pickle', 'wb') as handle:
            pickle.dump(self.configuration, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # ================================ #
    # Approximate Bayesian Computation #
    # ================================ #

    def draw_variables(self, N=10000):
        '''
        Assigns prior distribution to param_A(R_0/kappa), param_B(k/tau), sigma, alpha and gamma
        '''
        #param_A
        if self.configuration['param_A']['distribution'] == 'gamma':
            self.variable_couples['param_A'] = gamma.rvs(a=self.configuration['param_A']['gamma']['a'],
                                                    scale=self.configuration['param_A']['gamma']['b'], size=N)
        elif self.configuration['param_A']['distribution'] == 'normal':
            self.variable_couples['param_A'] = np.abs(norm.rvs(loc=self.configuration['param_A']['normal']['mean'],
                                                    scale=self.configuration['param_A']['normal']['std'], size=N))
        elif self.configuration['param_A']['distribution'] == 'uniform':
            self.variable_couples['param_A'] = np.random.uniform(low=self.configuration['param_A']['uniform']['low'],
                                                           high=self.configuration['param_A']['uniform']['high'], size=N)
        elif self.configuration['param_A']['distribution'] == 'loguniform':
            self.variable_couples['param_A'] = loguniform.rvs(self.configuration['param_A']['loguniform']['low'],
                                                              self.configuration['param_A']['loguniform']['high'], size=N)
        else:
            raise InvalidDistribution('param_A')

        #param_B
        if self.configuration['param_B']['distribution'] == 'uniform':
            self.variable_couples['param_B'] = np.random.uniform(low=self.configuration['param_B']['uniform']['low'],
                                                           high=self.configuration['param_B']['uniform']['high'], size=N)
        elif self.configuration['param_B']['distribution'] == 'normal':
            self.variable_couples['param_B'] = np.abs(norm.rvs(loc=self.configuration['param_B']['normal']['mean'],
                                                    scale=self.configuration['param_B']['normal']['std'], size=N))
        elif self.configuration['param_B']['distribution'] == 'gamma':
            self.variable_couples['param_B'] = gamma.rvs(a=self.configuration['param_B']['gamma']['a'],
                                                    scale=self.configuration['param_B']['gamma']['b'], size=N)

        elif self.configuration['param_B']['distribution'] == 'loguniform':
            self.variable_couples['param_B'] = loguniform.rvs(self.configuration['param_B']['loguniform']['low'],
                                                              self.configuration['param_B']['loguniform']['high'], size=N)
        elif self.degree_dist == 'Binom':
            self.variable_couples['param_B'] = np.around(self.variable_couples['param_B'])
        else:
            raise InvalidDistribution('param_B')

        #sigma
        if self.configuration['sigma']['distribution'] == 'uniform':
            self.variable_couples['sigma'] = np.random.uniform(low=self.configuration['sigma']['uniform']['low'],
                                                               high=self.configuration['sigma']['uniform']['high'], size=N)
        elif self.configuration['sigma']['distribution'] == 'loguniform':
            self.variable_couples['sigma'] = loguniform.rvs(self.configuration['sigma']['loguniform']['low'],
                                                              self.configuration['sigma']['loguniform']['high'], size=N)
        elif self.configuration['sigma']['distribution'] == 'gamma':
            self.variable_couples['sigma'] = gamma.rvs(a=self.configuration['sigma']['gamma']['a'],
                                                    scale=self.configuration['sigma']['gamma']['b'], size=N)
        elif self.configuration['sigma']['distribution'] == 'normal':
            self.variable_couples['sigma'] = np.abs(norm.rvs(loc=self.configuration['sigma']['normal']['mean'],
                                                    scale=self.configuration['sigma']['normal']['std'], size=N))
        #Param_C
        if self.configuration['param_C']['distribution'] == 'uniform':
            self.variable_couples['param_C'] = np.random.uniform(low=self.configuration['param_C']['uniform']['low'],
                                                               high=self.configuration['param_C']['uniform']['high'], size=N)
        elif self.configuration['param_C']['distribution'] == 'loguniform':
            self.variable_couples['param_C'] = loguniform.rvs(self.configuration['param_C']['loguniform']['low'],
                                                              self.configuration['param_C']['loguniform']['high'], size=N)
        elif self.configuration['param_C']['distribution'] == 'gamma':
            self.variable_couples['param_C'] = gamma.rvs(a=self.configuration['param_C']['gamma']['a'],
                                                    scale=self.configuration['param_C']['gamma']['b'], size=N)
        elif self.configuration['param_C']['distribution'] == 'normal':
            self.variable_couples['param_C'] = np.abs(norm.rvs(loc=self.configuration['param_C']['normal']['mean'],
                                                    scale=self.configuration['param_C']['normal']['std'], size=N))

        #gamma_shape
        if self.configuration['gamma_shape']['distribution'] == 'uniform':
            self.variable_couples['gamma_shape'] = np.random.uniform(low=self.configuration['gamma_shape']['uniform']['low'],
                                                               high=self.configuration['gamma_shape']['uniform']['high'], size=N)
        elif self.configuration['gamma_shape']['distribution'] == 'loguniform':
            self.variable_couples['gamma_shape'] = loguniform.rvs(self.configuration['gamma_shape']['loguniform']['low'],
                                                              self.configuration['gamma_shape']['loguniform']['high'], size=N)
        elif self.configuration['gamma_shape']['distribution'] == 'gamma':
            self.variable_couples['gamma_shape'] = gamma.rvs(a=self.configuration['gamma_shape']['gamma']['a'],
                                                    scale=self.configuration['gamma_shape']['gamma']['b'], size=N)
        elif self.configuration['gamma_shape']['distribution'] == 'normal':
            self.variable_couples['gamma_shape'] = np.abs(norm.rvs(loc=self.configuration['gamma_shape']['normal']['mean'],
                                                    scale=self.configuration['gamma_shape']['normal']['std'], size=N))


    def run_simulations(self):
        '''
        Computes the ABC according to parameters given by the input dictionnary

        '''
        print('Running ABC simulations...')
        self.temporary_data_array = np.zeros((10000, self.max_time+self.n_params), dtype=np.float32)  # create temporary array for result storing
        i = 0
        for simulation in tqdm(range(self.n_simulation), colour='red'):

            # Initialize each parameters
            params = [self.variable_couples['param_A'][i],self.variable_couples['param_B'][i],self.variable_couples['param_C'][i],
                        self.variable_couples['sigma'][i],self.variable_couples['gamma_shape'][i]]
            self.EDD_obj.set_params(params)




            # Initialize variables holding number of new cases and lists of infected times (total and new)
            new_cases = self.seed
            new_wave_times = [0] * self.seed
            all_infected_times = new_wave_times

            #Loop until we reached either 0 new cases or maximum cases
            while 0 < new_cases and len(all_infected_times) < self.max_cases:

                #Secondary cases from EDD
                secondary = self.EDD_obj.draw_random(new_cases)

                # Initialize lists of transmitters times and transmitted times i.e the times (in days) at which each transmitters were infected relative to
                # the begginning of the epidemic and the times (in days) at which each new infected were infected relative to their infectors
                rate = self.EDD_obj.gamma_shape / self.EDD_obj.sigma
                transmitted_times = gamma.rvs(a=self.EDD_obj.gamma_shape, scale=1 / rate, size=np.sum(secondary))
                transmitters_times = []
                for j, nb_receiver in enumerate(secondary):
                    transmitters_times += [new_wave_times[j]]*nb_receiver
                absolute_transmitted_times = np.array(transmitted_times)+np.array(transmitters_times)

                #Discard times overs maximum time studied
                new_wave_times = list(absolute_transmitted_times[absolute_transmitted_times < self.max_time])
                new_cases = len(new_wave_times) #Why the -1 (investigation needed)

                #Finally add new generation to the total number of cases and start over
                all_infected_times += new_wave_times

            histogram, bin_edges = np.histogram(all_infected_times, bins=np.arange(self.max_time + 1)) #Histogramme des temps d'infection
            epicurve = list(np.cumsum(histogram))

            self.temporary_data_array[i] = self.EDD_obj.params_values + epicurve


            if (i+1) % 10000 == 0 and i > 0:
                # saves and reset the temporay array every 10000 simulation, to optimize simulation speed
                self.compute_and_save_within_range_simulations()
                self.save_and_reset_temporary_array()

                self.draw_variables()

                i = -1
            i += 1
        self.delete_empty_rows()
        if len(self.temporary_data_array) > 0:
            self.compute_and_save_within_range_simulations()
            self.save_and_reset_temporary_array()
        header = self.param_names + ['T' + str(i) for i in range(1, self.max_time + 1)]
        self.dataframe = pd.DataFrame(data=np.loadtxt('{}/simulation_results.txt'.format(self.path)),columns=header)


    def delete_empty_rows(self):
        row_to_delete = []
        for i in range(len(self.temporary_data_array)):
            if np.all(self.temporary_data_array[i]==0):
                row_to_delete.append(i)
        for row in reversed(row_to_delete):
            self.temporary_data_array = np.delete(self.temporary_data_array, row, 0) #delete empty rows

    def save_and_reset_temporary_array(self):
        '''
        Writes all results from simulation in .txt files and resets the temporary data_array
        '''
        with open("{}/total_results.txt".format(self.path), "ab") as file:
            np.savetxt(fname=file, X=self.temporary_data_array)
        with open('{}/explored_parameters.txt'.format(self.path), 'ab') as file:
            np.savetxt(file, self.temporary_data_array[:, :self.n_params])
        self.temporary_data_array = np.zeros((10000, self.max_time + self.n_params), dtype=np.float16)

    def compute_and_save_within_range_simulations(self):
        '''
        Writes accepted results from simulation in .txt files
        '''
        data_array = self.temporary_data_array
        sim_lenght = self.max_time
        simulation_results_array = data_array[:, self.n_params: self.n_params+sim_lenght]
        selection_vector = self.selection_vector_generator(simulation_results_array, refine=False)
        with open("{}/simulation_results.txt".format(self.path), "ab") as file:
           np.savetxt(file, self.temporary_data_array[selection_vector, :self.max_time+self.n_params])

    def save_simulation_results(self, keep_rejected_simulations=False):
        '''
        Remove non accepted simulations data and transforms data in .csv
        '''
        if keep_rejected_simulations == False:
            os.remove(self.path+'/total_results.txt')

        saving_path = self.path+'/simulation_results.txt'
        self.dataframe.to_csv(saving_path)
        tqdm.write('Simulation results saved: {} '.format(saving_path))

    # ============================================ #
    # Projection calculation with standard methods #
    # ============================================ #

    def run_undirected_contagion_network_projection(self):
        '''
        Adds network style computed projections of epidemic final size to the dataframe (extra column)
        '''
        final_size_estimations = self.simulate_epidemic_size()
        self.S = final_size_estimations
        self.dataframe.insert(self.dataframe.shape[1], 'Undirected_final_size_estimation', final_size_estimations)

    def simulate_epidemic_size(self):
        '''
        Computes the projected final epidemic size (ratio) for each succesful simulation's set of parameters
        '''
        final_size_estimations = []
        for row in self.dataframe.itertuples():
            params = []
            for i in range(self.n_params):
                params.append(row[i+1])
            self.EDD_obj.set_params(params, simulation=False)
            final_size = self.EDD_obj.calculate_final_size()
            final_size_estimations.append(final_size)
        return final_size_estimations

    # ============================================ #
    #      Rejection of undesired simulations      #
    # ============================================ #



    def refine_tolerance(self):
        '''
        Redefines windows of tolerance for the simulation and writes a new .txt containing only refined succesful simulations
        '''
        tqdm.write('Rejecting invalid simulation according to refined tolerances...')
        self.refined_results_path, self.refined_dataframe = self.multiple_symmetric_selector(path=self.path)
        self.refined_dataframe.to_csv(self.refined_results_path)



    def multiple_symmetric_selector(self, path):
        '''
        Selects succesful simulations according to new windows of tolerance
        params: path(str): path to cvs file containing results
        Return: The saving path, the refined dataframe
        '''
        simulation, config = self.get_simulation_arrays()
        selection_vector = self.selection_vector_generator(simulation)
        new_data_array = self.dataframe.to_numpy()
        new_dataframe = pd.DataFrame(new_data_array[selection_vector.astype(bool), :], columns=self.dataframe.columns[:])
        saving_path = path+'/refined_simulation_results.txt'
        return saving_path, new_dataframe

    def get_simulation_arrays(self):
        '''
        Returns a tuple containing the simulation results array and then the parameters array
        '''
        data_array = self.dataframe.to_numpy()
        sim_lenght = self.max_time
        simulation_results_array = data_array[:, self.n_params: self.n_params+sim_lenght]
        configuration_array = data_array[:, :self.n_params]
        return simulation_results_array, configuration_array

    def selection_vector_generator(self, simulation_array, refine=True):
        '''
        Creates a boolean vector in which True corresponds to accepted simulations and false corresponds to rejected simulations
        Params: simulation_array: simulation array without parameters
        Return: List of boolean
        '''
        if refine:
            tolerances = self.tolerances[1]
            evaluation_days = np.array(self.tolerances[0])-1     
        else:
            tolerances = [self.tolerances[1][-1]]
            evaluation_days = np.array([self.max_time-1])
        selection_vector = np.array([],dtype=bool)
        simulation_array_at_eval_days = simulation_array[:, evaluation_days]
        real_curve_at_eval_days = self.real_epicurve[evaluation_days]
        real_array = np.array([real_curve_at_eval_days for i in range(len(simulation_array_at_eval_days))])
        tolerance_array = np.zeros((len(simulation_array_at_eval_days),len(evaluation_days)))
        for day, tol in enumerate(tolerances):
            tolerance_array[: , day] = np.isclose(simulation_array_at_eval_days[: , day], real_array[: , day],   rtol = tol/100)
        if refine:
            pbar3 = tqdm(total=len(tolerance_array),colour='red')
            for row in tolerance_array:
                pbar3.update()
                if all(row):
                    selection_vector = np.append(selection_vector, True)
                else:
                    selection_vector = np.append(selection_vector, False)
            if any(selection_vector):
                return selection_vector
            else:
                tqdm.write('No simulation accepted')
                sys.exit()
        else:
            for row in tolerance_array:
                if all(row):
                    selection_vector = np.append(selection_vector, True)
                else:
                    selection_vector = np.append(selection_vector, False)
            return selection_vector


    # ============================================ #
    #             Figure generation                #
    # ============================================ #


    def reset_and_load_results_arrays(self):
        self.simulation_array, self.post_config_array, self.projection_array = None, None, None
        self.simulation_array, self.post_config_array, self.projection_array = get_simulation_arrays_from_dataframe(self.refined_dataframe, self.n_params)

    def generate_prior_posterior_parameter_distributions_figure(self, fig_name='prior_&_posterior_param_distr.png'):
        prior_config_array = np.loadtxt(self.path + '/explored_parameters.txt')
        self.reset_and_load_results_arrays()
        draw_prior_posterior_distributions(prior_config_array, self.post_config_array, self.param_names)
        plt.savefig(self.path + '/' + fig_name, dpi=600, bbox_inches='tight')


    def generate_ABC_figure(self, fig_name='ABC_figure.png'):
        self.reset_and_load_results_arrays()
        draw_ABC_figure(self.simulation_array, self.post_config_array,  self.real_epicurve,
                        self.EDD_obj.params,self.max_time, self.tolerances, self.n_simulation)
        plt.savefig(self.path + '/' + fig_name, dpi=600, bbox_inches='tight')

    # ============================================ #
    #        Small components distribution         #
    # ============================================ #

    def compute_and_draw_small_size_components_distributions(self, color, s_max, fig_name='small_comp_distr', x_max=None, normalized=False):
        tqdm.write('Computing small components distribution...')
        self.reset_and_load_results_arrays()
        perc = self.compile_and_draw_small_component_distr(s_max=s_max, color=color, x_max=x_max, normalized=normalized)
        saving_path = self.path + '/' + fig_name
        plt.title('Degree distribution : {}'.format(self.degree_dist))

        np.savetxt(saving_path + '.txt', np.hstack((np.arange(start=1, stop=perc.shape[1]+1).reshape(perc.shape[1], 1), perc.T)),
                   header='ComponentSize     2.5     5     25     50     75     95     97.5    mean', delimiter='    ')
        plt.savefig(saving_path + '.png', dpi=600)


    def compile_and_draw_small_component_distr(self, s_max, color, x_max=None, normalized=False):
        '''
        Compile and draw small component using pyplot
        Params: s_max: maximum comp size evaluated (int), color(str): color of the drawn small comp distr,
                normalized: Boolean, divide results by [1-S] if True 
        Returns: Array containing distribution of probabilities of sizes of small-sized comp
        '''
        fig, ax = plt.subplots()
        distribution_array = np.zeros((self.post_config_array.shape[0], s_max-1))
        x_data = np.arange(1, s_max)
        pbar2 = tqdm(total=len(self.post_config_array),colour='red')
        for simul, param in enumerate(self.post_config_array):
            pbar2.update()
            #compute small comp distribution
            self.EDD_obj.set_params(param, simulation=False)
            if normalized:
                distribution_array[simul, :] = np.array(self.EDD_obj.small_comp_dist(s_max))/(1-self.S[simul])
            else:
                distribution_array[simul, :] = self.EDD_obj.small_comp_dist(s_max)
            label = None
            if simul == 1:
                label = 'simulated small comp size'
            #plot of each curves
            ax.semilogy(x_data, distribution_array[simul, :], alpha=0.15, color=color, label=label)

        #compute mean of infered parameters
        mean_distribution = np.mean(distribution_array, axis=0)

        #Label creation
        param_tilde = [self.EDD_obj.label[i]+ '{:.2f},  '.format(np.median(self.post_config_array[:, i])) for i in range(len(self.EDD_obj.label))]
        means = ''
        for i in param_tilde:
            means += i
            means += ' '
        label = f'Mean  {means} N={np.median(distribution_array.shape[0])}'

        #plot of mean curves
        ax.semilogy(x_data, mean_distribution, color='black', label=label)


        if x_max == None:
            x_limit = s_max
        else:
            x_limit = x_max
        percentile_and_mean_to_save = np.zeros((8, s_max-1))

        ymin = 0
        for index, perc in enumerate([2.5, 5, 25, 50, 75, 95, 97.5]):
            percentile_curve = np.percentile(distribution_array, q=perc, axis=0)
            if perc == 2.5:
                ymin = min(percentile_curve)
            percentile_and_mean_to_save[index, :] = percentile_curve
            ax.semilogy(x_data, percentile_curve, color='black', linestyle='dotted', linewidth=1)
            ax.annotate('{}'.format(perc), fontsize=8, xy=(x_data[int(x_limit * perc/110)], percentile_curve[int(x_limit * perc/110)]))

        percentile_and_mean_to_save[7, :] = mean_distribution
        plt.ylabel('Probability')
        if normalized:
            xlabel = "Normalized component size"
        else:
            xlabel = 'Component size'
        plt.xlabel(xlabel)
        plt.legend()
        plt.ylim([max(ymin,1e-8), 1])
        plt.xlim([0, x_limit])
        return percentile_and_mean_to_save


if __name__ == "__main__":
    print('\n')
