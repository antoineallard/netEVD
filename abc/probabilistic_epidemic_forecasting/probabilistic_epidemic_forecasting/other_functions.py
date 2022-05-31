"""
File containing many useful functions that couldn't find a family
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from numpy.lib.type_check import real
import seaborn as sns
import pandas as pd
import matplotlib.dates as mdates
import mpmath as mpm
import scipy as sc
import time
from math import factorial
from numba import njit
import math

def get_epicurve(dataframe):
    '''
    Takes a dataframe containing epidemic data and returns the epicurve
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
    return epicurve


def get_simulation_arrays_from_dataframe(dataframe, n_params):
    '''
    Takes the dataframe containing all information and splits it in three: simulated epicurve, infered parameters and prediction

    Params: dataframe: complete dataframe, n_params: number of variable parameters
    '''
    data_array = dataframe.to_numpy()
    simulation_results_array = data_array[:, n_params:-1]
    configuration_array = data_array[:, :n_params]
    prediction_array = data_array[:, -1]
    return simulation_results_array, configuration_array, prediction_array


def draw_prior_posterior_distributions(prior_configuration_array, posterior_configuration_array, param_names):
    '''
    Generates figure showing prior and posterior distributions of parameters
    params: prior_config_array, osterior_configuration_array, param_names(list): contains param_names(str)
    '''
    fig, ax = plt.subplots(2, len(param_names), figsize=(9, 7))
    colors = ['#1a10ae', '#db0073', '#ff743a', '#edde55']
    ax[0, 0].set_ylabel('PRIOR', fontsize=20)
    ax[1, 0].set_ylabel('POSTERIOR', fontsize=20)
    for index, param in enumerate(param_names):
        prior_distr, post_distr = prior_configuration_array[:, index], posterior_configuration_array[:, index]
        medians = [np.median(prior_distr), np.median(post_distr)]
        sns.kdeplot(prior_distr, shade=True, ax=ax[0, index], cut=0, alpha=0.50, color=colors[index])
        sns.kdeplot(post_distr, shade=True, ax=ax[1, index], cut=0, color=colors[index])
        ax[0, index].set_title(param, fontsize=20)
        for _index, median in enumerate(medians):
            ax[_index, index].axvline(x=median, linestyle='--', color='grey')
            ax[_index, index].text(x=median + median * 0.1, y=0.9, s=r'$\tilde{x}$' + '={:.2f}'.format(median),
                transform=transforms.blended_transform_factory(ax[_index, index].transData, ax[_index, index].transAxes))
            ax[_index, index].set_yticks([])


def draw_tolerance_windows(tolerances, real_epicurve):
    '''
    Adds symetric tolerance windows to ABC_figure
    Params: tolerances(list): all windows of tolerances e.g [[30],[15]] meaning one evaluation at day 30 with 15% relative tolerance
    '''
    for i in range(len(tolerances[0])):
        plt.errorbar(tolerances[0][i]-1, real_epicurve[tolerances[0][i]-1], yerr=0.01 * tolerances[1][i] * real_epicurve[tolerances[0][i]-1], 
                        ecolor='green', alpha=0.5, elinewidth=2, fmt='none', barsabove=True, capsize=3)

def draw_asymmetric_tolerance_windows(tolerances, real_epicurve):
    '''
    Adds asymetric tolerance windows to ABC_figure
    Params: tolerances(dict): all windows of tolerances e.g {30:(15,30), 45:(5,25)}
            ,real epicurve: array containing active cases per day (real data)
    '''
    for i in tolerances:
        down = real_epicurve[i-1]*(0.01 * tolerances[i][0] )   # -1 because days are counted from 0
        up =  real_epicurve[i-1] * (0.01 * tolerances[i][1])   # idem    
        plt.errorbar(i-1, real_epicurve[i-1], yerr=([down],[up]), 
                        ecolor='green', alpha=0.5, elinewidth=2, fmt='none', barsabove=True, capsize=3)

def draw_ABC_figure(simulation_array, config_array, real_epicurve, param_names, max_time, tolerances, N_simul, asym_tol=False):
    '''
    Plots raw data and simulation results. y: Number of cases x:Number of days since outbreak
    Params: simulation_array, config_array: parameters, real_epicurve, param_names(list): contains param_names(str), max_time(int): total time of epidemic study,
            tolerances(list): all windows of tolerances e.g [[30],[15]] meaning one evaluation at day 30 with 15% relative tolerance, N_simul(int): total number of simulations
    '''
    fig = plt.figure(constrained_layout=True, figsize=(12, 6.75))
    gs = fig.add_gridspec(2,len(param_names))
    fig_ax1 = fig.add_subplot(gs[0, :])
    for simulation_curve in simulation_array:
        fig_ax1.plot(np.arange(max_time), simulation_curve[:max_time], alpha=0.25,  color='grey', linewidth=.25)

    fig_ax1.text(0.55, 0.85, color='grey', s='+{:d}, -{:d} % tolerance'.format(int(tolerances[1][-1]), abs(int(tolerances[1][-1]))), transform=fig_ax1.transAxes)
    fig_ax1.fill_between(np.arange(0, tolerances[0][-1]), real_epicurve[:tolerances[0][-1]]*(tolerances[1][-1]/100+1),
                         real_epicurve[:tolerances[0][-1]]*(1-tolerances[1][-1]/100), color='turquoise')
    if asym_tol:
        draw_asymmetric_tolerance_windows(asym_tol, real_epicurve)
    draw_tolerance_windows(tolerances, real_epicurve)
    fig_ax1.set_title('{:d} retained trajectories \n out of {:d} '.format(np.shape(simulation_array)[0], N_simul))
    fig_ax1.set_ylabel('Cases')
    fig_ax1.set_xlabel('Days since outbreak', labelpad=-5)

    colors = ['#1a10ae', '#db0073', '#ff743a', '#edde55', '#1a10ae']
    for index, param in enumerate(param_names):
        fig_ax_param = fig.add_subplot(gs[1, index])
        param_post_distr = np.array(config_array[:, index], dtype=float)
        sns.kdeplot(param_post_distr, shade=True, cut=0, ax=fig_ax_param, color=colors[index])
        sns.rugplot(param_post_distr, ax=fig_ax_param, color='grey')
        fig_ax_param.set_xlabel(param)

def relaxation_method(func, u_0=None):
    '''
    Relaxation algorithm to find roots of function
    Params: func: python functions (most likely the same function with different parameters), 
            u_0: contains initial value
    Returns: a root of func
    '''
    u = u_0
    u_i = func(u_0)
    while (not np.isclose(u, u_i, rtol=0, atol=1e-10)):
        u = u_i
        u_i = func(u)
    return u

def relaxation_for_arrays(func, u_0=None):
    '''
    Relaxation algorithm to find root of function for multiple functions
    Params: func: python functions (most likely the same function with different parameters), 
            u_0(array): contains all initial values
    Returns: a root of func
    '''
    u = u_0
    u_i = func(u_0)
    while (not np.allclose(u.astype(complex), u_i.astype(complex), atol=1e-10)):
        u = u_i
        u_i = func(u)
    return u_i

def partitions(n, k, a):
    """
    Finds and yields all partitions of integer "n" containing "k" elements greater or equal to "a"
    """
    if k == 1 and a <= n :
        yield [n]
    elif n > 0 and k > 0:
        for x in range(a, n):
            for p in partitions(n-x, k-1, x):
                yield [x] + p

def prod(liste):
    '''
    Return the product of all elements in input list
    '''
    prod = 1
    for i in liste:
        prod *= i
    return prod

@njit
def fast_prod(array):
    '''
    Return the product of all elements in input list
    '''
    prod = 1.0
    for i in array:
        prod *= i
    return prod

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]

FACTORIAL_TABLE = np.array([
    1, 1, 2, 6, 24, 120, 720, 5040, 40320,
    362880, 3628800, 39916800, 479001600,
    6227020800, 87178291200, 1307674368000,
    20922789888000, 355687428096000, 6402373705728000,
    121645100408832000, 2432902008176640000], dtype='int64')

@njit
def fast_factorial(n):
    if n > 20:
        return math.gamma(n+1)
    return FACTORIAL_TABLE[n]

@njit
def find_k_fast(array):
    '''
    Returns a liste of permutations for the calculation of analytical SpecialPowExpcut small comp dist
    '''
    m = fast_factorial(len(array))
    permutations = 1
    permutations *= fast_factorial(np.count_nonzero(array==array[0]))
    array_copy = array 
    to_remove = array_copy[0]
    compteur = 0
    for i, value in enumerate(array_copy):
        if value == to_remove:
            array = np.delete(array,i-compteur)
            compteur += 1
    while len(array) >= 1:
        permutations *= fast_factorial(np.count_nonzero(array==array[0]))
        array_copy = array 
        to_remove = array_copy[0]
        compteur = 0
        for i, value in enumerate(array_copy):
            if value == to_remove:
                array = np.delete(array,i-compteur)
                compteur += 1
    return m/permutations

def find_k(liste):
    '''
    Returns the permutation coefficient for the calculation of analytical SpecialPowExpcut small components dist
    '''
    m = factorial(len(liste))
    permutations = 1
    permutations *= factorial(liste.count(liste[0]))
    liste = remove_values_from_list(liste, liste[0])
    while liste:
        permutations *= factorial(liste.count(liste[0]))
        liste = remove_values_from_list(liste, liste[0])
    return m/permutations

@njit
def fast_frac5(array, tau):
    frac5 = 0
    for i in array:
        frac5 += find_k_fast(i)/(fast_prod(i))**tau
    return frac5

@njit
def fast_binom(n, k):
    return fast_factorial(n) // fast_factorial(k) // fast_factorial(n - k)

@njit
def fast_frac4(s, partition, tau):
    frac4 = 0.0
    for m in range(1,s-1): # de 1 à s-2
        frac5 = fast_frac5(partition[m], tau)
        frac4 += frac5 * fast_binom(s, m)
