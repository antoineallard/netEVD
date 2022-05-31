'''
File containing multiple class containing all the methods for each excess degree distributions
'''

import numpy as np
import scipy as sp
import scipy.special as spec
import mpmath as mpm
import pickle
import os

#Local imports
from . import other_functions as fct



ParameterError = ValueError("Parameter out of range.")

class NBinom(object):
    """ Negative-binomial model."""
    def __init__(self, custom_p0=False):
        # init
        self.params = [r'$R_0$', r'$k$', r'$\sigma$', r'$\gamma_{shape}$']
        self.label = [r'$\tilde{R}_0$=', r'$\tilde{k}$=']
        self.params_values = [None,None,None,None]
        self.num_param = 4  # number of parameters
        self.k = None
        self.R0 = None
        self.custom_p0 = custom_p0
        self.p0 = None
        if custom_p0 is not None:
            self.p0 = custom_p0
        self.sigma = None
        self.gamma_shape = None

    def set_params(self, params, simulation=True):
        self.set_param_A(params[0])
        self.set_param_B(params[1])
        self.set_param_C(params[2])
        if not self.custom_p0:
            self.set_p0()
        if simulation:
            self.set_sigma(params[3])
            self.set_gamma_shape(params[4])

    def set_p0(self):
        if np.isclose(self.k, 1, atol=1e-8):
            self.p0 = 1 - np.log(1 + self.R0)
        else:
            self.p0 = 1 - (self.k * ((self.k / (self.k + self.R0)) ** (self.k - 1) - 1) / (1 - self.k))
        if self.p0 < 0:
            self.p0 = 0

    def set_param_A(self, R0):
        if R0 <= 0:
            raise(ParameterError)
        self.params_values[0] = R0
        self.R0 = R0

    def set_param_B(self, k):
        if k <= 0:
            raise(ParameterError)
        self.k = k
        self.params_values[1] = k

    def set_sigma(self, sig):
        self.sigma = sig
        self.params_values[2] = sig

    def set_gamma_shape(self, gamma):
        self.gamma_shape = gamma
        self.params_values[3] = gamma

    def set_param_C(self, alph):
        pass

    def G1(self, x, fft=True):
        return (1 + (self.R0 / self.k) * (1 - x)) ** (-self.k)

    def G0(self, x, fft=True):
        if np.isclose(self.k, 1, atol=1e-10):
            G0 = self.p0 + (1 - self.p0) *(1 - np.log(1 + self.R0 * (1 - x)) / np.log(1 + self.R0))
        else:
            G0 = self.p0 + (1 - self.p0) * ((1 - (1 - self.R0 * x / (self.R0 + self.k))**(1 -self.k)) / (1 - ( self.k / (self.R0 + self.k)) ** (1 - self.k)))
        return G0

    def draw_random(self, size):
        '''
        Draws time of infection for next wave of size wave_size from self.distribution
        '''
        prob_p = self.k / (self.R0 + self.k)
        drawn = sp.stats.nbinom.rvs(n=self.k, p=prob_p, size=size)
        return drawn

    def small_comp_dist(self, s_max):
        '''
        Computes analytically the small components distribution
        '''
        if self.k * s_max <100:
            distribution = [self.p0]
            for s in range(2,s_max):
                frac_1 = (1 - self.p0)
                frac_2 = (1 - self.k) / ((self.k / (self.R0 + self.k)) ** (self.k - 1) - 1)
                frac_3 = np.exp ( sp.special.loggamma(s * self.k + s - 2) - sp.special.loggamma(s * self.k) - sp.special.loggamma(s) )
                # frac_4 = (self.R0 / self.k) ** (s/2 )
                # frac_45 = (self.R0 / self.k) ** (s/2 - 1)
                # frac_5 = (1 + self.R0 / self.k) ** (-s * self.k - s + 2)
                # distribution.append(frac_1 * frac_2 * frac_3 * frac_4 * frac_5 * frac_45)
                frac_4 = np.exp( (s - 1) * np.log(self.R0 / self.k) - (s * self.k + s - 2) * np.log(1 + self.R0 / self.k) )
                distribution.append(frac_1 * frac_2 * frac_3 * frac_4)
        else:
            distribution = get_component_size_dist(self)[1:s_max]
        return distribution

    def calculate_final_size(self):
        '''
        Computes final size projection with standard network models
        '''
        u = fct.relaxation_method(self.G1, 0.634)
        R_infinity = 1 - self.G0(u)
        return R_infinity

#===================================================================//=============================================================================#

__location__ = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__))))

with open(os.path.join(__location__, "partitions.txt"), "rb") as fp:
    PARTITION_TABLE = pickle.load(fp)

class ShiftedPowExpcut(object):
    """Shifted PowExpcut model."""
    def __init__(self, custom_p0=False):

        self.params = [r'$\kappa$', r'$\tau$', r'$\sigma$', r'$\gamma_{shape}$']
        self.label = [r'$\tilde{\kappa}$=', r'$\tilde{\tau}$=']
        self.params_values = [None,None,None,None]
        self.num_param = 4  # number of parameters
        self.kappa = None
        self.tau = None
        self.custom_p0 = custom_p0
        self.p0 = None
        if custom_p0 is not None:
            self.p0 = custom_p0
        self.sigma = None
        self.gamma_shape = None

    def set_params(self, params, simulation=True):
        self.set_param_A(params[0])
        self.set_param_B(params[1])
        self.set_param_C(params[2])
        if simulation:
            self.set_sigma(params[3])
            self.set_gamma_shape(params[4])
        if not self.custom_p0:
            self.set_p0()

    def set_p0(self):
        Li = np.frompyfunc(lambda degree, x: mpm.fp.polylog(degree, x), 2, 1)
        self.p0 = 1 - (Li(self.tau - 1,np.exp(-1/self.kappa))
                    - Li(self.tau, np.exp(-1/self.kappa))) * Li(self.tau+1, np.exp(-1/self.kappa)) / Li(self.tau, np.exp(-1/self.kappa))**2
        if self.p0 < 0:
            self.p0 = 0
        if self.p0 > 1:
            self.p0 = 1

    def set_param_A(self, kappa):
        if kappa <= 0:
            raise(ParameterError)
        self.params_values[0] = kappa
        self.kappa = kappa

    def set_param_B(self, tau):
        if tau <= 0:
            raise(ParameterError)
        self.tau = tau
        self.params_values[1] = tau

    def set_sigma(self, sig):
        self.sigma = sig
        self.params_values[2] = sig

    def set_gamma_shape(self, gamma):
        self.gamma_shape = gamma
        self.params_values[3] = gamma

    def set_param_C(self, alph):
        pass


    def G1(self, x, fft=False):
        polylog_tau = np.frompyfunc(lambda z: mpm.fp.polylog(self.tau, z), 1, 1)
        if fft:
            return polylog_tau(x * np.exp(-1 / self.kappa)) / (float(mpm.fp.polylog(self.tau, np.exp(-1 / self.kappa))) * x)
        return (1/x) * (float(mpm.fp.polylog(self.tau, x*np.exp(-1 / self.kappa)))) * (1/float(mpm.fp.polylog(self.tau, np.exp(-1 / self.kappa))))

    def G0(self, x, fft=False):
        polylog_tau_p1 = np.frompyfunc(lambda z: mpm.fp.polylog(self.tau + 1, z), 1, 1)
        if fft:
            return self.p0 + (1 - self.p0) * ((polylog_tau_p1(x * np.exp(-1 / self.kappa)))) * (1/float(polylog_tau_p1(np.exp(-1 / self.kappa))))
        return self.p0 + (1 - self.p0) * (float(mpm.fp.polylog(self.tau + 1, x * np.exp(-1 / self.kappa)))) * (1/float(mpm.fp.polylog(self.tau + 1, np.exp(-1 / self.kappa))))


    def draw_random(self, size):
        '''
        Draws time of infection for next wave of size wave_size from self.distribution
        '''
        drawn = np.zeros(size)
        enough = np.array([False for i in range(size)])
        while not np.all(enough):
            r = np.random.uniform(size=size)
            k = (-self.kappa*np.log(1-r)*np.array([not i for i in enough])).astype(int)
            sel_vect = np.random.uniform(size=size) < k**self.tau
            drawn[sel_vect]= k[sel_vect]
            enough += sel_vect
        drawn = np.array([int(round(i)) for i in drawn])
        return drawn

    def small_comp_dist(self, s_max):
        '''
        Computes analytically the small components distribution
        '''
        if s_max > 45:
            eta_s = get_component_size_dist(self, 4096) # Numerical version (faster for s>45)
        if s_max < 46:
            eta_s = self.compute_analytical_distribution(s_max) #Analytical version

        return eta_s[1:s_max]


    def calculate_final_size(self):
        '''
        Computes final size projection with standard network models
        '''
        u = fct.relaxation_method(self.G1, 0.634)
        R_infinity = 1 - self.G0(u)
        return R_infinity



    def compute_analytical_distribution(self, s_max):
        '''
        Compute small composants size distribution analytically
        Params: s_max: maximum comp size evaluated (int)
        Returns: Array containing distribution of probabilities of sizes of small-sized comp
        '''
        Li = np.frompyfunc(lambda degree, x: mpm.fp.polylog(degree, x), 2, 1)
        distribution = [0,self.p0]
        for s in range(2,s_max):
            frac1 = (1-self.p0)/Li(self.tau + 1,np.exp(-1/self.kappa))
            frac2 = 1/(Li(self.tau,np.exp(-1/self.kappa)))**(s-1)
            frac3 = np.exp( (2-2*s) / self.kappa)
            if s==2:
                distribution.append(frac1*frac2*frac3)
            else:
                frac4 = 0
                for m in range(1,s-1): # de 1 à s-2
                    partition = np.array(PARTITION_TABLE[s,m])
                    frac5 = fct.fast_frac5(partition, self.tau)
                    frac4 += frac5 * spec.binom(s, m)
                prob = frac1*frac2*frac3*frac4/(s-1)
                distribution.append(prob)
        return distribution



#=================================================================================//==============================================================================================#

class Poisson(object):
    """Poisson distribution model."""
    def __init__(self, custom_p0=False):
        # init
        self.params = [r'$z$', r'$\sigma$', r'$\gamma_{shape}$']
        self.label = [r'$\tilde{z}$=']
        self.params_values = [None,None,None]
        self.num_param = 3  # number of parameters
        self.z = None
        self.p0 = None
        self.custom_p0 = custom_p0
        if custom_p0 is not None:
            self.p0 = custom_p0
        self.sigma = None
        self.gamma_shape = None

    def set_params(self, params, simulation=True):
        self.set_param_A(params[0])
        self.set_param_B(params[1])
        self.set_param_C(params[2])
        if not self.custom_p0:
            self.set_p0()
        if simulation:
            self.set_sigma(params[3])
            self.set_gamma_shape(params[4])

    def set_p0(self):
        self.p0 = np.exp(-self.z)

    def set_param_A(self, z):
        if z <= 0:
            raise(ParameterError)
        self.params_values[0] = z
        self.z = z

    def set_param_B(self, k):
        pass

    def set_sigma(self, sig):
        self.sigma = sig
        self.params_values[1] = sig

    def set_gamma_shape(self, gamma):
        self.gamma_shape = gamma
        self.params_values[2] = gamma

    def set_param_C(self, alph):
        pass


    def G1(self, x):
        return np.exp(self.z * (x - 1))

    def G0(self, x):
        return np.exp(self.z * (x - 1))

    def draw_random(self, size):
        '''
        Draws time of infection for next wave of size wave_size from self.distribution
        '''
        drawn = sp.stats.poisson.rvs(mu=self.z, size=size)
        return drawn

    def small_comp_dist(self, s_max):
        '''
        Computes analytically the small components distribution
        '''
        eta_s = [0]
        for s in range(1,s_max):
            eta = (self.z*s)**(s-1)*np.exp(-self.z*s)/fct.fast_factorial(s)
            eta_s.append(eta)
        return eta_s[1:s_max]

    def calculate_final_size(self):
        '''
        Computes final size projection with standard network models
        '''
        u = fct.relaxation_method(self.G1, 0.634)
        R_infinity = 1 - self.G0(u)
        return R_infinity

#===================================================================================//=============================================================================#

class Exponential(object):
    """Geometric distribution model with p=exp(-1/kappa)."""
    def __init__(self, custom_p0=False):
        # init
        self.params = [r'$\kappa$', r'$\sigma$', r'$\gamma_{shape}$']
        self.label = [r'$\tilde{\kappa}$=']
        self.params_values = [None,None,None]
        self.num_param = 3  # number of parameters
        self.kappa = None
        self.p0 = None
        self.custom_p0 = custom_p0
        if custom_p0 is not None:
            self.p0 = custom_p0
        self.sigma = None
        self.gamma_shape = None

    def set_params(self, params, simulation=True):
        self.set_param_A(params[0])
        self.set_param_B(params[1])
        self.set_param_C(params[2])
        if not self.custom_p0:
            self.set_p0()
        if simulation:
            self.set_sigma(params[3])
            self.set_gamma_shape(params[4])

    def set_p0(self):
        self.p0 = np.log(np.exp(1/self.kappa)-1)-(1/self.kappa)+1
        if self.p0 < 0:
            self.p0 = 0

    def set_param_A(self, kappa):
        if kappa <= 0:
            raise(ParameterError)
        self.params_values[0] = kappa
        self.kappa = kappa

    def set_param_B(self, k):
        pass

    def set_sigma(self, sig):
        self.sigma = sig
        self.params_values[1] = sig

    def set_gamma_shape(self, gamma):
        self.gamma_shape = gamma
        self.params_values[2] = gamma

    def set_param_C(self, alph):
        pass


    def G1(self, x):
        return (1-np.exp(-1/self.kappa))/(1-x*np.exp(-1/self.kappa))

    def G0(self, x):
        return self.p0 - (np.log( abs( x - np.exp( 1/self.kappa ) ) ) + (1/self.kappa) )

    def draw_random(self, size):
        '''
        Draws time of infection for next wave of size wave_size from self.distribution
        '''
        p = 1-np.exp(-1/self.kappa)
        drawn = sp.stats.geom.rvs(p=p, loc=-1, size=size)
        return drawn

    def small_comp_dist(self, s_max, resolution=4096):
        '''
        Computes analytically the small components distribution

        '''
        #eta_2 = pgf.get_component_size_dist( mo.Exponential(self.kappa,self.p0), 4096) #J-G. version
        eta_a = [0,self.p0]
        for s in range(2,s_max):
            frac1 = (1-self.p0)/(1-np.exp(1/self.kappa))
            frac1p5 = 1/(np.log(np.exp(1/self.kappa)-1)-1/self.kappa)
            frac2 = (1-np.exp(-1/self.kappa))**s / fct.fast_factorial(s-1)
            frac3 = spec.gamma(2*s-2)/spec.gamma(s)
            frac4 = np.exp((2-s)/self.kappa)
            eta_a.append(frac1*frac1p5*frac2*frac3*frac4)
        return eta_a[1:s_max]

    def calculate_final_size(self):
        '''
        Computes final size projection with standard network models
        '''
        u = fct.relaxation_method(self.G1, 0.634)
        R_infinity = 1 - self.G0(u)
        return R_infinity

#===================================================================//=============================================================================#

class Binom(object):
    """ Binomial model."""
    def __init__(self, custom_p0=False):
        # init
        self.params = [r'$R_0$', r'$k$', r'$\sigma$', r'$\gamma_{shape}$']
        self.label = [r'$\tilde{R}_0$=', r'$\tilde{k}$=']
        self.params_values = [None,None,None,None]
        self.num_param = 4  # number of parameters
        self.k = None
        self.R0 = None
        self.custom_p0 = custom_p0
        self.p0 = None
        if custom_p0 is not None:
            self.p0 = custom_p0
        self.sigma = None
        self.gamma_shape = None

    def set_params(self, params, simulation=True):
        self.set_param_A(params[0])
        self.set_param_B(params[1])
        self.set_param_C(params[2])
        if not self.custom_p0:
            self.set_p0()
        if simulation:
            self.set_sigma(params[3])
            self.set_gamma_shape(params[4])

    def set_p0(self):
        self.p0 = 1-(1-(1-self.R0/self.k)**(self.k+1))*(self.k/(self.k+1))
        if self.p0 < 0:
            self.p0 = 0

    def set_param_A(self, R0):
        if R0 <= 0:
            raise(ParameterError)
        self.params_values[0] = R0
        self.R0 = R0

    def set_param_B(self, k):
        if k <= 0:
            raise(ParameterError)
        self.k = k
        if self.R0/self.k > 1:
            raise(ParameterError)
        self.params_values[1] = k

    def set_sigma(self, sig):
        self.sigma = sig
        self.params_values[2] = sig

    def set_gamma_shape(self, gamma):
        self.gamma_shape = gamma
        self.params_values[3] = gamma

    def set_param_C(self, alph):
        pass

    def G1(self, x , fft=True):
        return (1 + (self.R0 / self.k) * (x - 1)) ** (self.k)

    def G0(self, x, fft=True):
        return self.p0 + self.k*((1+self.R0*(x-1)/self.k)**(self.k+1)-(1-self.R0/self.k)**(self.k+1))/(self.k+1)

    def draw_random(self, size):
        '''
        Draws time of infection for next wave of size wave_size from self.distribution
        '''
        prob_p = self.R0/self.k
        n = int(self.k)
        drawn = sp.stats.binom.rvs(n=n, p=prob_p, size=size)
        return drawn

    def small_comp_dist(self, s_max):
        '''
        Computes analytically the small components distribution
        '''
        #eta_s = pgf.get_component_size_dist( mo.NegativeBinomial(self.R0, self.k,self.p0), 4096) #J-G version
        if s_max*self.k < 150:
            distribution = [self.p0]
            for s in range(2,s_max):
                frac_1 = self.R0/fct.fast_factorial(s-1)
                frac_2 = (self.R0/self.k)**(s-2)
                frac_3 = sp.special.gamma((s * self.k) +1)/sp.special.gamma((s * self.k) - s + 3)
                frac_4 = (1-self.R0/self.k) ** (self.k*s-s+2)
                distribution.append(frac_1 * frac_2 * frac_3 * frac_4)
        else:
            distribution = get_component_size_dist(self)[1:s_max]
        return distribution

    def calculate_final_size(self):
        '''
        Computes final size projection with standard network models
        '''
        u = fct.relaxation_method(self.G1, 0.634)
        R_infinity = 1 - self.G0(u)
        return R_infinity

#===================================================================//=============================================================================#

class ShiftedPowerLaw(object):
    """Shifted PowExpcut model."""
    def __init__(self, custom_p0=False):

        self.params = [r'$\tau$', r'$\sigma$', r'$\gamma_{shape}$']
        self.label = [r'$\tilde{\tau}$=']
        self.params_values = [None,None,None]
        self.num_param = 3  # number of parameters
        self.tau = None
        self.p0 = None
        self.custom_p0 = custom_p0
        if custom_p0 is not None:
            self.p0 = custom_p0
        self.sigma = None
        self.gamma_shape = None

    def set_params(self, params, simulation=True):
        self.set_param_A(params[0])
        self.set_param_B(params[1])
        self.set_param_C(params[2])
        if simulation:
            self.set_sigma(params[3])
            self.set_gamma_shape(params[4])
        if not self.custom_p0:
            self.set_p0()

    def set_p0(self):
        self.p0 = 1 + spec.zeta(self.tau+1)*(spec.zeta(self.tau)-spec.zeta(self.tau-1))*(1/(spec.zeta(self.tau))**2)
        if self.p0 < 0:
            self.p0 = 0
        if self.p0 > 1:
            self.p0 = 1

    def set_param_A(self, tau):
        if tau < 2:
            raise(ParameterError)
        self.params_values[0] = tau
        self.tau = float(tau)

    def set_param_B(self, tau):
        pass

    def set_sigma(self, sig):
        self.sigma = sig
        self.params_values[1] = sig

    def set_gamma_shape(self, gamma):
        self.gamma_shape = gamma
        self.params_values[2] = gamma

    def set_param_C(self, alph):
        pass


    def G1(self, x, fft=False):
        polylog_tau = np.frompyfunc(lambda z: mpm.fp.polylog(self.tau, z), 1, 1)
        return polylog_tau(x) / (spec.zeta(self.tau) * x)


    def G0(self, x, fft=False):
        polylog_tau_p1 = np.frompyfunc(lambda z: mpm.fp.polylog(self.tau + 1, z), 1, 1)
        return self.p0 + (1-self.p0)*(polylog_tau_p1(x))/(spec.zeta(self.tau+1))


    def draw_random(self, size):
        '''
        Draws time of infection for next wave of size wave_size from self.distribution
        '''
        drawn = sp.stats.zipf.rvs(a=self.tau, loc=-1, size=size)
        return drawn

    def small_comp_dist(self, s_max):
        '''
        Computes analytically the small components distribution
        '''
        if s_max > 45:
            eta_s = get_component_size_dist(self, 4096) # Numerical version (faster for s>45)
        if s_max < 46:
            eta_s = self.compute_analytical_distribution(s_max) #Analytical version
        return eta_s[1:s_max]


    def calculate_final_size(self):
        '''
        Computes final size projection with standard network models
        '''
        u = fct.relaxation_method(self.G1, 0.634)
        R_infinity = 1 - self.G0(u)
        return R_infinity



    def compute_analytical_distribution(self, s_max):
        '''
        Compute small composants size distribution analytically
        Params: s_max: maximum comp size evaluated (int)
        Returns: Array containing distribution of probabilities of sizes of small-sized comp
        '''
        distribution = [0,self.p0]
        for s in range(2,s_max):
            frac1 = (spec.zeta(self.tau-1)-spec.zeta(self.tau))/(spec.zeta(self.tau))**(s+1)
            frac2 = 1/(s-1)
            if s==2:
                distribution.append(frac1*frac2)
            else:
                frac4 = 0
                for m in range(1,s-1): # de 1 à s-2
                    partition = np.array(PARTITION_TABLE[s,m])
                    frac5 = fct.fast_frac5(partition, self.tau)
                    frac4 += frac5 * spec.binom(s, m)
                prob = frac1*frac2*frac4
                distribution.append(prob)
        return distribution


# ============================================================================ #
#                             Common functions                                 #
# ============================================================================ #



def get_component_size_dist(model, resolution=4096):
    """Solve for the distribution of component sizes using FFTs."""

    # Sample of the unit circle
    z = np.exp(2j * np.pi  * np.arange(resolution) / resolution)

    # Compute the (evaluated) Kn PGF.
    A_n = get_An(z, model)
    K_n = get_Kn(z, A_n, model)

    # Extract coefficient with the inverse FFT.
    distribution = np.roll(np.absolute(np.fft.ifft(K_n))[::-1], 1)
    return distribution

def get_An(x, model):
    """Compute the (evaluated) An probability generating function.

    Solves for A_n by iteration of An(x) = x * G1(An(x)).
    """
    try:
        iter(x)
    except TypeError:
        x = np.array([x])
    An = 0.001 * np.ones(len(x))
    new_An = np.zeros(len(x))
    while True:
        new_An = x * model.G1(An, fft=True)
        if not np.allclose(An.astype(complex), new_An.astype(complex), atol=1e-8):
            An = new_An.copy()
        else:
            break
    return new_An

def get_Kn(x, An, model):
    """Compute the (evaluated) K_n probability generating function.

    Solves directly using Kn = x * G0(An(x)).
    """
    return x * model.G0(An, fft=True)
