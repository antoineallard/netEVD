

# Epidemiological parameters.
R0    = 1.2
k     = 0.04
sigma = 11
alpha = 2.5

# Simulation parameters
number_of_patient_zero = 1
max_time = 1e3
max_size = 1e6






# Initialize variables holding number of new cases and lists of infected times (total and new)
new_cases = number_of_patient_zero
new_wave_times = [0] * number_of_patient_zero
all_infected_times = new_wave_times

#Loop until we reached either 0 new cases or maximum cases
while 0 < new_cases and len(all_infected_times) < max_size:

    # Number of secondary cases for each newly infected individual.
    secondary = self.EDD_obj.draw_random(new_cases)

    # Initialize lists of transmitters times and transmitted times i.e the times (in days) at which each transmitters were infected relative to
    #   the beggining of the epidemic and the times (in days) at which each new infected were infected relative to their infectors.
    rate = alpha / sigma
    transmitted_times = gamma.rvs(a=alpha, scale=1/rate, size=np.sum(secondary))
    transmitters_times = []
    for j, nb_receiver in enumerate(secondary):
        transmitters_times += [new_wave_times[j]] * nb_receiver
    absolute_transmitted_times = np.array(transmitted_times) + np.array(transmitters_times)

    # Discard times overs maximum time studied.
    new_wave_times = list(absolute_transmitted_times[absolute_transmitted_times < max_time])
    new_cases = len(new_wave_times)

    #Finally add new generation to the total number of cases and start over
    all_infected_times += new_wave_times

if new_cases == 0:
    outbreak_size = len(all_infected_times)


else:
    if first_epicurve:
        first_epicurve = False
        histogram, bin_edges = np.histogram(all_infected_times, bins=np.arange(max_time + 1)) #Histogramme des temps d'infection
        epicurve = list(np.cumsum(histogram))

