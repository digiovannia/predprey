import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import json
import os

#############################################################################
################################# Functions #################################

def pred_rate(search, Nt, handling):
    '''
    Functional response
    '''
    return search*Nt/(1+search*Nt*handling)


def MPM_split(MPM):
    '''
    Given an MPM in the standard form for age structure,
    returns a tuple of the fecundity and survival matrices.

    Handles the case where an MPM includes two sexes as well.
    '''
    if np.any(MPM[1:].sum(axis=0) > 1):
        # Checks if there are two fecundity rows (this may
        # have some false negatives though, because it only
        # checks by making sure that the columns under the
        # fecundity row do not sum to more than 1),
        # in which case they will both be female
        # survival*fecundity.
        M = np.zeros(MPM.shape)
        S = np.zeros(MPM.shape)
        half = MPM.shape[0] // 2
        M[0] = MPM[0]
        M[half] = MPM[half] # second fecundity row
        S[1:half] = MPM[1:half]
        S[half+1:] = MPM[half+1:]
    else:
        M = np.zeros(MPM.shape)
        S = np.zeros(MPM.shape)
        M[0] = MPM[0]
        S[1:] = MPM[1:]
    return M, S

def logistic_scale(DNt, d, b, base, M):
    '''
    Computes the vector of scaling factors that are multiplied
    by base fecundity or survival to give the density-adjusted
    fecundity or survival, based on the assumption of logistically
    declining fecundity and survival probability (in general).
    
    Based on the assumption of exponentially declining fecundity
    and survival probability with density.
    '''
    p = (np.exp(d*(b-DNt))-1)/(np.exp(d*(b-1))-1)
    ratio = np.divide(base - M, M, where=M!=0)
    return np.divide(1, 1 + ratio*p, where=(1 + ratio*p)!=0)

def results(x):
    # Normalized
    all_res = []
    for p in parameters:
        ls = par_vals[p][x]
        new_ls = list(np.array(ls) / np.max(np.array(ls)))
        all_res.append(new_ls)
    return(np.array(all_res))

def title_format(z):
    words = z.split('_')
    cap = [w.capitalize() for w in words]
    title = ' '.join(cap)
    return title

def heatmap(z):
    r = results(z)
    cutoff = np.median(r)
    fig, ax = plt.subplots()
    im = ax.imshow(r)
    x = r.shape[1]
    y = r.shape[0]
    ax.set_xticks(np.arange(x))
    ax.set_xticklabels(['smallest'] + ['']*(num_vals-2) + ['largest'])
    ax.set_yticks(np.arange(y))
    ax.set_yticklabels(parameters)
    ax.set_title(title_format(z))
    for i in range(y):
        for j in range(x):
            number = round(r[i,j], 3)
            if number > cutoff:
                text = ax.text(j, i, str(number),
                               ha="center", va="center")
            else:
                text = ax.text(j, i, str(number),
                               ha="center", va="center", color="w")
    fig.tight_layout()
    return fig

#############################################################################
############################# Reading in data files #########################

verbose = bool(int(input('Print detailed output? (1 for yes, '
                         '0 for no): ')))
source = input(('Enter working directory path. This directory must'
                ' contain the base_mpms folder and parameters.txt'
                ' file: ')) #'C:/Users/antho/Downloads/'
os.chdir(source)
prey_filename = input('Enter prey filename, including .txt: ')
predator_filename = input('Enter predator filename, including .txt: ')
different_predation = bool(int(input('Do you want to supply a custom predation '
                                     'matrix? (1 for yes, 0 for no): ')))
spec_folder = prey_filename[:-4] + '_' + predator_filename[:-4]

# Creates a folder and all the necessary subfolders for the chosen
# prey-predator species pair, if they don't already exist.
if not os.path.exists('./' + spec_folder):
    os.makedirs('./' + spec_folder)
folders = ['/lost_to_dens', '/lost_to_pred', '/density_figures',
           '/dist_figures', '/pred_figures', '/total_figures',
           '/age_structures', '/heatmaps']
for fol in folders:
    if not os.path.exists('./' + spec_folder + fol):
        os.makedirs('./' + spec_folder + fol)
os.chdir('./' + spec_folder)

num_vals = int(input('Enter number of distinct values taken by each '
                     'parameter: '))
parameters = input('Enter names of the parameters whose variation '
                   'you want to compare, separated by spaces: ').split()

par_vals = {}
for p in parameters:
    par_vals[p] = {}

for p in parameters:
    max_min = input('Enter minimum and maximum values for %s, space-'
                    'separated: ' % p).split()
    values = list(np.linspace(float(max_min[0]), float(max_min[1]), num_vals))
    par_vals[p]['values'] = values

prey_MPM = np.array(pd.read_csv('../base_mpms/' + prey_filename,
                                sep='\t', header=None))
intrinsic_growth = np.linalg.eig(prey_MPM)[0].max()
# checking whether this population is expected to grow without density
# dep or pred
predator_MPM = np.array(pd.read_csv('../base_mpms/' + predator_filename,
                                    sep='\t', header=None))
prey_M, prey_S = MPM_split(prey_MPM)
predator_M, predator_S = MPM_split(predator_MPM)
with open('../parameters.txt', 'r') as read_file:
    params = json.load(read_file)
    
pos_integers = ['num_years', 'st_prey', 'st_predator']
fractions = ['baseline_prey', 'baseline_predator']
pos_reals = ['search', 'handling', 'resource',
             'predation_prop']
above_one = ['max_prey_surv_frac', 'max_predator_surv_frac',
             'max_prey_fec_frac', 'max_predator_fec_frac']
    
for par, val in params.items():
    if par in pos_integers:
        assert (type(val) == int and val > 0), 'Invalid %s value' % par
    if par in fractions:
        assert (0 < val < 1), 'Invalid %s value' % par
    if par in pos_reals:
        assert val > 0, 'Invalid %s value' % par
    if par in above_one:
        assert (val >= 1), 'Invalid %s value' % par
    if par == 'noise_sd':
        assert val >= 0, 'Invalid %s value' % par

#############################################################################
########### Defining vectors based on the parameter dictionary ##############

for p in parameters:
    for val in par_vals[p]['values']:
        params[p] = val
        name = 'lost_to_pred/pred_' + p + '_' + str(val).replace('.', ',')
        predation_file = open(name + '.txt', 'w')
        name = 'lost_to_dens/dens_' + p + '_' + str(val).replace('.', ',')
        density_file = open(name + '.txt', 'w')
        name = 'age_structures/age_str_' + p + '_' + str(val).replace('.', ',')
        age_file = open(name + '.txt', 'w')
        predation_data = [[0]*prey_M.shape[0]]
        density_data = [[0]*prey_M.shape[0]]

        q = 2/params['resource'] # per-prey term for predation risk, before
                                # functional response; the more prey resource, which
                                # is assumed proportional to prey area, the lower
                                # the intrinsic predation risk
        if different_predation:
            Q = np.array(pd.read_csv('predation_matrix.txt', sep='\t', header=None))
        else:
            Q = np.ones((prey_M.shape[0], predator_M.shape[0]))
        # Simplifying assumption that predation rates do not differ between ages
        Q = q*Q
        # Simplifying assumption that density dependence strengths do not differ
        # between ages:
        prey_dens_dep_fec = np.array([params['prey_ddf']]*prey_M.shape[0])
        prey_dens_dep_surv = np.array([params['prey_dds']]*prey_M.shape[0])
        predator_dens_dep_fec = np.array([params['predator_ddf']]*predator_M.shape[0])
        predator_dens_dep_surv = np.array([params['predator_dds']]*predator_M.shape[0])

        prey_base_fec = prey_M[0]
        predator_base_fec = predator_M[0]
        prey_base_surv = np.append(np.diagonal(prey_S[1:]), 0)
        predator_base_surv = np.append(np.diagonal(predator_S[1:]), 0)

        # max_prey_surv_frac and max_predator_surv_frac set the limits of
        # the survival probabilities at lowest density, with the restriction
        # that probabilities must be <= 1
        prey_max_surv_factor = np.minimum(np.divide(1, prey_base_surv,
                                                    where=prey_base_surv!=0),
                                        params['max_prey_surv_frac'])
        prey_max_surv = prey_base_surv*prey_max_surv_factor
        predator_max_surv_factor = np.minimum(np.divide(1, predator_base_surv,
                                                where=predator_base_surv!=0),
                                            params['max_predator_surv_frac'])
        predator_max_surv = predator_base_surv*predator_max_surv_factor
        prey_max_fec = prey_base_fec*params['max_prey_fec_frac']
        predator_max_fec = predator_base_fec*params['max_predator_fec_frac']

        preyMPM_list = []

        #############################################################################
        ################### Computing age structures over time ######################

        Qmat = Q.copy()
        prey_ddf = prey_dens_dep_fec.copy()
        prey_dds = prey_dens_dep_surv.copy()
        predator_ddf = predator_dens_dep_fec.copy()
        predator_dds = predator_dens_dep_surv.copy()

        # Initializing prey and predator vectors, assuming populations
        # start at age 0
        prey_structure = np.zeros(prey_M.shape[0])
        prey_structure[0] = params['st_prey']
        predator_structure = np.zeros(predator_M.shape[0])
        predator_structure[0] = params['st_predator']

        age_file.write(' '.join([str(i) for i in list(prey_structure)]) + '\n')

        # Lists for storing results
        total_prey_pop = [params['st_prey']]
        total_predator_pop = [params['st_predator']]
        prey_structure_list = [list(prey_structure)]
        predator_structure_list = [list(predator_structure)] 

        if verbose:
            print('Starting Populations')
            print('Prey:')
            print(prey_structure)
            print('Predators:')
            print(predator_structure)
            print('\n')

        for year in range(params['num_years']):
            Nt = np.sum(prey_structure)
            predNt = np.sum(predator_structure)
            max_predation = pred_rate(params['search'], Nt, params['handling'])
            if Nt == 0: #implies max_predation is 0 as well
                predK = 0 # no resource so starvation probability must be at max
            else:
                predK = Nt/max_predation
            if predK == 0:
                pred_ratio = np.Inf
            else:
                pred_ratio = predNt/predK
            # The following variables are measures of population saturation,
            # to determine density effects.
            prey_DNt = max(min(1, 1-Nt/params['resource']), 0)
            predator_DNt = max(min(1, 1-pred_ratio), 0)

            Sd = prey_S.copy()
            # Adding noise to first-year prey survival
            Sd[1,0] = min(max(0, Sd[1,0]+np.random.normal(0, params['noise_sd'])), 1)

            # Computing density effect matrices for both populations
            prey_dens_fec = np.diag(logistic_scale(prey_DNt, prey_ddf,
                                                    params['baseline_prey'],
                                                    prey_base_fec,
                                                    prey_max_fec))
            prey_dens_surv = np.diag(logistic_scale(prey_DNt, prey_dds,
                                                    params['baseline_prey'],
                                                    prey_base_surv,
                                                    prey_max_surv))
            predator_dens_fec = np.diag(logistic_scale(predator_DNt,
                                                    predator_ddf,
                                                    params['baseline_predator'],
                                                    predator_base_fec,
                                                    predator_max_fec))
            predator_dens_surv = np.diag(logistic_scale(predator_DNt,
                                                        predator_dds,
                                                        params['baseline_predator'],
                                                        predator_base_surv,
                                                        predator_max_surv))                                                                              

            Mw = predator_M.copy()
            Mw = Mw*max_predation*params['predation_prop'] # Predator fecundity
                                                    # increases with kill rate

            # Computing density-modified fecundity and survival matrices
            Md = prey_M.dot(prey_dens_fec)
            Mw = Mw.dot(predator_dens_fec)
            Sd = Sd.dot(prey_dens_surv)
            Sw = predator_S.dot(predator_dens_surv)

            # To determine the predation effect on survival probability, kill rate
            # is multiplied by intrinsic predation risk (f), and the effect on each
            # age class of prey is the sum of the predator population weighted by
            # age-specific predation effect terms in Ffunc.
            Qfunc = max_predation*Qmat # Functional response
            predation_matrix = np.diag(Qfunc.dot(predator_structure))

            # The prey MPM is modified by density and predation effects, and the
            # predator MPM is modified by density effects. We then project the
            # prey and predator age vectors with these modified MPMs, rounding
            # negative numbers up to 0 if necessary.
            prey_MPM_prime = (Md+Sd).dot(np.eye(prey_S.shape[0])-predation_matrix)
            preyMPM_list.append(prey_MPM_prime)
            predator_MPM_prime = Mw + Sw
            new_prey_structure = prey_MPM_prime.dot(prey_structure)
            new_prey_structure = np.maximum(new_prey_structure, 0)
            new_predator_structure = predator_MPM_prime.dot(predator_structure)
            new_predator_structure = np.maximum(new_predator_structure, 0)

            # To calculate the number of prey lost to density effects, we compute
            # the prey vector we would have expected if there were no density
            # effect (fix prey_DNt = 1), and subtract the actual prey vector.
            prey_dens_fec = np.diag(logistic_scale(1, prey_ddf, params['baseline_prey'],
                                                prey_base_fec, prey_max_fec))
            prey_dens_surv = np.diag(logistic_scale(1, prey_dds, params['baseline_prey'],
                                                    prey_base_surv, prey_max_surv))
            Md_sd = prey_M.dot(prey_dens_fec)
            Sd_sd = prey_S.dot(prey_dens_surv) 
            MPMsd = (Md_sd+Sd_sd).dot(np.eye(prey_S.shape[0])-predation_matrix)
            sans_density = MPMsd.dot(prey_structure)
            lost_to_density = sans_density - new_prey_structure

            # To calculate the number of prey lost to predation, we compute the
            # prey vector expected assuming no predation effect, and subtract
            # the actual prey vector.
            sans_predation = (Md + Sd).dot(prey_structure)
            lost_to_pred = sans_predation - new_prey_structure

            # Update the vectors
            prey_structure = new_prey_structure
            predator_structure = new_predator_structure

            total_prey_pop.append(np.sum(prey_structure))
            total_predator_pop.append(np.sum(predator_structure))
            prey_structure_list.append(list(prey_structure))
            predator_structure_list.append(list(predator_structure))

            predation_data.append(list(lost_to_pred))
            density_data.append(list(lost_to_density))
            predation_file.write(' '.join([str(i) for i in list(lost_to_pred)]) + '\n')
            density_file.write(' '.join([str(i) for i in list(lost_to_density)]) + '\n')
            age_file.write(' '.join([str(i) for i in list(prey_structure)]) + '\n')

            if verbose:
                print('Year', str(year+1))
                print('Prey Age Structure:')
                print(list(prey_structure))
                print('Prey Lost to Predation:')
                print(list(lost_to_pred))
                print('Prey Lost to Density:')
                print(list(lost_to_density))
                print('Modified Prey MPM:')
                print(prey_MPM_prime)
                print('Predators:')
                print(list(predator_structure))
                print('\n')

        fig = plt.figure()
        plt.plot(total_prey_pop)
        plt.title('Total Prey Population Over Time')
        name = 'total_figures/prey_pop_' + p + '_' + str(val).replace('.', ',')
        fig.savefig(name + '.png')
        plt.cla()

        # Age distributions
        with np.errstate(divide='ignore', invalid='ignore'):
            x = range(params['num_years']+1)
            y = np.array(prey_structure_list).T / np.array(total_prey_pop)
            plt.stackplot(x, y)
        plt.title('Prey Age Structure Over Time')
        name = 'dist_figures/prey_rel_age_' + p + '_' + str(val).replace('.', ',')
        fig.savefig(name + '.png')
        plt.cla()

        y = np.array(predation_data).T
        plt.stackplot(x, y)
        plt.title('Prey Killed by Predation')
        name = 'pred_figures/by_predation_' + p + '_' + str(val).replace('.', ',')
        fig.savefig(name + '.png')
        plt.cla()

        y = np.array(density_data).T
        plt.stackplot(x, y)
        plt.title('Prey Killed by Density')
        name = 'density_figures/by_density_' + p + '_' + str(val).replace('.', ',')
        fig.savefig(name + '.png')
        plt.cla()

        print(p + ' == ' + str(val).replace('.', ','))
        print('\n')
        print('Prey Equilibrium Age Structure(s):')
        if total_prey_pop[-1] > 0:
            prey_eq_age = np.array(prey_structure_list[-1])/total_prey_pop[-1]
            prey_eq_last_age = np.array(prey_structure_list[-2])/total_prey_pop[-2]
            if not np.allclose(prey_eq_age, prey_eq_last_age):
                # There may be fluctuations in the age distribution. If the last
                # two age distributions differ substantially, check if every other
                # age distribution is equal. If so, the "equilibrium" is a pair
                # of distributions.
                prey_step = np.array(prey_structure_list[-3])/total_prey_pop[-3]
                prey_ls_step = np.array(prey_structure_list[-4])/total_prey_pop[-4]
                if np.allclose(prey_eq_age, prey_step):
                    if np.allclose(prey_eq_last_age, prey_ls_step):
                        prey_eq_age = [prey_eq_age, prey_eq_last_age]
                else:
                    prey_eq_age = 'No Equilibrium'
        else:
            prey_eq_age = 'Extinct'
        print(prey_eq_age)

        predation_file.close()
        density_file.close()
        age_file.close()
        plt.close(fig)

    with open('../parameters.txt', 'r') as read_file:
        params = json.load(read_file)


for p in parameters:
    x = 'age_structures/'
    percent_age_0_at_equilibrium = []
    for v in par_vals[p]['values']:
        V = str(v).replace('.', ',')
        data = pd.read_csv(x + 'age_str_' + p + '_' + V + '.txt',
                           header=None, sep=' ')
        end_row = data.iloc[-1,]
        eq_dist = end_row / (end_row.sum())
        percent_age_0_at_equilibrium.append(eq_dist[0])
    par_vals[p]['percent_age_0_at_equilibrium'] = percent_age_0_at_equilibrium

for p in parameters:
    x = 'lost_to_dens/'
    total_deaths_by_density = []
    for v in par_vals[p]['values']:
        V = str(v).replace('.', ',')
        data = pd.read_csv(x + 'dens_' + p + '_' + V + '.txt',
                           header=None, sep=' ')
        total_died = data.values.sum()
        total_deaths_by_density.append(total_died)
    par_vals[p]['total_deaths_by_density'] = total_deaths_by_density

for p in parameters:
    x = 'lost_to_pred/'
    total_deaths_by_predation = []
    for v in par_vals[p]['values']:
        V = str(v).replace('.', ',')
        data = pd.read_csv(x + 'pred_' + p + '_' + V + '.txt',
                           header=None, sep=' ')
        total_died = data.values.sum()
        total_deaths_by_predation.append(total_died)
    par_vals[p]['total_deaths_by_predation'] = total_deaths_by_predation

for p in parameters:
    x = 'lost_to_pred/'
    y = 'lost_to_dens/'
    fraction_deaths_by_predation = []
    for v in par_vals[p]['values']:
        V = str(v).replace('.', ',')
        data_x = pd.read_csv(x + 'pred_' + p + '_' + V + '.txt',
                             header=None, sep=' ')
        data_y = pd.read_csv(y + 'dens_' + p + '_' + V + '.txt',
                             header=None, sep=' ')
        total_pred = data_x.values.sum()
        total_dens = data_y.values.sum()
        z = total_pred / (total_pred + total_dens)
        fraction_deaths_by_predation.append(z)
    par_vals[p]['fraction_deaths_by_predation'] = fraction_deaths_by_predation

h1 = heatmap('fraction_deaths_by_predation')
h1.savefig('heatmaps/fraction_deaths_by_predation_' + '_'.join(parameters))
h2 = heatmap('total_deaths_by_predation')
h2.savefig('heatmaps/total_deaths_by_predation_' + '_'.join(parameters))
h3 = heatmap('total_deaths_by_density')
h3.savefig('heatmaps/total_deaths_by_density_' + '_'.join(parameters))
h4 = heatmap('percent_age_0_at_equilibrium')
h4.savefig('heatmaps/percent_age_0_at_equilibrium_' + '_'.join(parameters))
