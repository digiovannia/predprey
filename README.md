# predprey

This program uses matrix population models (MPMs) to compute the expected long-term dynamics of age structure within a predator-prey system. Expanding on the model presented in Jensen and Miller (2001), predprey applies discrete-time updates to the age structures of a predator population and prey population via matrix operations. Given two species whose MPMs are available from the COMADRE database, and a particular parameter a user wishes to vary, the script predprey.py produces (in the specified file format):

1) Prey age structures at each of a specified number of years (.txt)
2) Graphs of the total prey population over time (.png)
3) Graphs of the relative age distribution over time (.png)
4) Number of prey of each age killed by predation each year (.txt)
5) Number of prey of each age killed by density effects/starvation each year (.txt)
6) Graphs of prey lost to predation over time (.png)
7) Graphs of prey lost to density over time (.png)

Alternatively, a user can run the script predprey-single.py, which prints summaries of the information above and displays each plot, corresponding to a simulation based on the default parameters. To compare the effects of variation of several parameters on the following outcomes of interest, the script predprey-heatmap.py produces all of the files above for each parameter, plus a set of heatmaps (.png files) for:

1) Total number of prey that die of predation over the specified number of years
2) Total number of prey that die of density effects
3) Fraction of prey that die of predation
4) Percent of prey in the equilibrium population that are age 0

## Assumptions

1) The rate of predation follows a Type II functional response. Future iterations of this program may benefit from the use of other functional response curves according to the particular predator species studied.

2) Each MPM is either age- or stage-structured. If an MPM has columns that sum to >1 below the first row (that is, the fecundity row), this MPM is assumed to be a two-sex MPM. This assumption may not hold for some stage-structured MPMs, however.

3) Fecundity and survival probabilities for a given population decline with the density of that population according to a logistic function of the form:

f(x) = 1/[1 + p*exp(d*(b-x))-1)]

  Where x is a quantity inverse to the population density. Generally, we assume fecundity and survival follow exponential decline, which approximately holds for certain values of the parameters in the functional form above. However, the logistic function form is used to allow for more flexibility, based on the user's assumptions about the baseline population from which the MPM was derived.

4) The amount of food ("resource") available to the prey population is constant, as in Jensen and Miller (2001). Considering that this quantity may vary significantly with environmental influences and competition with other species of the same trophic level, future versions of this program could be improved with non-constant functions for resource abundance.

5) The strength of density effects on fecundity and survival probabilities is the same for all age/stage classes, as in Jensen and Miller (2001). This simplifying assumption allows for the program to easily adapt to MPMs of different sizes, however, relaxing this assumption may permit more realistic analysis (e.g. plausibly, the oldest or youngest individuals in a population may be more severely excluded from access to food when intraspecific competition is high).

6) Predator fecundity increases with the predation rate, as defined by the type II functional response. This follows the assumption of Rosenzweig and MacArthur (1963).

7) By default, the intrinsic predation risk and functional response do not depend on prey age or predator age. This simplifies the parametrization, but a more realistic iteration of our model would account for the greater vulnerability of different ages of prey and greater hunting ability/fitness of different ages of predators. We include a functionality for customizing the predation matrix if the user desires.

8) Predators will always consume the maximum number of prey that they feasibly can, given the constraints of the functional response equation noted in (6). This allows for the possibility that a predator species over-hunts, leaving the prey population extinct. Plausibly, however, a more accurate model would account for adjustments made by more intelligent predators to ration their limited resources.

9) The prey do not compete with any other species for their food resource, and they have only one predator species. 

## Parameters

num_years:  number of years for which the age structures of the populations are computed

st_prey:  starting prey population, beginning at age/stage 0

st_predator:  starting predator population, beginning at age/stage 0

search:  search rate or area of discovery in functional response; as this parameter increases, the predation rate as a function of prey density more quickly approaches its maximum value

handling:  handling time for predator; the reciprocal of this parameter determines the upper bound for the predation rate

resource:  amount of food available to the prey species; this is proportional to a carrying capacity term in standard logistic growth models, however the equilibrium population density of the prey in the absence of predation may not equal this parameter exactly

max_prey_surv_frac:  factor by which the base survival probabilities for the prey population are multiplied to give the maximum survival probabilities (i.e. in the limit of zero density)

max_predator_surv_frac:  factor by which the base survival probabilities for the predator population are multiplied to give the maximum survival probabilities

max_prey_surv_frac:  factor by which the base fecundity values for the prey population are multiplied to give the maximum fecundity

max_predator_surv_frac:  factor by which the base fecundity values for the predator population are multiplied to give the maximum fecundity

prey_ddf:  strength of the effect of density on prey fecundity; corresponds to a steeper exponential decline of fecundity with increasing density

predator_ddf:  strength of the effect of density on predator fecundity

prey_dds:  strength of the effect of density on prey survival

predator_dds:  strength of the effect of density on predator survival

baseline_prey:  fraction of the maximum population at which the base prey MPM is fixed; e.g. baseline = 1 implies the base MPM represents the population at maximum density, baseline = 0.5 corresponds to a population at half maximum density, etc.

baseline_predator:  baseline fraction for predator MPM

predation_prop:  factor by which functional response is scaled before predator fecundity is multiplied by it; i.e. this parameter controls the strength of the effect of the kill rate of predators on their fecundity

noise_sd:  standard deviation of a mean-zero normal random variable added to the age 0 survival probability of the prey at each time step; if 0, there is no noise

## Guide

The predprey.py script reads from the parameters.txt input file and the base_mpms folder, which must be contained within your working directory. Optionally, it will also read from the predation_matrix.txt file if the custom predation matrix functionality is activated. parameters.txt must be in JSON format, containing key-value pairs such that all keys match those listed under "Parameters" above. This repo contains a sample parameters.txt file. Each base MPM is stored in the base_mpms folder in a text file with the format matching the example provided in this repository, that is, with entries on the same row separated by tabs.

When you execute the command "python predprey.py", you will be prompted to specify:

1) The name of the parameter whose value you would like to vary.
2) A string containing the minimum, maximum (exclusive), and increment over which the parameter value will vary. These three values must
be space-separated.
3) Whether you would like detailed output to be printed to the terminal. If not, only the parameter values and corresponding equilibrium prey age distribution will be printed.
4) The path for your working directory.
5) File name for the prey species.
6) File name for the predator species.
7) Whether you want to supply a custom predation matrix.

Based on 5 and 6, the script will create a folder (if one doesn't already exist) for the desired prey-predator species pair, which will contain subfolders for each of the 8 results specified in the intro paragraph. For each value of the varying parameter within the specified range, text files and graphs for the projection of the prey and predator age structures over time are generated, based on the equations detailed in the equations.tex file in this repo.

Running predprey-single.py (which provides results only for the default parameters) includes only prompts 3-7, and in addition to prompts 3-7, predprey-heatmap.py also includes the following prompts:

1) Number of values taken by each varied parameter.
2) Names of all the parameters the user wishes to vary.
3) For each of the above parameters, minimum and maximum values for the ranges over which they will vary.

## References

1) A. Jensen, Miller D. Age structured matrix predation model for the dynamics of wolf and deer populations. Ecol. Model., 141 (1) (2001), pp. 299-305.
2) M. L. Rosenzweig and R. H. MacArthur, "Graphical Representation and Stability Conditions of Predator-Prey Interactions," The American Naturalist 97, no. 895 (Jul. - Aug., 1963): 209-223.
