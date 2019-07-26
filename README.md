# predprey

This program uses matrix population models (MPMs) to compute the expected long-term dynamics of age structure within a predator-prey system. Given two species whose MPMs are available from the COMADRE database, the script predprey.py:

1) ...

## Assumptions

1) The rate of predation follows a Type II functional response. Future iterations of this program may benefit from the use of other functional response curves according to the particular predator species studied.

2) Each MPM is either age- or stage-structured. If an MPM has columns that sum to >1 below the first row (that is, the fecundity row), this MPM is assumed to be a two-sex MPM. This assumption may not hold for some stage-structured MPMs, however.

3) Fecundity and survival probabilities for a given population decline with the density of that population according to a logistic function of the form:
f(x) = 1/[1 + p*exp(d*(b-x))-1)]
Generally, we assume fecundity and survival follow exponential decline, which holds for certain values of the parameters in the functional form above. However, the logistic function form is used to allow for more flexibility, based on the user's assumptions about the baseline population from which the MPM was derived.

4) The amount of food ("resource") available to the prey population is constant. Considering that this quantity may vary significantly with environmental influences and competition with other species of the same trophic level, future versions of this program could be improved with non-constant functions for resource abundance.

5) The strength of density effects on fecundity and survival probabilities is the same for all age/stage classes. This simplifying assumption allows for the program to easily adapt to MPMs of different sizes, however, relaxing this assumption may permit more realistic analysis (e.g. plausibly, the oldest or youngest individuals in a population may be more severely excluded from access to food when intraspecific competition is high).

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

baseline:  fraction of the maximum population at which the base MPM is fixed; e.g. baseline = 1 implies the base MPM represents the population at maximum density,
baseline = 0.5 corresponds to a population at half maximum density, etc.
