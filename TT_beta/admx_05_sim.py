import msprime as msp
import numpy as np
import TTo_func
from time import time

# Population 1, 2 and the outgroup
P1, P2, O = 0, 1, 2

def set_up_pops(nS, tS, n0):
    # Samples from population 1
    samples = [msp.Sample(population=P1, time=tS[0])]*(1*nS[0])
    # Samples from population 2
    samples.extend([msp.Sample(population=P2, time=tS[0])]*(1*nS[0]))
    # Samples from outgroup
    samples.extend([msp.Sample(population=O, time=tS[0])]*(1*n0))
    return samples

def set_up_demography(t12, t120, tBT, N_bn, N):
    # determine topology of tree (divergences)
    divergence = [msp.MassMigration(time=t12,source=P2,destination=P1,proportion=1),
                  msp.MassMigration(time=t120,source=O,destination=P1,proportion=1)]
    # admixture times and proportion
    #admix = [msp.MassMigration(time=ta,source=P1,destination=P2,proportion=f)]
    admix2 = [msp.MigrationRateChange(time=0.95*t12, rate=0.5, matrix_index=(P1,P2)),
              msp.MigrationRateChange(time=t12, rate=0, matrix_index=(P1,P2))]
    # sum up and sort the demography
    demography = divergence + admix2
    #Version without admixture
    #demography = divergence

    return sorted(demography, key = lambda x: x.time)


# Model parameters
L = 10000000 #length of sequence
Ne = 10000.0 #effective population size
mu = 2.5e-8 #mutation rate
r = 1.25e-8 #recombination rate
generation_time = 25
t12 = 50000/generation_time #50000 years ago in generations
t120 = 100000/generation_time #80000 years ago in generations
tBT = [0.45*t12, 0.55*t12] # start and end date of bottleneck
N_bn = 2000 # population size during bottleneck
#ta = 45000/generation_time
nS=[2] # number of samples per population
n0=1 # number of samples from the outgroup
tS=[0] #all samples taken in modern times?
#f=0.05
N=[Ne]*3 #effective population sizes of the three pops
seed=None
migration_matrix = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]] # original migration matrix

n_reps = 10000
start_time = time()

samples = set_up_pops(nS,tS,n0)
demography = set_up_demography(t12, t120, tBT, N_bn, Ne)
pops = [msp.PopulationConfiguration(initial_size = n) for n in N]

# dd = msp.DemographyDebugger(population_configurations=pops,
                            # demographic_events=demography,
                            # migration_matrix=migration_matrix)
# dd.print_history()

sims = msp.simulate(samples=samples,Ne=N[0],population_configurations=pops,
                    demographic_events = demography, mutation_rate = mu,
                    length=L, recombination_rate = r, random_seed = seed,
                    num_replicates=n_reps)
                    #migration_matrix=migration_matrix)

# input: variant.genotypes from TreeSequence object
def count_m(geno):
    pop1 = geno[0]+geno[1]
    pop2 = geno[2]+geno[3]

    if pop1+pop2==1:
        if pop1==1:
            return 0
        else:
            return 1
    if pop1+pop2==2:
        if pop1==2:
            return 2
        elif pop2==2:
            return 3
        else:
            return 4
    if pop1+pop2==3:
        if pop1==2:
            return 5
        else:
            return 6
    if pop1+pop2==0 or pop1+pop2==4:
        return 7
    print("ERROR: Something went wrong!")
    return(1000)

norm_arr = np.empty((n_reps,4))
bound_arr = np.empty((n_reps,4))
conv_arr = np.empty((n_reps,4))
TTo_arr = np.empty((n_reps,2))
TT_arr = np.empty((n_reps,2))


def convert_to_years(input):

    if input[0] != 'NaN':
        res = g*np.array(input)/mu
    else:
        res = np.empty(len(input))
        res[:] = np.NaN

    return res


# run simulations
for i, sim in enumerate(sims):
    # sims = msp.simulate(samples=samples,Ne=N[0],population_configurations=pops,
                        # demographic_events = demography, mutation_rate = mu,
                        # length=L, recombination_rate = r, random_seed = seed,
                        # migration_matrix=migration_matrix)

    m = [0,0,0,0,0,0,0,0]
    m_cond = [0,0,0,0,0,0,0,0]

    for variant in sim.variants():
        geno = variant.genotypes
        # conditional m
        if geno[4]==1:
            index = count_m(geno)
            m_cond[index] = m_cond[index]+1
        index = count_m(geno)
        m[index] = m[index]+1

    m[7] = L - sum(m[0:7])

    results = TTo_func.estimate_param(m,m_cond)
    results_b = TTo_func.estimate_param(m,m_cond,type='bound')
    results_c = TTo_func.estimate_param(m,m_cond,type='conv')
    results_classic = TTo_func.estimate_param_classic(m,m_cond)
    results_TT = TTo_func.estimate_param_TT(m)
    g = generation_time

    res = convert_to_years(results[8:12])
    res_b = convert_to_years(results_b[8:12])
    res_c = convert_to_years(results_c[8:12])
    res_classic = convert_to_years(results_classic[16:18])
    res_TT = convert_to_years(results_TT[3:5])

    norm_arr[i] = res
    bound_arr[i] = res_b
    conv_arr[i] = res_c
    TTo_arr[i] = res_classic
    TT_arr[i] = res_TT


# Save for future plot
np.save('output/norm_10k_admx_05.npy',norm_arr)
np.save('output/bound_10k_admx_05.npy',bound_arr)
np.save('output/conv_10k_admx_05.npy',conv_arr)
np.save('output/TTo_10k_admx_05.npy',TTo_arr)
np.save('output/TT_10k_admx_05.npy',TT_arr)

elapsed_time = (time() - start_time)/60
print('Time elapsed in minutes: ', elapsed_time)
