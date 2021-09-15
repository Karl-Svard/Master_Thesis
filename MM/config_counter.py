import msprime as msp
import numpy as np
import equations_for_a_b_tauA_tauB as eq


# Population 1, 2 and the outgroup
P1, P2, O = 0, 1, 2

def set_up_pops(nS, tS, n0):
    # Samples from population 1
    samples = [msp.Sample(population=P1, time=tS[0])]*(1*nS[0])
    # Samples from population 2
    samples.extend([msp.Sample(population=P2, time=tS[0])]*(1*nS[1]))
    # Samples from outgroup
    #samples.extend([msp.Sample(population=O, time=tS[0])]*(1*n0))
    return samples

def set_up_demography(t12, t120, f):
    # determine topology of tree (divergences)
    divergence = [msp.MassMigration(time=t12,source=P2,destination=P1,proportion=1)]
                  #msp.MassMigration(time=t120,source=O,destination=P1,proportion=1)]
    # admixture times and proportion
    #admix = [msp.MassMigration(time=ta,source=P1,destination=P2,proportion=f)]
    admix2 = [msp.MigrationRateChange(time=0.95*t12, rate=f, matrix_index=(P1,P2)),
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
#r=0
generation_time = 25
t12 = 50000/generation_time #50000 years ago in generations
t120 = 100000/generation_time #80000 years ago in generations
#ta = 45000/generation_time
nS=[4,4] # number of samples per population
n0=1 # number of samples from the outgroup
tS=[0] #all samples taken in modern times?
f=0.2
N=[Ne]*2 #effective population sizes of the three pops
seed=None
migration_matrix = [
    [0, 0],
    [0, 0]] # original migration matrix

reps = 1000

samples = set_up_pops(nS,tS,n0)
demography = set_up_demography(t12, t120, f)
pops = [msp.PopulationConfiguration(initial_size = n) for n in N]

#dd = msp.DemographyDebugger(population_configurations=pops,
#                            demographic_events=demography,
#                            migration_matrix=migration_matrix)
#dd.print_history()

sims = msp.simulate(samples=samples,Ne=N[0],population_configurations=pops,
                    demographic_events = demography, mutation_rate = mu,
                    length=L, recombination_rate = r, random_seed = seed,
                    migration_matrix=migration_matrix, num_replicates=reps)

#configs = [[(1, 3, 1, 3), (1, 3, 2, 2), (1, 3, 3, 1), (1, 3, 4, 0), (2, 2, 1, 3),
#            (2, 2, 2, 2), (2, 2, 3, 1), (2, 2, 4, 0), (3, 1, 1, 3), (3, 1, 2, 2),
#            (3, 1, 3, 1), (3, 1, 4, 0), (4, 0, 1, 3), (4, 0, 2, 2), (4, 0, 3, 1)],
#           [(0, 4, 0, 4), (0, 4, 1, 3), (0, 4, 2, 2), (0, 4, 3, 1), (0, 4, 4, 0),
#            (1, 3, 0, 4), (2, 2, 0, 4), (3, 1, 0, 4), (4, 0, 0, 4), (4, 0, 4, 0)]]

configs = eq.get_sample_sets(4,4)


# Rewrite later so that the contigs are calculated from the function
def config_counts(configs, sims, nr_samp):
    nr_conf = len(configs[0])
    counts = [0]*nr_conf

    for variant in sims.variants():
        geno = variant.genotypes
        p1_der = sum(geno[0:nr_samp[0]])
        p2_der = sum(geno[nr_samp[0]:])

        obs_conf = (p1_der, nr_samp[0]-p1_der, p2_der, nr_samp[1]-p2_der)
        if obs_conf in configs[0]:
            index = configs[0].index(obs_conf)
            counts[index] = counts[index] + 1

    return counts

res_list = np.empty((reps,4))
bval_list = np.empty((reps))
for i, sim in enumerate(sims):
    M = config_counts(configs, sim, [4,4])
    print(M)
    nA=4
    nB=4
    res_list[i] = eq.find_optimum(nA,nB,configs[0],M).x
    bval_list[i] = eq.f_alt(res_list[i])

print(res_list.mean(axis=0))
print(bval_list.mean())

np.save('output/mm_results_2admx.npy',res_list)
np.save('output/bval_2admx.npy',bval_list)
