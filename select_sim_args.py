"""
Uses latin hypercube sampling to decide which simulation datasets we should input for Random Forest Feature Selection
# Name	Cover	FragLenAve	FragLenStd	PopSize	TreeLen	ART_Profile	CodonSites	SameSeed	WindowSize	MinWinWidth	MinWinDepth	MinQual	TipSwap	Ns
"""
from pyDOE import *  # handles latin hypercube sampling
import csv
import os
from collections import namedtuple
import scipy.stats
import random
import math

MAX_SEED = math.pow(2, 32)-1  # internally asg_driver.py used numpy random generator which can only take up to 32bit seeds

# a csv that specifies the ranges for fields dictating a simulated dataset
# Use the ranges to randomly pick field values for each simulated dataset
SIM_DATA_DIR =  os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations/data"
SIM_FIELD_RANGE_CSV = "./sim_config/sim_args_range.csv"
SIM_ARGS_TSV = "./sim_config/sim_args.tsv"  # a tsv that specifies the actual field values for each simulated dataset
NUM_DATASETS = 12

LONGSHOT_ART_PROFILE = os.path.abspath(os.path.dirname(__file__) + os.sep + "art_profiles" + os.sep +  "longshot_ART_profile.")
CODON_SITES = 300
POPN_SIZE = 1000
SIM_DATASET_NAME_PREFIX = "Sim"

CONSTANT_FIELDS = ["Name",
                   "PopSize",
                   "CodonSites",
                   "ART_Profile",
                   # Technically this is not a constant field, but we don't want to use the latin hypercube sampling to create it
                   # because we don't care how correlated the seed is to other fields
                   "Seed"
]

# Dataset simulation argument values that will be set according to random selection,
# but not through latin hypercube sampling
RANDOM_FIELDS = ["TotalMutRates"]

# Dataset simulation argument values that will be randomly set according to latin hypercube sampling algorithm
LHS_FIELDS = ["Cover",  # fraction of individuals covered by a read
              "FragLenAve",
              "FragLenStd",
              "WindowSize",
              "MinWinWidth",  # fraction of window size required to have true base in sequence
              "MinWinDepth",  # fraction of indiv, or fraction of coverage, whichever is smaller
              "MinQual",
              "NumBreakpoints",   # total recombination breakpoints in genome
              "Generations",
              "SelectionRate"
]

LHS_MUTRATE_FIELDS = ["TreeLen"]  # total substitutions/site for full population tree

FieldRange = namedtuple("FieldRange", field_names=["Min", "Max", "Is_Int"])

# FieldName	Min	Max Is_Int
rand_field_ranges = dict()
with open(SIM_FIELD_RANGE_CSV, 'rU') as fh_in:
    reader = csv.DictReader(filter(lambda row: not row.startswith("#"), fh_in))
    for row in reader:
        field = row["FieldName"]
        rand_field_ranges[field] = FieldRange(Min=float(row["Min"]), Max=float(row["Max"]), Is_Int=int(row["Is_Int"]) != 0)

# Double check that we don't have typos in our args and that we have accounted for all the simulations args
extra_rand_fields = [x for x in rand_field_ranges.keys() if x not in LHS_FIELDS + RANDOM_FIELDS + LHS_MUTRATE_FIELDS]
if len(extra_rand_fields) > 0:
    raise ValueError("There are specified ranges for simulation arguments that shouldn't be randomized in " + SIM_FIELD_RANGE_CSV +
                     ": " + str(extra_rand_fields))

missing_rand_fields = [x for x in LHS_FIELDS+RANDOM_FIELDS+LHS_MUTRATE_FIELDS if x not in rand_field_ranges.keys()]
if len(missing_rand_fields) > 0:
    raise ValueError("There are no specified ranges for simulation arguments that should be randomized in " + SIM_FIELD_RANGE_CSV +
                     ": " + str(missing_rand_fields))


# We want to randomize the total distinct mutation rates.
# The problem is that we want each mutation rate needs to be sampled from a latin hypercube,
# and the dimensions of the latin hypercube change with the value of the total distinct mutation rates.
# Thus we can't sample the total distinct mutation rates in the latin hypercube.
# Instead, we randomly select the total distinct mutation rates first.
# Then we build the latin hypercube corresonding to the total distinct mutation rates.

if "TotalMutRates" not in rand_field_ranges:
    raise ValueError("Must specify range of  TotalMutRates in " + SIM_FIELD_RANGE_CSV)

# [dataset index] [fields sampled from latin hypercube for that dataset]
dataset_lhs = []

# dataset index ==> total mutation rates for that dataset
dataset_total_mut_rates = []

for dataset_idx in xrange(NUM_DATASETS):
    total_mut_rates = random.randint(rand_field_ranges["TotalMutRates"].Min, rand_field_ranges["TotalMutRates"].Max)
    dataset_total_mut_rates.extend([total_mut_rates])

    total_fields = len(LHS_FIELDS) + total_mut_rates
    # 1 x total_fields numpy array of values in [0, 1]
    dataset_lhs.append(lhs(n=total_fields, samples=1)[0])

print dataset_lhs


# Create a sim_args_tsv
with open(SIM_ARGS_TSV, 'w')  as fh_out:
    writer = csv.DictWriter(fh_out, fieldnames=CONSTANT_FIELDS+LHS_FIELDS+LHS_MUTRATE_FIELDS, delimiter="\t")
    writer.writeheader()

    for dataset_idx in xrange(NUM_DATASETS):
        outrow = dict()

        # Use the same dataset name prefix and append seed to make a unique id
        seed = random.randint(0, MAX_SEED)
        dataset_name = SIM_DATASET_NAME_PREFIX +  str(seed)

        # Keep regenerating a uid if the dataset folder with the same uid exists
        dataset_dir = SIM_DATA_DIR + os.sep + dataset_name
        while os.path.exists(dataset_dir):
            seed = random.randint(0, MAX_SEED)
            dataset_name = SIM_DATASET_NAME_PREFIX +  str(seed)

        outrow["Name"] = dataset_name
        outrow["Seed"] = seed

        # Simulation arguments that will be constant between all datasets
        outrow["PopSize"] = POPN_SIZE
        outrow["ART_Profile"] = LONGSHOT_ART_PROFILE
        outrow["CodonSites"] = CODON_SITES



        # Output all Dataset simulation arguments whose values are sampled from latin hypercube, except TreeLen
        lhs_vals = dataset_lhs[dataset_idx]
        for field_idx, rand_field_name in enumerate(LHS_FIELDS):
            # The latin hypercube specifies a scaling in between [0, 1]
            # Convert the scaling from [0, 1] to the field's range [min, max]
            scale = lhs_vals[field_idx]
            field_val_dist = rand_field_ranges[rand_field_name].Max - rand_field_ranges[rand_field_name].Min
            scaled_field_val = rand_field_ranges[rand_field_name].Min + (field_val_dist * scale)


            if rand_field_ranges[rand_field_name].Is_Int:
                outrow[rand_field_name] = int(round(scaled_field_val))
            else:
                outrow[rand_field_name] = scaled_field_val


        # TreeLen is a special Field
        # There will be total_mut_rates samples in the dataset's latin hypercube for TreeLen.
        # We output comma-separated list of each TreeLen value into the TSV
        total_mut_rates = dataset_total_mut_rates[dataset_idx]
        if (len(lhs_vals) - len(LHS_FIELDS)) != total_mut_rates:
            raise ValueError("We expect " + str(total_mut_rates) + " for dataset " + dataset_name  + " with 0based index " +
                             str(dataset_idx) + " but " + str(len(lhs_vals) - len(LHS_FIELDS)) + " were generated")

        mutation_rates = []
        for field_idx_offset in range(total_mut_rates):
            # The latin hypercube specifies a scaling in between [0, 1]
            # Convert the scaling from [0, 1] to the field's range [min, max]
            scale = lhs_vals[len(LHS_FIELDS) + field_idx_offset ]
            field_val_dist = rand_field_ranges["TreeLen"].Max - rand_field_ranges["TreeLen"].Min
            scaled_field_val = rand_field_ranges["TreeLen"].Min + (field_val_dist * scale)
            mutation_rates.extend([scaled_field_val])

        outrow["TreeLen"] = ",".join([str(x) for x in mutation_rates])


        # MinWinDepth is a special field.
        # We can't expect the window to cover more depth than is available in reads that fit the breadth threshold.
        # We also don't need the window to cover more than number of individuals in the population,
        #   because perfect reads covering the same individual will be duplicates of each other.
        # The highest possible threshold for min reads per window is the min of both.
        # In the sim_args_ranges.tsv file, MinWinDepth is set to the fraction of that highest possible threshold.


        # Find the fraction of reads that exceed the window breadth threshold
        min_win_width_bp = int(outrow["MinWinWidth"] * outrow["WindowSize"])
        frag_ave = outrow["FragLenAve"]
        frag_std = outrow["FragLenStd"]
        prob_exceed_width_thresh = 1 - scipy.stats.norm.cdf(x=min_win_width_bp, loc=frag_ave, scale=frag_std)

        # Find the highest possible window depth threshold that makes sense
        popsize = outrow["PopSize"]
        cov_depth = outrow["Cover"]
        if cov_depth >= 1:
            highest_win_depth_thresh = prob_exceed_width_thresh * popsize
        else:
            highest_win_depth_thresh = prob_exceed_width_thresh * cov_depth * popsize

        fraction_win_depth_thresh = outrow["MinWinDepth"]
        min_win_read_depth = int(fraction_win_depth_thresh * highest_win_depth_thresh)
        outrow["MinWinDepth"] = min_win_read_depth

        writer.writerow(outrow)













