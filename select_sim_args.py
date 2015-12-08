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

MATE_LEN_BP = 251

# a csv that specifies the ranges for fields dictating a simulated dataset
# Use the ranges to randomly pick field values for each simulated dataset
SIM_DATA_DIR =  os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations/data"
SIM_FIELD_RANGE_CSV = "./sim_config/sim_args_range.csv"
SIM_ARGS_TSV = "./sim_config/sim_args.tsv"  # a tsv that specifies the actual field values for each simulated dataset
NUM_DATASETS = 50

# Use path relative to directory of simulated dataset config file
# /Umberjack_Benchmark/simulations/data/SimDataset/SimDataset.config
LONGSHOT_ART_PROFILE = os.pardir + os.sep + os.pardir + os.sep + os.pardir + os.sep + "art_profiles" + os.sep +  "longshot_ART_profile."
CODON_SITES = 300
POPN_SIZE = 1000
SELECTION_RATE = 0.01
NUM_MUTATION_RATES = 3
SIM_DATASET_NAME_PREFIX = "Sim"
FRAG_STD = 100

CONSTANT_FIELDS = ["Name",
                   "PopSize",
                   "CodonSites",
                   "ART_Profile",
                   # Technically this is not a constant field, but we don't want to use the latin hypercube sampling to create it
                   # because we don't care how correlated the seed is to other fields
                   "Seed",
                   "SelectionRate",
                   "FragLenStd"
]

# Dataset simulation argument values that will be randomly set according to latin hypercube sampling algorithm
LHS_FIELDS = ["Cover",  # fraction of individuals covered by a read
              "FragLenAve",
              "RecomboRate",  # Recombinations Per Genome Codon  (Only calculated once for the entire genome, not repeated for each individual)
              "Generations",
              "MutationRate1", "MutationRate2", "MutationRate3",
              "UmberjackConfig"
]

# Randomly use one of two umberjack configurations
UMBERJACK_CONFIGS = [{
                        "WindowSize": 150,
                        "MinWinWidth": 0.7,
                        "MinWinDepth": 10,
                        "MinQual": 15
                    },
                    {
                        "WindowSize": 300,
                        "MinWinWidth": 0.875,
                        "MinWinDepth": 10,
                        "MinQual": 20
                    }]

#LHS_MUTRATE_FIELDS = ["MutationRate"]  # Nucleotide substitutions across tree per generation per individual


INPUT_COLS = LHS_FIELDS
OUTPUT_COLS = CONSTANT_FIELDS + ["FragLenAve",
                                 "Cover",
                                 "NumBreakpoints",
                                 "Generations",
                                 "TreeLen",
                                 "WindowSize",
                                 "MinWinWidth",
                                 "MinWinDepth",
                                 "MinQual"]

FieldRange = namedtuple("FieldRange", field_names=["Min", "Max", "Is_Int"])

# Read in the ranges for the latin hypercube samples
# FieldName	Min	Max Is_Int
rand_field_ranges = dict()
with open(SIM_FIELD_RANGE_CSV, 'rU') as fh_in:
    reader = csv.DictReader(filter(lambda row: not row.startswith("#"), fh_in))
    for row in reader:
        field = row["FieldName"]
        rand_field_ranges[field] = FieldRange(Min=float(row["Min"]), Max=float(row["Max"]), Is_Int=int(row["Is_Int"]) != 0)

# Double check that we don't have typos in our args and that we have accounted for all the simulations args
extra_rand_fields = [x for x in rand_field_ranges.keys() if x not in INPUT_COLS]
if len(extra_rand_fields) > 0:
    raise ValueError("There are specified ranges for simulation arguments that shouldn't be randomized in " + SIM_FIELD_RANGE_CSV +
                     ": " + str(extra_rand_fields))

missing_rand_fields = [x for x in INPUT_COLS if x not in rand_field_ranges.keys()]
if len(missing_rand_fields) > 0:
    raise ValueError("There are no specified ranges for simulation arguments that should be randomized in " + SIM_FIELD_RANGE_CSV +
                     ": " + str(missing_rand_fields))

# Don't bother building multiple latin hypercubes, each representing different number of mutation rates.
# We will always have 3 mutation rates.
# lhs() needs to do all the sampling at once to ensure that the overlap in field-space is minimal.
# That means the dimensions of the hypercube must be constant.
total_fields = len(LHS_FIELDS)
# Dimensions:  [total datasets ] x [total fields]
# [dataset index] [field index] = numpy array element representing field value in [0, 1]
dataset_lhs = lhs(n=total_fields, samples=NUM_DATASETS, criterion="maximin")
print dataset_lhs


# Create a sim_args_tsv
is_append = False
fieldnames = None
if os.path.exists(SIM_ARGS_TSV) and os.path.getsize(SIM_ARGS_TSV):
    is_append = True
    with open(SIM_ARGS_TSV, 'r')  as fh_in:
        reader = csv.DictReader(fh_in, delimiter="\t")
        fieldnames = reader.fieldnames

with open(SIM_ARGS_TSV, 'a')  as fh_out:
    fieldnames = (fieldnames or OUTPUT_COLS)
    writer = csv.DictWriter(fh_out, fieldnames=OUTPUT_COLS, delimiter="\t")
    if not is_append:
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
        outrow["SelectionRate"] = SELECTION_RATE
        outrow["FragLenStd"] = FRAG_STD

        # Output all Dataset simulation arguments whose values are sampled from latin hypercube
        lhs_vals = dataset_lhs[dataset_idx]
        mutation_rates = []
        for field_idx, rand_field_name in enumerate(LHS_FIELDS):
            # The latin hypercube specifies a scaling in between [0, 1]
            # Convert the scaling from [0, 1] to the field's range [min, max]
            scale = lhs_vals[field_idx]
            field_val_dist = rand_field_ranges[rand_field_name].Max - rand_field_ranges[rand_field_name].Min
            scaled_field_val = rand_field_ranges[rand_field_name].Min + (field_val_dist * scale)

            if rand_field_ranges[rand_field_name].Is_Int:
                scaled_field_val = int(round(scaled_field_val))



            if rand_field_name in OUTPUT_COLS:
                outrow[rand_field_name] = scaled_field_val

            # Override Cover
            if rand_field_name == "Cover":
                outrow["Cover"] = 2 ** scaled_field_val

            elif rand_field_name == "Generations":
                outrow["Generations"] = int(round(2 ** scaled_field_val))

            # Keep track of all the mutation rates to output TreeLen field
            elif rand_field_name.find("MutationRate") >= 0:
                mutation_rates.extend([scaled_field_val])

            # Override Umberjack config fields
            elif rand_field_name == "UmberjackConfig":
                outrow.update(UMBERJACK_CONFIGS[scaled_field_val])

            elif rand_field_name == "RecomboRate":  # recombo rate is fraction of genome codons with recombination breakpoints
                outrow["NumBreakpoints"] = int(round(scaled_field_val * CODON_SITES))

        # Override TreeLen
        # Mutation rates are total substitutions/site across phylogeny per generation per individual
        outrow["TreeLen"] = ",".join([str(rate * POPN_SIZE * outrow["Generations"]) for rate in mutation_rates])


        # MinWinDepth is a special field.
        # We can't expect the window to cover more depth than is available in reads that fit the breadth threshold.
        # We also don't need the window to cover more than number of individuals in the population,
        #   because perfect reads covering the same individual will be duplicates of each other.
        # The highest possible threshold for min reads per window is the min of both.
        # In the sim_args_ranges.tsv file, MinWinDepth is set to the fraction of that highest possible threshold.

        # Find the fraction of fragments that exceed the window breadth threshold
        # min_win_width_bp = int(outrow["MinWinWidth"] * outrow["WindowSize"])
        # frag_ave = outrow["FragLenAve"]
        # frag_std = outrow["FragLenStd"]
        #
        # # When the fragment size is longer than 2xmate length, there will be a gap between the mates.
        # # That gap shouldn't be bigger than the max allowable gap
        # max_gap_bp = outrow["WindowSize"] - min_win_width_bp
        # max_frag_len_fit_win = 2 * MATE_LEN_BP + max_gap_bp
        #
        #
        # # Find the probability of a window starting at a position in the fragments such that the fragment
        # # passes the window breadth threshold.
        # # P(fragment fits window) = sum over fragment size {P(fragment fits window | fragment size) * P (fragment size)}
        # p_frag_fit_win = 0.0
        # for frag_size in xrange(min_win_width_bp, max_frag_len_fit_win):
        #     # the fraction of window start positions wrt fragment that allow the fragment to fit in the window
        #     p_frag_fits_win_given_fragsize = (frag_size - min_win_width_bp + 1.0)/float(frag_size)
        #     p_fragsize = scipy.stats.norm.pdf(x=frag_size, loc=frag_ave, scale=frag_std)
        #     p_frag_fit_win += p_frag_fits_win_given_fragsize * p_fragsize
        #
        #
        # # Find the highest possible window depth threshold that makes sense
        # popsize = outrow["PopSize"]
        # cov_depth = outrow["Cover"]  # how many times each genome is covered per individual
        #
        #
        # # At 1x coverage, we expect a fragment covering each site for each individual
        # if cov_depth >= 1:
        #     highest_win_depth_thresh = p_frag_fit_win * popsize
        # else:
        #     highest_win_depth_thresh = p_frag_fit_win * cov_depth * popsize
        #
        # fraction_win_depth_thresh = outrow["MinWinDepth"]
        # min_win_read_depth = int(fraction_win_depth_thresh * highest_win_depth_thresh)
        # outrow["MinWinDepth"] = min_win_read_depth

        writer.writerow(outrow)













