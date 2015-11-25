"""
Runs SlidingWindow/test/simulations/sim_pipeline.py in mpi mode.
Reads in the Umberjack_Benchmark/sim_config/sim_args.tsv file to find the datasets.
Assumes all simulated dataset configs will be under Umberjack_Benchmark/simulations/data/<datasetname>/dataset.config

"""
from mpi4py import MPI
import os
import csv
import subprocess
import sys
import logging
import config.settings

config.settings.setup_logging()

LOGGER = logging.getLogger(__name__)


# Get the simulation datasets from Umberjack_Benchmark/sim_config/sim_args.tsv
SIM_ARGS_TSV = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "sim_config" + os.sep + "sim_args.tsv")
SIM_DATA_DIR = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations" + os.sep + "data")
SIM_PIPELINE_EXE = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + os.pardir + os.sep +
                                   "SlidingWindow/test/simulations/sim_pipeline.py")

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

LOGGER.debug("I am rank " + str(my_rank) + " of " + str(nprocs))

if len(sys.argv) < 2:
    sim_args_tsv = SIM_ARGS_TSV
else:
    sim_args_tsv = sys.argv[1]

with open(sim_args_tsv, 'rU') as fh_in:
    reader = csv.DictReader(fh_in, delimiter="\t")
    for i, row in enumerate(reader):
        name = row["Name"]
        dataset_config_file = SIM_DATA_DIR + os.sep + name + os.sep + name + ".config"

        if i % nprocs != my_rank or my_rank > i:
            continue

        subprocess.check_call(["python", SIM_PIPELINE_EXE, dataset_config_file])