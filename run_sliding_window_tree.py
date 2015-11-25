import os
import subprocess
import glob
import re
import logging
import random
import sys
import Bio.SeqIO as SeqIO
import csv
from pool import pool_traceback
import config.settings as settings
from collections import namedtuple
import math

settings.setup_logging()

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = 1
MAX_SEED = math.pow(2, 32)-1  # internally asg_driver.py used numpy random generator which can only take up to 32bit seeds

SIM_ARGS_TSV = "./sim_config/sim_args.tsv"  # a csv that specifies which simulations to run

SIM_DATA_DIR =  os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations/data"
SIM_OUT_DIR =   os.path.dirname(os.path.realpath(__file__)) + os.sep +"simulations/out"

MAPQ_CUTOFF = 20  # alignment quality cutoff
# We set this to 1 because we know that our error rate is good.  We only want breadth to determine sliding windows.
MAX_PROP_N = 1  # maximum proportion of N bases in MSA-aligned sequence



CONCURRENT_MPIRUN = 3
THREADS_PER_WINDOW = 4
WINDOW_PROCS = 3
WINDOW_SLIDE = 30
PROCS = 20

MASK_STOP_CODON = True
REMOVE_DUPLICATES = True
KEEP_INSERTS = False

REF = "consensus"




#TEST_PREFIX_FORMAT = "{}.cov{}.frag{}_{}.indiv{}.codon{}.scale{}"
PopnGroup = namedtuple("PopnGroup", field_names=["dataset", "art_profile", "config_file", "indiv", "codonsites",
                                                         "cov_depth", "scales", "frag_ave", "frag_std",
                                                         "breakpoints", "num_breakpoints", "seed",
                                                         "selection_rate", "generations"])

UmberjackGroup = namedtuple("UmberjackGroup", field_names=["window_size", "min_win_width", "min_win_depth", "min_qual"])


def get_umberjack_outdir_from_simargs_tsv(popn_group, umberjack_group):
    """
    Gets the umberjack output directory for the dataset as specifed in run_sliding_window_tree.py in the input_csv file
    :param popn_group:
    :param umberjack_group:
    :return:
    """

    OUT_DIR =   (SIM_OUT_DIR + os.sep +
                 "window{}.breadth{}.depth{}.qual{}".format(umberjack_group.window_size,
                                                            umberjack_group.min_win_width,
                                                            umberjack_group.min_win_depth,
                                                            umberjack_group.min_qual))

    sample_out_dir = (OUT_DIR + os.sep + popn_group.dataset )

    return sample_out_dir


def get_sim_dataset_dir(popn_group):
    """
    For the given simulated dataset, return its directory
    :param popn_group:
    :return:
    """
    return SIM_DATA_DIR + os.sep + popn_group.dataset


def do_sliding_window(outdir, input_csv,
                      window_size, window_depth_cutoff, window_breadth_cutoff, min_qual,
                      machine_file):



    call = ["mpirun",
            "--machinefile", machine_file,
            "--output-filename", outdir + os.sep + "logs" + os.sep + os.path.basename(outdir) + ".log",
            "python",
            os.path.dirname(os.path.realpath(__file__)) + os.sep + "../SlidingWindow/umberjack.py",
            "--out_dir", outdir,
            "--map_qual_cutoff", str(MAPQ_CUTOFF),
            "--read_qual_cutoff", str(min_qual),
            "--max_prop_n", str(MAX_PROP_N),
            "--window_size", str(window_size),
            "--window_slide", str(WINDOW_SLIDE),
            "--window_breadth_cutoff", str(window_breadth_cutoff),
            "--window_depth_cutoff", str(window_depth_cutoff),
            "--threads_per_window", str(THREADS_PER_WINDOW),
            "--mpi",
            "--mode", "DNDS",
            "--input_csv", input_csv
    ]
    if MASK_STOP_CODON:
        call = call + ["--mask_stop_codon"]
    if REMOVE_DUPLICATES:
        call = call + ["--remove_duplicates"]
    if KEEP_INSERTS:
        call = call + ["--insert"]


    LOGGER.debug("About to execute cmd:" + " ".join(call))
    try:
        subprocess.check_call(call, shell=False, env=os.environ)
    except Exception, e:
        LOGGER.error("Failure in " + " ".join(call) + "\n" + e.message)
        LOGGER.exception(e.message)

    # rconfig_file = os.path.dirname(os.path.realpath(__file__)) + os.sep +"simulations" + os.sep + "R" + os.sep + "umberjack_unit_test.config"
    # with open(rconfig_file, 'w') as fh_out_config:
    #     fh_out_config.write("ACTUAL_DNDS_FILENAME=" + output_csv + "\n")
    #     fh_out_config.write("EXPECTED_DNDS_FILENAME=" + expected_dnds_filename + "\n")
    #     fh_out_config.write("EXPECTED_DNDS_START_NUC_POS=" + str(START_NUCPOS) + "\n")
    #     fh_out_config.write("EXPECTED_DNDS_END_NUC_POS=" + str(ref_len) + "\n")
    #     fh_out_config.write("INDELIBLE_DNDS_FILENAME=" + indelible_dnds_filename + "\n")
    #
    # Rscript_wdir =  os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations" + os.sep + "R")
    # subprocess.check_call(["Rscript", "-e",
    #                        ("library(knitr); " +
    #                         "setwd('{}'); ".format(Rscript_wdir) +
    #                         "spin('umberjack_unit_test.R', knit=FALSE); " +
    #                         "knit2html('./umberjack_unit_test.Rmd', stylesheet='./markdown_bigwidth.css')")],
    #                       shell=False, env=os.environ)
    # shutil.copy(Rscript_wdir + os.sep + "umberjack_unit_test.html",
    #             outdir + os.sep + "umberjack_unit_test.html")



# def do_collate(outdir, output_csv, ref_fasta, full_popn_fasta,expected_dnds_filename, indelible_dnds_filename,
#                full_popn_conserve_csv, orig_conserve_csv, aln_conserve_csv):
#
#     # collect_stats.collect_dnds(output_dir=outdir, output_csv_filename=output_csv, full_popn_fasta=full_popn_fasta)
#
#     # Rscript_wdir =  os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "R")
#     #
#     # rcollate_config_file = Rscript_wdir + os.sep + "aggreg_window.config"
#     # with open(rcollate_config_file, 'w') as fh_out_config:
#     #     fh_out_config.write("COLLATE_DNDS_FILENAME=" + output_csv + "\n")
#     #     fh_out_config.write("EXPECTED_DNDS_FILENAME=" + expected_dnds_filename + "\n")
#     #     fh_out_config.write("EXPECTED_DNDS_START_NUC_POS=" + str(START_NUCPOS) + "\n")
#     #     fh_out_config.write("EXPECTED_DNDS_END_NUC_POS=" + str(ref_len) + "\n")
#     #     fh_out_config.write("INDELIBLE_DNDS_FILENAME=" + indelible_dnds_filename + "\n")
#     #     fh_out_config.write("SMOOTH_DIST=" + str(SMOOTH_DIST) + "\n")
#     #     fh_out_config.write("FULL_POPN_CONSERVE_CSV=" + full_popn_conserve_csv + "\n")
#     #     fh_out_config.write("ORIG_CONSERVE_CSV=" + orig_conserve_csv + "\n")
#     #     fh_out_config.write("ALN_CONSERVE_CSV=" + aln_conserve_csv + "\n")
#
#
#     # subprocess.check_call(["Rscript", "-e",
#     #                        ("library(knitr); " +
#     #                         "setwd('{}'); ".format(Rscript_wdir) +
#     #                         "spin('aggreg_window.R', knit=FALSE); " +
#     #                         "knit2html('./aggreg_window.Rmd', stylesheet='./markdown_bigwidth.css')")],
#     #                       shell=False, env=os.environ)
#     # shutil.copy(Rscript_wdir + os.sep + "aggreg_window.html",
#     #             outdir + os.sep + "aggreg_window.html")

def filter_fasta(in_fasta, keep_names, out_fasta):
    """
    Keeps only the sequences from fasta specified in keep_names
    :param str in_fasta: file path to input fasta
    :param set keep_names:  Set of sequences names to keep
    :param str out_fasta:  filepath to output fasta
    :return:
    """
    orig_ids = set()
    with open(out_fasta, 'w') as fh_out:
        for record in SeqIO.parse(in_fasta, "fasta"):
            orig_ids.add(record.id)
            if record.id in keep_names:
                fh_out.write(">" + record.id + "\n")
                fh_out.write(str(record.seq) + "\n")

    for keep_name in keep_names:
        if keep_name not in orig_ids:
            print "!!!Record " + keep_name + " from keep_names not in " + in_fasta


def get_id_set(fasta):
    """
    :param str fasta:  filepath to fasta
    :return:  Set of names from a fasta
    """
    ids = set()
    for record in SeqIO.parse(fasta, "fasta"):
        ids.add(record.id)
    return ids


def downsample_errfree_window_by_typical(typical_outdir, errfree_outdir, downsample_errfree_outdir):
    """
    For the same window, removes the error free sequences that don't exist in the typical sequences.

    :param str typical_outdir:  path to typical error sequence window slice fastas
    :param str errfree_outdir: path to error free sequence window slice fastas
    :param str downsample_errfree_outdir: path to downsampled error free sequence window slice fastas
    """
    if not os.path.exists(downsample_errfree_outdir):
        os.makedirs(downsample_errfree_outdir)

    for typical_slice_fasta_filename in glob.glob(typical_outdir + os.sep + "*.fasta"):
        print "Handling " + typical_slice_fasta_filename
        errfree_slice_fasta_basename = os.path.basename(typical_slice_fasta_filename).replace(".reads", ".reads.errFree")
        errfree_slice_fasta_filename = errfree_outdir + os.sep + errfree_slice_fasta_basename
        typical_ids = get_id_set(typical_slice_fasta_filename)
        downsample_errfree_slice_fasta_filename = downsample_errfree_outdir + os.sep + errfree_slice_fasta_basename
        filter_fasta(errfree_slice_fasta_filename, typical_ids, downsample_errfree_slice_fasta_filename)



def downsample_Nmask_errfree_window_by_typical(typical_outdir, errfree_outdir, downsample_errfree_outdir):
    """
    For the same window, removes the error free sequences that don't exist in the typical sequences.
    N-masks the same bases from the remaining error free sequences at the same positions as the typical sequences

    :param str typical_outdir:  path to typical error sequence window slice fastas
    :param str errfree_outdir: path to error free sequence window slice fastas
    :param str downsample_errfree_outdir: path to downsampled error free sequence window slice fastas
    """
    if not os.path.exists(downsample_errfree_outdir):
        os.makedirs(downsample_errfree_outdir)

    for typical_slice_fasta_filename in glob.glob(typical_outdir + os.sep + "*.fasta"):
        print "Handling " + typical_slice_fasta_filename
        errfree_slice_fasta_basename = os.path.basename(typical_slice_fasta_filename).replace(".reads", ".reads.errFree")
        errfree_slice_fasta_filename = errfree_outdir + os.sep + errfree_slice_fasta_basename
        downsample_errfree_slice_fasta_filename = downsample_errfree_outdir + os.sep + errfree_slice_fasta_basename
        with open(downsample_errfree_slice_fasta_filename, 'w') as fh_out:
            errfree_recdict = SeqIO.to_dict(SeqIO.parse(errfree_slice_fasta_filename, "fasta"))
            for typical_rec in SeqIO.parse(typical_slice_fasta_filename, "fasta"):
                if not errfree_recdict.get(typical_rec.id):
                    print "!!!Record " + typical_rec.id + " from typical reads not in " + errfree_slice_fasta_filename
                errfree_seq = errfree_recdict[typical_rec.id].seq.tomutable()
                for match in re.finditer(r"N", str(typical_rec.seq)):
                    errfree_seq[match.start()] = "N"
                errfree_recdict[typical_rec.id].seq = errfree_seq
                fh_out.write(">" + typical_rec.id + "\n")
                fh_out.write(str(errfree_recdict[typical_rec.id].seq) + "\n")


def downsample_Nmask_pad_errfree_window_by_typical(typical_outdir, errfree_outdir, downsample_errfree_outdir):
    """
    For the same window, removes the error free sequences that don't exist in the typical sequences.
    N-masks the same bases from the remaining error free sequences at the same positions as the typical sequences.
    Left and right pads the same bases from the remaining error free sequences at the same pad positions as the typical sequences.

    :param str typical_outdir:  path to typical error sequence window slice fastas
    :param str errfree_outdir: path to error free sequence window slice fastas
    :param str downsample_errfree_outdir: path to downsampled error free sequence window slice fastas
    """
    if not os.path.exists(downsample_errfree_outdir):
        os.makedirs(downsample_errfree_outdir)

    for typical_slice_fasta_filename in glob.glob(typical_outdir + os.sep + "*.fasta"):
        print "Handling " + typical_slice_fasta_filename
        errfree_slice_fasta_basename = os.path.basename(typical_slice_fasta_filename).replace(".reads", ".reads.errFree")
        errfree_slice_fasta_filename = errfree_outdir + os.sep + errfree_slice_fasta_basename
        downsample_errfree_slice_fasta_filename = downsample_errfree_outdir + os.sep + errfree_slice_fasta_basename
        with open(downsample_errfree_slice_fasta_filename, 'w') as fh_out:
            errfree_recdict = SeqIO.to_dict(SeqIO.parse(errfree_slice_fasta_filename, "fasta"))
            for typical_rec in SeqIO.parse(typical_slice_fasta_filename, "fasta"):
                if not errfree_recdict.get(typical_rec.id):
                    print "!!!Record " + typical_rec.id + " from typical reads not in " + errfree_slice_fasta_filename
                errfree_seq = errfree_recdict[typical_rec.id].seq.tomutable()

                # Find the position of the first non-gap char.  Anything before this is a left-pad gap.
                typical_truebase_start = re.search(r"[^\-]", str(typical_rec.seq)).start()
                # Find the position of the last non-gap char.  Anything after this is a right-pad gap.
                typical_truebase_end = re.search(r"[^\-][\-]*$", str(typical_rec.seq)).start()

                # Check if the error free read is shorter than the typical read.  This is not expected
                # Find the position of the first non-gap char.  Anything before this is a left-pad gap.
                errfree_truebase_start = re.search(r"[^\-]", str(errfree_seq)).start()
                # Find the position of the last non-gap char.  Anything after this is a right-pad gap.
                errfree_truebase_end = re.search(r"[^\-][\-]*$", str(errfree_seq)).start()
                if errfree_truebase_start > typical_truebase_start or errfree_truebase_end < typical_truebase_end:
                    print "WARN:  read " + typical_rec.id + " typical is longer than error free"


                errfree_seq[:typical_truebase_start] = "-"*typical_truebase_start  # replace with leftpads
                errfree_seq[typical_truebase_end+1:] = "-"*(len(str(errfree_seq))-typical_truebase_end-1)  # replace with rightpads

                # N-mask
                for match in re.finditer(r"N", str(typical_rec.seq)):
                    errfree_seq[match.start()] = "N"
                errfree_recdict[typical_rec.id].seq = errfree_seq
                fh_out.write(">" + typical_rec.id + "\n")
                fh_out.write(str(errfree_recdict[typical_rec.id].seq) + "\n")


def gen_sim_data_helper(kw_namedtuple):
    """
    Helper for passing arguments to gen_sim_data() when calling in parallel
    :param kwargs:
    :return:
    """
    kwargs = kw_namedtuple._asdict()
    return gen_sim_data(**kwargs)


# Name	Cover	ReadLenAve	ReadLenStd	PopSize	Scale	Error	CodonSites	WindowSize	MinWinWidth	MinWinDepth	TipSwap	SameTree	Ns
def gen_sim_data(config_file,
                 art_profile, indiv, codonsites, cov_depth, scales,
                 frag_ave, frag_std, breakpoints, num_breakpoints,
                 selection_rate, generations,
                 seed=None, **kwargs):
    """
    Generate the config file and create the simulated dataset.

    All paths specified in the config file are relative to the config file's directory.

    :param str config_file: config file to output to
    :param str art_profile:  ART read simulator quality profile tsv
    :param int indiv:  total individals in population
    :param int codonsites: number of codon sites.
    :param int cov_depth:  fold read coverage  Must be >= 1
    :param [float] scales:  number of different mutation scalings.  Tree length will be scaled by this number in INDELible
    :param int frag_ave:  ART read fragment ave size
    :param int frag_std: ART read fragment standard deviation size
    :param [int] breakpoints:  list of 1-based codon position recombination breakpoints to explicitly set.  If not set, then they will be randomly chosen.
    :param int num_breakpoints:  total recombination breakpoints
    :param float selection_rate:  rate at which dominant mutation is selected in asg_driver.py
    :param int generations:  total generations in asg_driver.py
    :param int seed:  seed.  Should be < 2^32
    :param tuples kwargs:  extra arguments to pass

    """

    sim_outdir = os.path.dirname(config_file)
    filename_prefix = os.path.basename(config_file).split(".config")[0]
    if not os.path.exists(sim_outdir):
        os.makedirs(sim_outdir)

    dnds_tsv = sim_outdir + os.sep + "subs" + os.sep + filename_prefix + ".dnds.tsv"
    if os.path.exists(dnds_tsv) and os.path.getsize(dnds_tsv):
        LOGGER.warn("Not regenerating dataset for " + config_file)
    else:
        if not seed:
            seed = random.randint(0, MAX_SEED)

        LOGGER.debug("Creating simulation config " + config_file + " with seed " + str(seed))
        if codonsites % len(scales) != 0:
            raise ValueError("Number of scalings should divide evenly into number of codon sites")
        #codons_per_block = codonsites / len(scales)
        with open(config_file, 'w') as fh_out:
            fh_out.write("[sim]\n")
            fh_out.write("FILENAME_PREFIX={}\n".format(filename_prefix))
            fh_out.write("NUM_INDIV={}\n".format(indiv))
            fh_out.write("SEED={}\n".format(seed))
            fh_out.write("NUM_CODON_SITES={}\n".format(codonsites))

            fh_out.write("NUM_BREAKPOINTS={}\n".format(num_breakpoints))
            if breakpoints:
                fh_out.write("BREAKPOINTS={}\n".format(",".join([str(x) for x in breakpoints])))

            fh_out.write("SELECTION_RATE={}\n".format(selection_rate))
            fh_out.write("GENERATIONS={}\n".format(generations))

            fh_out.write("INDELIBLE_BIN_DIR=../../../../SlidingWindow/test/simulations/bin/indelible/indelible_1.03/linux_x64\n")
            fh_out.write("INDELIBLE_SCALING_RATES={}\n".format(",".join([str(x) for x in scales])))
            fh_out.write("ART_BIN_DIR = ../../../../SlidingWindow/test/simulations/bin/art/art_3.19.15_adapter/linux_x64\n")
            fh_out.write("ART_QUAL_PROFILE_TSV1 = {}R1.txt\n".format(art_profile))
            fh_out.write("ART_QUAL_PROFILE_TSV2 = {}R2.txt\n".format(art_profile))
            fh_out.write("ART_FOLD_COVER={}\n".format(cov_depth))
            fh_out.write("ART_MEAN_FRAG = {}\n".format(frag_ave))
            fh_out.write("ART_STDEV_FRAG = {}\n".format(frag_std))
            fh_out.write("ART_READ_LENGTH = 251\n")
            fh_out.write("ART_INSERT_RATE1 = 0.00045\n")
            fh_out.write("ART_INSERT_RATE2 = 0.00045\n")
            fh_out.write("ART_DEL_RATE1 = 0.00045\n")
            fh_out.write("ART_DEL_RATE2 = 0.00045\n")
            fh_out.write("ART_QUAL_SHIFT1 = 0\n")
            fh_out.write("ART_QUAL_SHIFT2 = 0\n")
            fh_out.write("PICARD_BIN_DIR = ../../../../SlidingWindow/test/simulations/bin/picard/picard_1.129\n")
            fh_out.write("BWA_BIN_DIR = ../../../../SlidingWindow/test/simulations/bin/bwa/bwa_0.7.12/linux_x64\n")
            fh_out.write("PROCS = {}\n".format(PROCS))
            fh_out.write("FASTTREE_EXE = ../../../../SlidingWindow/test/simulations/bin/fasttree/fasttree_2.1.7/linux_x64/FastTree\n")
            fh_out.write("HYPHY_EXE = ../../../../SlidingWindow/test/simulations/bin/hyphy/hyphy_2.2.3/linux_x64/HYPHYMP\n")
            fh_out.write("HYPHY_BASEPATH = ../../../../SlidingWindow/test/simulations/bin/hyphy/hyphy_2.2.3/res/TemplateBatchFiles\n")

        # simulations/data/dataset/subs/dataset.dnds.tsv

        sim_pipeline_exe = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + os.pardir + os.sep + "SlidingWindow/test/simulations/sim_pipeline.py")
        subprocess.check_call(["python", sim_pipeline_exe, config_file])





def process_window_size_helper((popn_umberjack_index, concur_mpi, umberjack_group, popn_groups)):
    """
    Helper function to call process_window_size with dict as args
    :param args:
    :return:
    """
    # use 0-based popn_umberjack_index to decide which machines to run on
    machine_file_index = popn_umberjack_index % concur_mpi

    machine_file = SIM_OUT_DIR + os.sep + "machine{}.txt".format(machine_file_index)

    # Convert UmberjackGroup namedtuple into dict so that we can convert into keyword args
    kwargs = umberjack_group._asdict()

    # Tack on the datasetes keyword arg, which is list of all the prefixes of the populations that we need to run through these umberjack settings
    test_prefixes = [popn_group.dataset for popn_group in popn_groups]
    kwargs["test_prefixes"] = test_prefixes
    kwargs["machine_file"] = machine_file

    return process_window_size(**kwargs)



def process_window_size(test_prefixes, window_size, min_win_width, min_win_depth, min_qual, machine_file, seed=None, **kwargs):
    """
    Umberjack output for the given configurations.
    Umberjack does a good job of  parallelizing multiple sams.  Let it do most of the work.
    We just need to ensure that all the sams that we want to process with the same umberjack settings
    (window size, min breadth, min depth, min qual) are together in the test_prefixes.
    :param str test_prefixes:
    :param int window_size:
    :param float min_win_width:
    :param int min_win_depth:
    :param int min_qual:
    :param float tip_swap:
    :param float missing:
    :param int seed:
    :param kwargs:
    :return:
    """


    UMBERJACK_OUT_DIR =  SIM_OUT_DIR + os.sep + "window{}.breadth{}.depth{}.qual{}".format(window_size, min_win_width, min_win_depth, min_qual)
    if not os.path.exists(UMBERJACK_OUT_DIR):
        os.makedirs(UMBERJACK_OUT_DIR)

    # typical reads
    umberjack_input_csv = UMBERJACK_OUT_DIR + os.sep + "umberjack_input.csv"

    with open(umberjack_input_csv, 'w') as fh_out:
        writer = csv.DictWriter(fh_out, fieldnames=["File", "Ref", "OutputPrefix"])
        writer.writeheader()

        for dataset in test_prefixes:

            # Run sliding windows and collect stats
            #FULL_POPN_FASTA =  SIM_DATA_DIR + os.sep + test_prefix + os.sep + "fullpopn" + os.sep + test_prefix + "_TRUE.fasta"
            #FULL_POPN_TREE =  SIM_DATA_DIR + os.sep + test_prefix + os.sep + "topology" + os.sep  + test_prefix + ".nwk"
            REFERENCE_FASTA =  SIM_DATA_DIR + os.sep + dataset + os.sep + dataset + ".consensus.fasta"

            FULL_POPN_CONSERVE_CSV = REFERENCE_FASTA.replace(".consensus.fasta", ".conserve.csv")
            INDELIBLE_DNDS_FILENAME = SIM_DATA_DIR + "/" + dataset + os.sep + dataset + ".rates.csv"
            EXPECTED_DNDS_FILENAME = REFERENCE_FASTA.replace("consensus.fasta", "dnds.tsv")

            # simulations/data/mydataset/aln/mydataset.reads.bwa.sort.query.sam
            SAM_FILENAME = SIM_DATA_DIR + os.sep + dataset + "/aln/" + dataset + ".reads.bwa.sort.query.sam"
            ORIG_CONSERVE_CSV = SIM_DATA_DIR + os.sep + dataset + "/reads/" + dataset + ".reads.conserve.csv"
            ALN_CONSERVE_CSV = SIM_DATA_DIR + os.sep + dataset + "/aln/" + dataset + ".reads.bwa.conserve.csv"

            ERR_FREE_ALN_CONSENSUS_SAM_FILENAME = SAM_FILENAME.replace(".reads.", ".reads.errFree.")
            ERR_FREE_ORIG_CONSERVE_CSV = ORIG_CONSERVE_CSV.replace(".reads", ".reads.errFree")
            ERR_FREE_ALN_CONSERVE_CSV = ALN_CONSERVE_CSV.replace(".reads", ".reads.errFree")

            outrow=dict(File=SAM_FILENAME, Ref="consensus", OutputPrefix=UMBERJACK_OUT_DIR + os.sep + dataset + os.sep + dataset)

            writer.writerow(outrow)


            sample_ref_out_dir = (UMBERJACK_OUT_DIR + os.sep + os.sep + dataset +  os.sep +"consensus")


            ACTUAL_DNDS_FILENAME = (sample_ref_out_dir + os.sep + dataset + ".dnds.csv")
            # COLLATE_ACT_DNDS_FILENAME = OUT_DIR + os.sep + "collate_dnds.csv"


            # ERR_FREE_OUT_DIR =   SIM_OUT_DIR + os.sep + test_prefix + os.sep + REF + os.sep + "window{}.breadth{}.depth{}.errFree".format(window_size, breadth, depth)
            # ERR_FREE_ACTUAL_DNDS_CSV = ERR_FREE_OUT_DIR + os.sep + 'actual_dnds_by_site.csv'
            # COLLATE_ACT_ERRFREE_DNDS_FILENAME = ERR_FREE_OUT_DIR + os.sep + "collate_dnds.csv"
            #
            # umberjack_html =   OUT_DIR + os.sep + "umberjack_unit_test.html"
            # errfree_umberjack_html =   ERR_FREE_OUT_DIR + os.sep + "umberjack_unit_test.html"
            # if os.path.exists(umberjack_html) and os.path.getsize(umberjack_html) and os.path.exists(errfree_umberjack_html) and os.path.getsize(errfree_umberjack_html):
            #     LOGGER.warn("Not redoing simulations for " + umberjack_html)
            #     return


    do_sliding_window(outdir=UMBERJACK_OUT_DIR, input_csv=umberjack_input_csv,
                      window_size=window_size, window_depth_cutoff=int(min_win_depth), window_breadth_cutoff=min_win_width, min_qual=min_qual,
                      machine_file=machine_file)



    # do_collate(outdir=OUT_DIR, output_csv=COLLATE_ACT_DNDS_FILENAME,
    #            ref_fasta=REFERENCE_FASTA,
    #            full_popn_fasta=FULL_POPN_FASTA,
    #            expected_dnds_filename=EXPECTED_DNDS_FILENAME,
    #            indelible_dnds_filename=INDELIBLE_DNDS_FILENAME,
    #            full_popn_conserve_csv=FULL_POPN_CONSERVE_CSV,
    #            orig_conserve_csv=ORIG_CONSERVE_CSV, aln_conserve_csv=ALN_CONSERVE_CSV)



    LOGGER.debug("Done umberjack and collectdnds for " + UMBERJACK_OUT_DIR)



def get_popn_scales(scale_str):
    """
    From the SIM_ARGS_CSV, gets the genome block mutation scalings per population.

    we allow reproducing the same population at different mutation scalings,
    or generating the genome by concatenating blocks at different mutation scalings to
    Each population is separated by semicolon.
    EG) this indicates that there should be 3 populations, each scaled at 1, 5, 10
    1; 5; 10

    Each genome block is separated by a comma.  Multiple genome blocks per population enclosed by square brackets.
    EG) This indicates that there should be 2 populations,
    one population with concatenated genome blocks scaled at 0.5, 1, 5, 10, 20, 10, 5, 1, 0.5
    and another population with concatenated genome blocks scaled at 0.5 and 5
    [0.5, 1,  5, 10, 20, 10, 5, 1, 0.5]; [0.5, 5]
    :return [[float]] :  for each population returns a list of its genome block scalints
    """
    popn_scales = []
    for popn_scale_str in scale_str.split(";"):
        popn_scale_str = popn_scale_str.lstrip().rstrip()
        popn_scale_str = popn_scale_str.replace("[", "")
        popn_scale_str = popn_scale_str.replace("]", "")
        popn_scale_str = popn_scale_str.replace(" ", "")

        genome_block_scales = [float(block_scale_s) for block_scale_s in popn_scale_str.split(",")]

        popn_scales.append(genome_block_scales)

    return popn_scales


def get_popn_scales_str(scale_str):
    """
    From the SIM_ARGS_CSV, gets the genome block mutation scalings per population.

    we allow reproducing the same population at different mutation scalings,
    or generating the genome by concatenating blocks at different mutation scalings to
    Each population is separated by semicolon.
    EG) this indicates that there should be 3 populations, each scaled at 1, 5, 10
    1; 5; 10

    Each genome block is separated by a comma.  Multiple genome blocks per population enclosed by square brackets.
    EG) This indicates that there should be 2 populations,
    one population with concatenated genome blocks scaled at 0.5, 1, 5, 10, 20, 10, 5, 1, 0.5
    and another population with concatenated genome blocks scaled at 0.5 and 5
    [0.5, 1,  5, 10, 20, 10, 5, 1, 0.5]; [0.5, 5]
    :return [[str]] :  for each population returns a list of its genome block scalings as strings
    """
    popn_scales = []
    for popn_scale_str in scale_str.split(";"):
        popn_scale_str = popn_scale_str.lstrip().rstrip()
        popn_scale_str = popn_scale_str.replace("[", "")
        popn_scale_str = popn_scale_str.replace("]", "")
        popn_scale_str = popn_scale_str.replace(" ", "")

        genome_block_scales = popn_scale_str.split(",")

        popn_scales.append(genome_block_scales)

    return popn_scales


def parse_sim_args_tsv(sim_args_tsv=None):
    """
    Reads in the sim_args.tsv file and gathers all the simulated population groups, and umberjack settings groups.

    :return set, dict:  set of PopnGroup namedtuples, dict of {UmberjackGroup : [PopnGroup's that use those UmberjackGroup settings]}
    """
    # Process all the population creation together. That can most likely be done on my machine whereas umberjack should be done on cluster.
    popn_groups = set()
    umberjack_group_to_popn = dict()

    if not sim_args_tsv:
        sim_args_tsv = SIM_ARGS_TSV

    with open(sim_args_tsv, 'rU') as fh_simargs:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', fh_simargs), delimiter="\t")

        for row in reader:
            seed = row["Seed"]

            # columns that dictate how the population sequencing is formed
            dataset = row["Name"]
            cov_depth = float(row["Cover"])
            fraglen_ave = int(row["FragLenAve"])
            fraglen_std = int(row["FragLenStd"])
            art_profile = row["ART_Profile"]
            popsize = int(row["PopSize"])
            codonsites = int(row["CodonSites"])

            breakpoints_s = row["Breakpoints"] if "Breakpoints" in row else None
            if breakpoints_s:
                breakpoints = [int(x) for x in breakpoints_s.split(",")]
            else:
                breakpoints = []

            num_breakpoints = int(row["NumBreakpoints"])
            selection_rate = float(row["SelectionRate"])
            generations = int(row["Generations"])

            # columns allowed multiple values that dictate how the population sequencing is formed
            scale_str = row["TreeLen"]

            # umberjack settings
            window_size = int(row["WindowSize"])
            min_win_depth = int(row["MinWinDepth"])
            min_qual = int(row["MinQual"])
            min_win_width = float(row["MinWinWidth"])

            # clean the data from errant spaces so when we split the items, we can use the resulting strings in filenames
            scale_str = scale_str.replace(" ", "")

            for genome_block_scales_strs in get_popn_scales_str(scale_str):

                config_file = SIM_DATA_DIR + os.sep + dataset + os.sep + dataset + ".config"

                pgroup = PopnGroup(dataset=dataset,
                                   art_profile=art_profile,
                                   config_file=config_file,
                                   indiv=popsize,
                                   codonsites=codonsites,
                                   cov_depth=cov_depth,
                                   scales=tuple(genome_block_scales_strs),
                                   frag_ave=fraglen_ave,
                                   frag_std=fraglen_std,
                                   breakpoints=tuple(breakpoints),
                                   num_breakpoints=num_breakpoints,
                                   seed=seed,
                                   selection_rate=selection_rate,
                                   generations=generations)

                popn_groups.add(pgroup)

                ugroup = UmberjackGroup(
                    window_size=window_size,
                    min_qual=min_qual,
                    min_win_depth=min_win_depth,
                    min_win_width=min_win_width)

                if ugroup in umberjack_group_to_popn:
                    umberjack_group_to_popn[ugroup].add(pgroup)
                else:
                    umberjack_group_to_popn[ugroup] = set([pgroup])

    return popn_groups, umberjack_group_to_popn


if __name__ == "__main__":

    sim_args_tsv = sys.argv[1]
    concur_mpi = CONCURRENT_MPIRUN
    if len(sys.argv) >= 2:
        concur_mpi = int(sys.argv[2])

    popn_groups, umberjack_group_to_popn = parse_sim_args_tsv(sim_args_tsv=sim_args_tsv)

    thepool =  pool_traceback.LoggingPool(processes=concur_mpi)
    for result in thepool.imap_unordered(gen_sim_data_helper, popn_groups):
            pass

    thepool.close()


    # Let umberjack do the mpi distribution.  Just give it a list of samfiles
    # For each umberjack mpi run, we need to cluster the simulations with the same wnidow breadth and window depth,
    # since each umberjack run allows only 1 breadth & wdepth thresholds, but allows for multiple samfiles.
    thepool =  pool_traceback.LoggingPool(processes=concur_mpi)
    args_itr = []
    for i, (umbjerack_group, pgroups_per_umberjack) in enumerate(umberjack_group_to_popn.iteritems()):
        args_itr.append((i, concur_mpi, umbjerack_group, pgroups_per_umberjack))
    for result in thepool.imap_unordered(process_window_size_helper, args_itr):
        pass
    thepool.close()