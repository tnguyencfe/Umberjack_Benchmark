import os
import subprocess
import shutil
import glob
import re
import logging
import random
import sys
import Bio.SeqIO as SeqIO
import Bio.Phylo as Phylo
import csv
from pool import pool_traceback
import Utility
import collect_stats
import config.settings as settings
import sam.sam_handler
from test_topology import TestTopology
import tempfile


settings.setup_logging()

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = 1


SIM_ARGS_TSV = "./sim_config/sim_args.tsv"  # a csv that specifies which simulations to run

SIM_DATA_DIR =  os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations/data"
SIM_OUT_DIR =   os.path.dirname(os.path.realpath(__file__)) + os.sep +"simulations/out"

MAPQ_CUTOFF = 20  # alignment quality cutoff
# We set this to 1 because we know that our error rate is good.  We only want breadth to determine sliding windows.
MAX_PROP_N = 1  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]


CONCURRENT_MPIRUN = 3
THREADS_PER_WINDOW = 4
WINDOW_PROCS = 3
START_NUCPOS = 1
SMOOTH_DIST=10
WINDOW_SLIDE = 30
PROCS = 20

MASK_STOP_CODON = True
REMOVE_DUPLICATES = True
KEEP_INSERTS = False

REF = "consensus"

# def do_sliding_window(outdir, output_csv, samfilename_list, expected_dnds_filename, indelible_dnds_filename,
#                       window_size, window_depth_cutoff, window_breadth_cutoff, min_qual):
def do_sliding_window(outdir, output_csv, samfilename_list,
                      window_size, window_depth_cutoff, window_breadth_cutoff, min_qual):
    # ref_len = Utility.get_seq2len(ref_fasta)[REF]


    call = ["mpirun",
            "--machinefile", SIM_OUT_DIR + os.sep + "machine_window" + str(window_size) + ".txt",
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
            "--output_csv_filename", output_csv,
            "--sam_filename_list", samfilename_list
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



def do_collate(outdir, output_csv, ref_fasta, full_popn_fasta,expected_dnds_filename, indelible_dnds_filename,
               full_popn_conserve_csv, orig_conserve_csv, aln_conserve_csv):

    collect_stats.collect_dnds(output_dir=outdir, output_csv_filename=output_csv, full_popn_fasta=full_popn_fasta)

    # Rscript_wdir =  os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "R")
    #
    # rcollate_config_file = Rscript_wdir + os.sep + "aggreg_window.config"
    # with open(rcollate_config_file, 'w') as fh_out_config:
    #     fh_out_config.write("COLLATE_DNDS_FILENAME=" + output_csv + "\n")
    #     fh_out_config.write("EXPECTED_DNDS_FILENAME=" + expected_dnds_filename + "\n")
    #     fh_out_config.write("EXPECTED_DNDS_START_NUC_POS=" + str(START_NUCPOS) + "\n")
    #     fh_out_config.write("EXPECTED_DNDS_END_NUC_POS=" + str(ref_len) + "\n")
    #     fh_out_config.write("INDELIBLE_DNDS_FILENAME=" + indelible_dnds_filename + "\n")
    #     fh_out_config.write("SMOOTH_DIST=" + str(SMOOTH_DIST) + "\n")
    #     fh_out_config.write("FULL_POPN_CONSERVE_CSV=" + full_popn_conserve_csv + "\n")
    #     fh_out_config.write("ORIG_CONSERVE_CSV=" + orig_conserve_csv + "\n")
    #     fh_out_config.write("ALN_CONSERVE_CSV=" + aln_conserve_csv + "\n")

    
    # subprocess.check_call(["Rscript", "-e",
    #                        ("library(knitr); " +
    #                         "setwd('{}'); ".format(Rscript_wdir) +
    #                         "spin('aggreg_window.R', knit=FALSE); " +
    #                         "knit2html('./aggreg_window.Rmd', stylesheet='./markdown_bigwidth.css')")],
    #                       shell=False, env=os.environ)
    # shutil.copy(Rscript_wdir + os.sep + "aggreg_window.html",
    #             outdir + os.sep + "aggreg_window.html")

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
def gen_sim_data(config_file, art_profile, indiv, codonsites, cov_depth, scales,
                 readfrag_ave, readfrag_std,
                 seed=None, **kwargs):
    """

    :param str config_file: config file to output to
    :param int indiv:  total individals in population
    :param int codonsites: number of codon sites.  This must be perfectly divisible by the number of scales
    :param int cov_depth:  fold read coverage
    :param [str] scales:  number of different mtuations scalings.  a different population will be generated at each mutation scaling and
    then sites from the different populations will be combined in randomized manner into single population.  But each scaling will be represented equally in the combo population.
    """

    sim_outdir = os.path.dirname(config_file)
    filename_prefix = os.path.basename(config_file).split(".config")[0]
    if not os.path.exists(sim_outdir):
        os.makedirs(sim_outdir)
    if os.path.exists(config_file) and os.path.getsize(config_file):
        LOGGER.warn("Not regenerating " + config_file)
    else:
        if not seed:
            seed = random.randint(0, sys.maxint)

        LOGGER.debug("Creating simulation config " + config_file + " with seed " + str(seed))
        with open(config_file, 'w') as fh_out:
            fh_out.write("[sim]\n")
            fh_out.write("FILENAME_PREFIX={}\n".format(filename_prefix))
            fh_out.write("NUM_INDIV={}\n".format(indiv))
            fh_out.write("SEED={}\n".format(seed))
            fh_out.write("NUM_CODON_SITES={}\n".format(codonsites))
            fh_out.write("CODONS_PER_BLOCK={}\n".format(codonsites))
            fh_out.write("INDELIBLE_BIN_DIR=../../bin/indelible/indelible_1.03/linux_x64\n")
            fh_out.write("INDELIBLE_SCALING_RATES={}\n".format(",".join([str(x) for x in scales])))
            fh_out.write("ART_BIN_DIR = ../../bin/art/art_3.19.15_adapter/linux_x64\n")
            fh_out.write("ART_QUAL_PROFILE_TSV1 = {}R1.txt\n".format(art_profile))
            fh_out.write("ART_QUAL_PROFILE_TSV2 = {}R2.txt\n".format(art_profile))
            fh_out.write("ART_FOLD_COVER={}\n".format(cov_depth))
            fh_out.write("ART_MEAN_FRAG = {}\n".format(readfrag_ave))
            fh_out.write("ART_STDEV_FRAG = {}\n".format(readfrag_std))
            fh_out.write("ART_READ_LENGTH = 250\n")
            fh_out.write("ART_INSERT_RATE1 = 0.00045\n")
            fh_out.write("ART_INSERT_RATE2 = 0.00045\n")
            fh_out.write("ART_DEL_RATE1 = 0.00045\n")
            fh_out.write("ART_DEL_RATE2 = 0.00045\n")
            fh_out.write("ART_QUAL_SHIFT1 = 0\n")
            fh_out.write("ART_QUAL_SHIFT2 = 0\n")
            fh_out.write("PICARD_BIN_DIR = ../../bin/picard/picard_1.129\n")
            fh_out.write("BWA_BIN_DIR = ../../bin/bwa/bwa_0.7.12/linux_x64\n")
            fh_out.write("PROCS = {}\n".format(PROCS))
            fh_out.write("FASTTREE_EXE = ../../bin/fasttree/fasttree_2.1.7/linux_x64/FastTree\n")
            fh_out.write("HYPHY_EXE = ../../bin/hyphy/hyphy_2.2.3/linux_x64/HYPHYMP\n")
            fh_out.write("HYPHY_BASEPATH = ../../bin/hyphy/hyphy_2.2.3/res/TemplateBatchFiles\n")

    # small.cov5.indiv1000.codon500.window350.breadth0.6.depth100.0/mixed/aln/small.cov5.indiv1000.codon500.window350.breadth0.6.depth100.0.mixed.reads.cov.html
    cov_html = sim_outdir + os.sep + "mixed" + os.sep + "aln" + os.sep + filename_prefix + ".mixed.reads.cov.html"
    if os.path.exists(cov_html) and os.path.getsize(cov_html):
        LOGGER.warn("Not remaking simluated data for " + cov_html)
    else:
        sim_pipeline_exe = os.path.dirname(os.path.realpath(__file__)) + "/simulations/sim_pipeline.py"
        subprocess.check_call(["python", sim_pipeline_exe, config_file])


def swap_window_tips(full_popn_treefile, samfile, samref,
                     outdir, window_size, window_slide, breadth, depth, fraction_tip_swap, seed=None):
    """
    Create all the window fastas and window trees.  Shuffle the tip names of up to fraction_tip_swap tips.
    :param str full_popn_read_fasta:
    :param str full_popn_treefile:
    :param str samfile:
    :param str samref:
    :param str outdir:
    :param int window_size:
    :param int window_slide:
    :param float breadth:
    :param int depth:
    :param float fraction_tip_swap:
    :param int seed:
    :return:
    """
    if seed is None:
        seed = random.randint(0, sys.maxint)
        randomizer = random.Random(seed)

    genome_size = sam.sam_handler.get_reflen(sam_filename=samfile, ref=samref)
    # windows are in 1-based nucleotide positions
    for window_start in range(1, genome_size-window_size+2, window_slide):
        window_end = window_start + window_size - 1

        window_prefix = outdir + os.sep + os.path.basename(samfile).replace(".sam", "") + ".{}_{}".format(window_start, window_end)
        window_fasta = window_prefix + ".fasta"
        total_tips = sam.sam_handler.create_msa_slice_from_sam(sam_filename=samfile, ref=samref,
                                                  out_fasta_filename=window_fasta,
                                                  mapping_cutoff=MAPQ_CUTOFF,
                                                  read_qual_cutoff=READ_QUAL_CUTOFF,
                                                  max_prop_N=1,
                                                  breadth_thresh=breadth, start_pos=0, end_pos=0,
                                                  do_insert_wrt_ref=False, do_mask_stop_codon=True,
                                                  do_remove_dup=True, ref_len=0)


        # To duplicate Umberjack functionality, don't bother making window tree if there aren't enough reads
        if total_tips < depth:
            continue

        # For the reads that made it into the window, make tree that follows the full population tree topology.
        correct_treefile_prefix = window_prefix + ".correct.nwk"
        correct_tree_tmpfile = tempfile.NamedTemporaryFile(mode="w+", delete=False, prefix=correct_treefile_prefix)
        correct_tree_tmpfile.close()

        TestTopology.prune_copies_by_seq(in_treefile=full_popn_treefile, out_treefile=correct_tree_tmpfile.name, fastafile=window_fasta)


        # Randomly select subset of tips to shuffle
        total_tip_swap = round(fraction_tip_swap * total_tips)
        orig_selected_tipnums = randomizer.sample(range(0, total_tips), total_tip_swap)

        # Randomly shuffle selected tips
        shuffled_tipnums =  orig_selected_tipnums.copy()  # make a copy because random.shuffle() shuffles list in place
        randomizer.shuffle(shuffled_tipnums)

        window_tree = Phylo.read(correct_tree_tmpfile.name, "newick")
        tips = window_tree.get_terminals()
        orig_selected_tipnames = [ tips[i].name for i in orig_selected_tipnums]  # keep copy of original names because we will overwrite tip names
        for i, swap_tipnum in enumerate(shuffled_tipnums):
            tips[swap_tipnum].name = orig_selected_tipnames[i].name

        window_treefile = window_prefix + ".nwk"
        Phylo.write(window_tree, window_treefile, "newick")


def process_window_size_helper(kw_namedtuple):
    """
    Helper function to call process_window_size with dict as args
    :param args:
    :return:
    """
    kwargs = kw_namedtuple._asdict()
    return process_window_size(**kwargs)


def process_window_size(test_prefixes, window_size, breadth, depth, min_qual, tip_swap, missing, seed=None, **kwargs):
    """
    Umberjack output for the given configurations.
    :return:
    """
    OUT_DIR =   (SIM_OUT_DIR + os.sep +
                 "window{}.breadth{}.depth{}.qual{}.tipswap{}.missing{}".format(window_size, breadth, depth, min_qual, tip_swap, missing))

    sam_filename_list = OUT_DIR + os.sep + "umberjack_sams.txt"

    with open(sam_filename_list, 'w') as fh_out_samlist:

        for test_prefix in test_prefixes:
            # Run sliding windows and collect stats
            FULL_POPN_FASTA =  SIM_DATA_DIR + os.sep + test_prefix + "/mixed/" + test_prefix + ".mixed.fasta"
            FULL_POPN_TREE =  SIM_DATA_DIR + os.sep + test_prefix + "/mixed/" + test_prefix + ".mixed.nwk"
            REFERENCE_FASTA =  SIM_DATA_DIR + os.sep + test_prefix + "/mixed/" + test_prefix + ".mixed.consensus.fasta"

            FULL_POPN_CONSERVE_CSV = REFERENCE_FASTA.replace(".consensus.fasta", ".conserve.csv")
            INDELIBLE_DNDS_FILENAME = SIM_DATA_DIR + "/" + test_prefix + "/mixed/" + test_prefix + ".mixed.rates.csv"
            EXPECTED_DNDS_FILENAME = REFERENCE_FASTA.replace("consensus.fasta", "dnds.tsv")

            SAM_FILENAME = SIM_DATA_DIR + os.sep + test_prefix + "/mixed/aln/" + test_prefix + ".mixed.reads.consensus.bwa.sort.query.sam"
            ORIG_CONSERVE_CSV = SIM_DATA_DIR + os.sep + test_prefix + "/mixed/reads/" + test_prefix + ".mixed.reads.conserve.csv"
            ALN_CONSERVE_CSV = SIM_DATA_DIR + os.sep + test_prefix + "/mixed/aln/" + test_prefix + ".mixed.reads.consensus.bwa.conserve.csv"

            ERR_FREE_ALN_CONSENSUS_SAM_FILENAME = SAM_FILENAME.replace(".reads.", ".reads.errFree.")
            ERR_FREE_ORIG_CONSERVE_CSV = ORIG_CONSERVE_CSV.replace(".reads", ".reads.errFree")
            ERR_FREE_ALN_CONSERVE_CSV = ALN_CONSERVE_CSV.replace(".reads", ".reads.errFree")

            # simulations/data/small.cov0.5.indiv1000.codon3000.scale1/mixed/aln/small.cov0.5.indiv1000.codon3000.scale1.mixed.reads.consensus.bwa.sort.query.sam
            fh_out_samlist.write(SAM_FILENAME + "\n")

            sample_ref_out_dir = (OUT_DIR + os.sep +
                                  test_prefix + ".mixed.reads.consensus.bwa.sort.query" + os.sep +
                                  "consensus")

            ACTUAL_DNDS_FILENAME = (sample_ref_out_dir + os.sep +
                                   test_prefix + ".mixed.reads.consensus.bwa.sort.query.dnds.csv")
            # COLLATE_ACT_DNDS_FILENAME = OUT_DIR + os.sep + "collate_dnds.csv"

            # If we want to test recombination, we randomly swap tips in the window tree.
            # Create the window fasta and window trees.  Then call umberjack.
            # Expect umberjack to avoid recreating the window trees since they already exist.


            if tip_swap:
                samref = Utility.get_fasta_headers(REFERENCE_FASTA)[0]
                swap_window_tips(full_popn_treefile=FULL_POPN_TREE, samfile=SAM_FILENAME, samref=samref,
                 outdir=sample_ref_out_dir, window_size=window_size, window_slide=WINDOW_SLIDE,
                 breadth=breadth, depth=depth, fraction_tip_swap=tip_swap, seed=seed)



            # ERR_FREE_OUT_DIR =   SIM_OUT_DIR + os.sep + test_prefix + os.sep + REF + os.sep + "window{}.breadth{}.depth{}.errFree".format(window_size, breadth, depth)
            # ERR_FREE_ACTUAL_DNDS_CSV = ERR_FREE_OUT_DIR + os.sep + 'actual_dnds_by_site.csv'
            # COLLATE_ACT_ERRFREE_DNDS_FILENAME = ERR_FREE_OUT_DIR + os.sep + "collate_dnds.csv"
            #
            # umberjack_html =   OUT_DIR + os.sep + "umberjack_unit_test.html"
            # errfree_umberjack_html =   ERR_FREE_OUT_DIR + os.sep + "umberjack_unit_test.html"
            # if os.path.exists(umberjack_html) and os.path.getsize(umberjack_html) and os.path.exists(errfree_umberjack_html) and os.path.getsize(errfree_umberjack_html):
            #     LOGGER.warn("Not redoing simulations for " + umberjack_html)
            #     return


    # # typical reads
    do_sliding_window(outdir=OUT_DIR, output_csv="",
                      samfilename_list=sam_filename_list,
                      window_size=window_size, window_depth_cutoff=int(depth), window_breadth_cutoff=breadth, min_qual=min_qual)



    # do_collate(outdir=OUT_DIR, output_csv=COLLATE_ACT_DNDS_FILENAME,
    #            ref_fasta=REFERENCE_FASTA,
    #            full_popn_fasta=FULL_POPN_FASTA,
    #            expected_dnds_filename=EXPECTED_DNDS_FILENAME,
    #            indelible_dnds_filename=INDELIBLE_DNDS_FILENAME,
    #            full_popn_conserve_csv=FULL_POPN_CONSERVE_CSV,
    #            orig_conserve_csv=ORIG_CONSERVE_CSV, aln_conserve_csv=ALN_CONSERVE_CSV)



    LOGGER.debug("Done umberjack and collectdnds for " + OUT_DIR)



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



if __name__ == "__main__":

    TEST_PREFIX_FORMAT = "{}.cov{}.frag{}_{}.indiv{}.codon{}.scale{}"

    # Process all the population creation together. That can most likely be done on my machine whereas umberjack should be done on cluster.
    with open(SIM_ARGS_TSV, 'rU') as fh_simargs:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', fh_simargs), delimiter="\t")
        #Name	Cover	ReadLenAve	ReadLenStd	PopSize	Scale	ART_Profile	CodonSites	SameSeed	WindowSize	MinWinWidth	MinWinDepth	MinQual	TipSwap	Ns

        sim_args = set()

        from collections import namedtuple

        PopnGroup = namedtuple("PopnGroup", field_names=["test_prefix", "art_profile", "config_file", "indiv", "codonsites",
                                                         "cov_depth", "scales", "readfrag_ave", "readfrag_std", "seed"])

        UmberjackGroup = namedtuple("UmbjerackGroup", field_names=["window_size", "min_qual", "min_win_width", "tip_swap", "missing"])

        umberjack_group_to_args = dict()

        # Simulations are grouped in the CSV for compact display.  Expand the arguments for each separate simulation.
        # Only some columns allow for grouping of multiple simulations:
        # - QualityShift
        # - MinWinWidth
        # - TipSwap
        # - Ns
        for row in reader:
            # columns only allowed one value that dictate how the population sequencing is formed
            prefix = row["Name"]
            cov_depth = row["Cover"]
            fraglen_ave = row["FragLenAve"]
            fraglen_std = row["FragLenStd"]
            art_profile = row["ART_Profile"]
            popsize = row["PopSize"]
            codonsites = row["CodonSites"]
            is_same_seed = row["SameSeed"].lstrip().rstrip() == "1"

            # columns only allowed one value that dicate how umberjack is run
            window_size = row["WindowSize"]
            min_win_depth = row["MinWinDepth"]

            # columns allowed multiple values that dictate how the population sequencing is formed
            scale_str = row["Scale"]

            # columns allowed multiple values that dictate how umberjack is run
            min_quals = row["MinQual"]
            min_win_widths = row["MinWinWidth"]
            tip_swaps = row["TipSwap"]
            ns = row["Ns"]


            # indicates whether we need to keep the same asg_driver.py seed for all simulations in the group
            if is_same_seed:
                seed = random.randint(0, sys.maxint)
            else:
                seed = None


            for genome_block_scales_strs in get_popn_scales_str(scale_str):

                genome_block_scales_str = "_".join(genome_block_scales_strs)
                test_prefix = TEST_PREFIX_FORMAT.format(prefix,
                                                        cov_depth,
                                                        fraglen_ave,
                                                        fraglen_std,
                                                        popsize,
                                                        codonsites,
                                                        genome_block_scales_str)


                for min_qual in min_quals.split(","):
                    config_file = SIM_DATA_DIR + os.sep + test_prefix + os.sep + test_prefix + ".config"

                    for min_win_width in min_win_widths:
                        for tip_swap in tip_swaps:
                            for missing in ns:
                                pgroup = PopnGroup(test_prefix=test_prefix,
                                               art_profile=art_profile,
                                           config_file=config_file,
                                           indiv=popsize,
                                           codonsites=codonsites,
                                           cov_depth=cov_depth,
                                           scales=tuple(genome_block_scales_strs),
                                           readfrag_ave=fraglen_ave,
                                           readfrag_std=fraglen_std,
                                           seed=seed)

                                sim_args.add(pgroup)

                                ugroup = UmberjackGroup(
                                    window_size=window_size,
                                    min_qual=min_qual,
                                    min_win_width=min_win_width,
                                    tip_swap=tip_swap,
                                    missing=missing)

                                if ugroup in umberjack_group_to_args:
                                    umberjack_group_to_args[ugroup].add(pgroup)
                                else:
                                    umberjack_group_to_args[ugroup] = set(pgroup)

    thepool =  pool_traceback.LoggingPool(processes=8)
    for result in thepool.imap_unordered(gen_sim_data_helper, sim_args):
            pass

    thepool.close()


    # Let umberjack do the mpi distribution.  Just give it a list of samfiles
    # For each umberjack mpi run, we need to cluster the simulations with the same wnidow breadth and window depth,
    # since each umberjack run allows only 1 breadth & wdepth thresholds, but allows for multiple samfiles.
    thepool =  pool_traceback.LoggingPool(processes=8)
    for umbjerack_group, group_sim_args in umberjack_group_to_args.iteritems():
        for result in thepool.imap_unordered(process_window_size_helper, group_sim_args):
            pass
    thepool.close()