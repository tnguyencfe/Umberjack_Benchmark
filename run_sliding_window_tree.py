import os
import subprocess
import shutil
import glob
import re
import logging
import random
import sys
import Bio.SeqIO as SeqIO

from pool import pool_traceback
import Utility
import collect_stats
import config.settings as settings


settings.setup_logging()

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = 1


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

def do_sliding_window(outdir, output_csv, samfilename, ref_fasta, expected_dnds_filename, indelible_dnds_filename,
                      window_size, window_depth_cutoff, window_breadth_cutoff):
    ref_len = Utility.get_seq2len(ref_fasta)[REF]

    call = ["mpirun",
            "--machinefile", SIM_OUT_DIR + os.sep + "machine_window" + str(window_size) + ".txt",
            "--output-filename", outdir + os.sep + "logs" + os.sep + os.path.basename(outdir) + ".log",
            "python",
            os.path.dirname(os.path.realpath(__file__)) + os.sep + "../SlidingWindow/umberjack.py",
            "--out_dir", outdir,
            "--map_qual_cutoff", str(MAPQ_CUTOFF),
            "--read_qual_cutoff", str(READ_QUAL_CUTOFF),
            "--max_prop_n", str(MAX_PROP_N),
            "--window_size", str(window_size),
            "--window_slide", str(WINDOW_SLIDE),
            "--window_breadth_cutoff", str(window_breadth_cutoff),
            "--window_depth_cutoff", str(window_depth_cutoff),
            "--start_nucpos", str(START_NUCPOS),
            "--threads_per_window", str(THREADS_PER_WINDOW),
            "--mpi",
            "--mode", "DNDS",
            "--output_csv_filename", output_csv,
            "--sam_filename_list", samfilename_list,
            "--ref", REF
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

    # umberjack.eval_windows_async(ref=REF, sam_filename=samfilename,
    #                                                        out_dir=outdir,
    #                                                        map_qual_cutoff=MAPQ_CUTOFF,
    #                                                        read_qual_cutoff=READ_QUAL_CUTOFF,
    #                                                        max_prop_n=MAX_PROP_N,
    #                                                        start_nucpos=START_NUCPOS,
    #                                                        end_nucpos=ref_len,
    #                                                        window_size=window_size,
    #                                                        window_depth_cutoff=window_depth_cutoff,
    #                                                        window_breadth_cutoff=window_breadth_cutoff,
    #                                                        threads_per_window=THREADS_PER_WINDOW,
    #                                                        concurrent_windows=WINDOW_PROCS,
    #                                                        output_csv_filename=output_csv,
    #                                                        window_slide=WINDOW_SLIDE,
    #                                                        insert=KEEP_INSERTS,
    #                                                        mask_stop_codon=MASK_STOP_CODON,
    #                                                        remove_duplicates=REMOVE_DUPLICATES)

    rconfig_file = os.path.dirname(os.path.realpath(__file__)) + os.sep +"simulations" + os.sep + "R" + os.sep + "umberjack_unit_test.config"
    with open(rconfig_file, 'w') as fh_out_config:
        fh_out_config.write("ACTUAL_DNDS_FILENAME=" + output_csv + "\n")
        fh_out_config.write("EXPECTED_DNDS_FILENAME=" + expected_dnds_filename + "\n")
        fh_out_config.write("EXPECTED_DNDS_START_NUC_POS=" + str(START_NUCPOS) + "\n")
        fh_out_config.write("EXPECTED_DNDS_END_NUC_POS=" + str(ref_len) + "\n")
        fh_out_config.write("INDELIBLE_DNDS_FILENAME=" + indelible_dnds_filename + "\n")

    Rscript_wdir =  os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations" + os.sep + "R")
    subprocess.check_call(["Rscript", "-e",
                           ("library(knitr); " +
                            "setwd('{}'); ".format(Rscript_wdir) +
                            "spin('umberjack_unit_test.R', knit=FALSE); " +
                            "knit2html('./umberjack_unit_test.Rmd', stylesheet='./markdown_bigwidth.css')")],
                          shell=False, env=os.environ)
    shutil.copy(Rscript_wdir + os.sep + "umberjack_unit_test.html",
                outdir + os.sep + "umberjack_unit_test.html")



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


def gen_sim_data(config_file, indiv, codonsites, cov, scales):
    """

    :param str config_file: config file to output to
    :param int indiv:  total individals in population
    :param int codonsites: number of codon sites.  This must be perfectly divisible by the number of scales
    :param int cov:  fold read coverage
    :param [int] scales:  number of different mtuations scalings.  a different population will be generated at each mutation scaling and
    then sites from the different populations will be combined in randomized manner into single population.  But each scaling will be represented equally in the combo population.
    """
    sim_outdir = os.path.dirname(config_file)
    filename_prefix = os.path.basename(config_file).split(".config")[0]
    if not os.path.exists(sim_outdir):
        os.makedirs(sim_outdir)
    if os.path.exists(config_file) and os.path.getsize(config_file):
        LOGGER.warn("Not regenerating " + config_file)
    else:
        seed = random.randint(0, sys.maxint)


        readfrag_ave = 87
        readfrag_std = 125
        LOGGER.debug("Creating simulation config " + config_file + " with seed " + str(seed))
        with open(config_file, 'w') as fh_out:
            fh_out.write("[sim]\n")
            fh_out.write("FILENAME_PREFIX={}\n".format(filename_prefix))
            fh_out.write("NUM_INDIV={}\n".format(indiv))
            fh_out.write("SEED={}\n".format(seed))
            fh_out.write("NUM_CODON_SITES={}\n".format(codonsites))
            fh_out.write("INDELIBLE_BIN_DIR=../../bin/indelible/indelible_1.03/linux_x64\n")
            fh_out.write("INDELIBLE_SCALING_RATES={}\n".format(",".join([str(x) for x in scales])))
            fh_out.write("ART_BIN_DIR = ../../bin/art/art_3.09.15/linux_x64\n")
            fh_out.write("ART_QUAL_PROFILE_TSV1 = ../../bin/art/art_3.09.15/Illumina_profiles/EmpMiSeq250R1.txt\n")
            fh_out.write("ART_QUAL_PROFILE_TSV2 = ../../bin/art/art_3.09.15/Illumina_profiles/EmpMiSeq250R2.txt\n")
            fh_out.write("ART_FOLD_COVER={}\n".format(cov))
            fh_out.write("ART_MEAN_FRAG = {}\n".format(260))
            fh_out.write("ART_STDEV_FRAG = {}\n".format(75))
            fh_out.write("ART_READ_LENGTH = 250\n")
            fh_out.write("ART_INSERT_RATE1 = 0.00045\n")
            fh_out.write("ART_INSERT_RATE2 = 0.00045\n")
            fh_out.write("ART_DEL_RATE1 = 0.00045\n")
            fh_out.write("ART_DEL_RATE2 = 0.00045\n")
            fh_out.write("ART_QUAL_SHIFT1 = 2\n")
            fh_out.write("ART_QUAL_SHIFT2 = 2\n")
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


def process_window_size((test_prefix, indiv, window_size)):
    """
    Umberjack output for the given configurations.
    :return:
    """
     # Run sliding windows and collect stats
    FULL_POPN_FASTA =  SIM_DATA_DIR + os.sep + test_prefix + "/mixed/" + test_prefix + ".mixed.fasta"
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

    for breadth in [0.7, 0.8, 0.9]:
        for depth in [0.01*indiv, 0.1*indiv]:
            OUT_DIR =   SIM_OUT_DIR + os.sep + test_prefix + os.sep + REF + os.sep + "window{}.breadth{}.depth{}".format(window_size, breadth, depth)
            ACTUAL_DNDS_FILENAME = OUT_DIR + os.sep + 'actual_dnds_by_site.csv'
            COLLATE_ACT_DNDS_FILENAME = OUT_DIR + os.sep + "collate_dnds.csv"


            ERR_FREE_OUT_DIR =   SIM_OUT_DIR + os.sep + test_prefix + os.sep + REF + os.sep + "window{}.breadth{}.depth{}.errFree".format(window_size, breadth, depth)
            ERR_FREE_ACTUAL_DNDS_CSV = ERR_FREE_OUT_DIR + os.sep + 'actual_dnds_by_site.csv'
            COLLATE_ACT_ERRFREE_DNDS_FILENAME = ERR_FREE_OUT_DIR + os.sep + "collate_dnds.csv"

            umberjack_html =   OUT_DIR + os.sep + "umberjack_unit_test.html"
            errfree_umberjack_html =   ERR_FREE_OUT_DIR + os.sep + "umberjack_unit_test.html"
            if os.path.exists(umberjack_html) and os.path.getsize(umberjack_html) and os.path.exists(errfree_umberjack_html) and os.path.getsize(errfree_umberjack_html):
                LOGGER.warn("Not redoing simulations for " + umberjack_html)
                return


            # # typical reads
            do_sliding_window(outdir=OUT_DIR, output_csv=ACTUAL_DNDS_FILENAME,
                              samfilename=SAM_FILENAME, ref_fasta=REFERENCE_FASTA,
                              expected_dnds_filename=EXPECTED_DNDS_FILENAME,
                              indelible_dnds_filename=INDELIBLE_DNDS_FILENAME,
                              window_size=window_size, window_depth_cutoff=int(depth), window_breadth_cutoff=breadth)
            # do_collate(outdir=OUT_DIR, output_csv=COLLATE_ACT_DNDS_FILENAME,
            #            ref_fasta=REFERENCE_FASTA,
            #            full_popn_fasta=FULL_POPN_FASTA,
            #            expected_dnds_filename=EXPECTED_DNDS_FILENAME,
            #            indelible_dnds_filename=INDELIBLE_DNDS_FILENAME,
            #            full_popn_conserve_csv=FULL_POPN_CONSERVE_CSV,
            #            orig_conserve_csv=ORIG_CONSERVE_CSV, aln_conserve_csv=ALN_CONSERVE_CSV)


            # #errfree reads
            # do_sliding_window(outdir=ERR_FREE_OUT_DIR, output_csv=ERR_FREE_ACTUAL_DNDS_CSV,
            #                   samfilename=ERR_FREE_ALN_CONSENSUS_SAM_FILENAME, ref_fasta=REFERENCE_FASTA,
            #                   expected_dnds_filename=EXPECTED_DNDS_FILENAME,
            #                   indelible_dnds_filename=INDELIBLE_DNDS_FILENAME,
            #                   window_size=window_size, window_depth_cutoff=int(depth), window_breadth_cutoff=breadth)
            # do_collate(outdir=ERR_FREE_OUT_DIR, output_csv=COLLATE_ACT_ERRFREE_DNDS_FILENAME,
            #            ref_fasta=REFERENCE_FASTA,
            #            full_popn_fasta=FULL_POPN_FASTA,
            #            expected_dnds_filename=EXPECTED_DNDS_FILENAME,
            #            indelible_dnds_filename=INDELIBLE_DNDS_FILENAME,
            #            full_popn_conserve_csv=FULL_POPN_CONSERVE_CSV,
            #            orig_conserve_csv=ERR_FREE_ORIG_CONSERVE_CSV, aln_conserve_csv=ERR_FREE_ALN_CONSERVE_CSV)

            LOGGER.debug("Done umberjack and collectdnds for " + OUT_DIR)



if __name__ == "__main__":

    pool = pool_traceback.LoggingPool(CONCURRENT_MPIRUN)
    TEST_PREFIX_FORMAT = "small.cov{}.indiv{}.codon{}.scale{}"
    for cov in [0.5, 1, 2]:
        for indiv in [1000, 2000]:
            for codonsites in [3000]:
                # 90 days to 10 years, average 1 year apart

                for scales in [[1], [5], [10], [1, 5, 10]]:
                    scales_join = "_".join([str(x) for x in scales])
                    test_prefix = TEST_PREFIX_FORMAT.format(cov, indiv, codonsites, scales_join)
                    config_file = SIM_DATA_DIR + os.sep + test_prefix + os.sep + test_prefix + ".config"
                    LOGGER.debug("Handling simulated config " + config_file)

                    gen_sim_data(config_file, indiv, codonsites, cov, scales)


                    # for result in pool.imap_unordered(process_window_size, [(test_prefix, indiv, 200),
                    #                                                         (test_prefix, indiv, 300),
                    #                                                         (test_prefix, indiv, 350)]):
                    #     pass

    pool.close()