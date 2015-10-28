# Check what error rates we get for different quality profiles for ART read generator
import subprocess
import os
import logging
import sys
import config.settings as settings
import random
import fnmatch


LOGGER = logging.getLogger(__name__)
settings.setup_logging()

CWD = os.path.dirname(__file__)
# Special ART binary that produces adapter contamination when fragment size is shorter than read length
ADAPTER_ENABLED_ART_EXE = CWD + "/thirdparty/art_v3.19.15/art_src_ChocolateCherryCake_Linux/art_illumina"
# out of box ART exe
ART_EXE = CWD + "/simulations/bin/art/art_3.09.15/linux_x64/art_illumina"


def get_path_str(path, pardir):
        """
        If absolute path, then returns the path as is.
        If relative path, then returns absolute path of concatenated pardir/path
        :param str path:  absolute or relative file or directory path
        :param str pardir: parent directory to concatenate to path if path is relative directory
        :return str: absolute resolved path
        """
        if not os.path.isabs(path):
            return os.path.join(pardir, path)
        else:
            return path


def create_population_fasta():
    """
    Subset of sim_pipeline.py.  Just creates population fasta using INDELible.
    :return:
    """

    OUTDIR =  CWD + "/simulations/data/benchmark_art_profile"  # Output directory for simulated data

    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)

    # Generate Tree
    SEED = random.randint(0, sys.maxint)
    FILENAME_PREFIX = "bench_prof"
    NUM_CODON_SITES = 360
    NUM_INDIV = 1000
    INDELIBLE_SCALING_RATES = 1

    popn_fasta =   (OUTDIR + os.sep + str(INDELIBLE_SCALING_RATES) + os.sep +
                    "{}.{}_TRUE.fasta".format(FILENAME_PREFIX, INDELIBLE_SCALING_RATES))

    LOGGER.info("Creating population fasta " + popn_fasta + " Using seed=" + str(SEED))


    treefile = OUTDIR + os.sep + FILENAME_PREFIX + ".nwk"
    renamed_treefile = OUTDIR + os.sep + FILENAME_PREFIX + ".rename.nwk"
    if os.path.exists(treefile) and os.path.getsize(treefile) and os.path.exists(renamed_treefile) and os.path.getsize(renamed_treefile):
        LOGGER.warn("Not regenerating trees {} and {}".format(treefile, renamed_treefile) )
    else:
        asg_driver_exe = CWD + "/simulations/asg_driver.py"
        asg_driver_cmd = ["python", asg_driver_exe,
                          OUTDIR + os.sep + FILENAME_PREFIX,
                          str(NUM_INDIV),
                          str(SEED)]
        LOGGER.debug("About to execute " + " ".join(asg_driver_cmd))
        subprocess.check_call(asg_driver_cmd, env=os.environ)
        LOGGER.debug("Finished execute ")


        # Relabel tree nodes to more manageable names.  Reformat tree so that indelible can handle it.
        relabel_phylogeny_exe = CWD + "/simulations/relabel_phylogeny.py"
        relabel_phylogeny_cmd = ["python", relabel_phylogeny_exe,
                                 treefile]
        LOGGER.debug("About to execute " + " ".join(relabel_phylogeny_cmd))
        subprocess.check_call(relabel_phylogeny_cmd, env=os.environ)
        LOGGER.debug("Finished execute ")




    # Use Indelible to create population sequences at different scaling factors (ie mutation rates)
    INDELIBLE_BIN_DIR = CWD + "/simulations/bin/indelible/indelible_1.03/linux_x64"

    if os.path.exists(popn_fasta) and os.path.getsize(popn_fasta):
        LOGGER.warn("Not regenerating population fasta " + popn_fasta)
    else:
        batch_indelible_exe = CWD + "/simulations/indelible/batch_indelible.py"
        indelible_cmd = ["python", batch_indelible_exe,
                         renamed_treefile,  # full filepath to tree
                         str(INDELIBLE_SCALING_RATES),
                         str(SEED),  # random seed
                         str(NUM_CODON_SITES), # number of codon sites in genome
                         OUTDIR,  # indelible output file directory
                         FILENAME_PREFIX,  # Indelible output filename prefix
                         INDELIBLE_BIN_DIR]  # indelible bin dir
        LOGGER.debug("About to execute " + " ".join(indelible_cmd))
        subprocess.check_call(indelible_cmd, env=os.environ)
        LOGGER.debug("Finished execute ")


    return popn_fasta




def create_art_reads(art_exe, ref_fasta, profile_tsv1, profile_tsv2,
                            ir, ir2, dr, dr2, qs, qs2, fold_cover,
                            frag_mean, frag_std, output_prefix, seed):
    """
    Calculates the error rate that ART read simulator produces given the ART fastq profile.
    :return:
    """


    ART_READ_LENGTH = 250


    ART_CMD=[art_exe,
         "-na", # don't output alignment file
         "-ef", # create both error-free and normal reads
         "-sam", # create sam alignment
         "-p",  # paired end,
         "-rs", str(seed),
         "-1", profile_tsv1, # 1st read quality  profile
         "-2", profile_tsv2,  # 2nd read quality profile
         "-d", "read", # read id prefix
         "-i",  ref_fasta, # dna reference fasta
         "-ir", str(ir), # 1st read insertion rate
         "-ir2",  str(ir2), # 2nd read insertion rate
         "-dr",  str(dr), # 1st read deletion rate
         "-dr2",  str(dr2), # 2nd read deletion rate
         "-qs", str(qs),  # Bump up the quality scores of every base in 1st mate so that average error rate = 0.006
         "-qs2", str(qs2),  # Bump up the quality scores of every base in 2nd mate so that average error rate = 0.006
         "-l",  str(ART_READ_LENGTH), # length of read
         "-f", str(fold_cover), # fold coverage
         "-m",  str(frag_mean), # mean fragment size
         "-s",  str(frag_std), # std dev fragment size
         "-o",  output_prefix  # output prefix, including directory
         ]


    logfile = output_prefix + ".art.log"
    with open(logfile, 'w') as fh_log:
        LOGGER.debug( "Logging to " + logfile)
        LOGGER.debug( "About to execute " + " ".join(ART_CMD))
        subprocess.check_call(ART_CMD, env=os.environ, stdout=fh_log, stderr=fh_log)


def cmp_err_rates(outcsv, sams):
    """
    Calls seq_err_rate.py to make a single csv of error rates of the different profiles.
    :param outcsv: output csv
    :param [str] sams:  list of fullpaths to sams created by ART using different profiles
    :return:
    """
    cmd = ["python", CWD + "/seq_err_rate.py",
           outcsv,
           str(min(9, len(sams))), # number of processors
           ",".join(sams)]
    LOGGER.debug("About to execute " + " ".join(cmd))
    subprocess.check_call(cmd)


if __name__ == "__main__":
    is_skip_read_gen = False
    if len(sys.argv) > 1 and sys.argv[1] == "-s":
        print "Skipping art read generation"
        is_skip_read_gen = True


    popn_fasta = create_population_fasta()

    profile_sams = []

    # indel rates taken from http://www.nature.com/nbt/journal/v31/n4/fig_tab/nbt.2522_T1.html
    INDEL_RATE = 0.00045

    TRIALS = 3

    reads_dir = os.path.dirname(popn_fasta) +"/reads"

    if not is_skip_read_gen:

        for i in range(0, TRIALS):
            seed = random.randint(0, sys.maxint)

            for args in [dict(profname="longshot",
                                 profile_tsv1=CWD + "/art_profiles/longshot_ART_profile.R1.txt",
                                 profile_tsv2=CWD + "/art_profiles/longshot_ART_profile.R2.txt",
                                 quality_shift = 0,
                                 fold_cover = 9,
                                 frag_mean = 136,
                                 frag_std = 94),
                        dict(profname="longshot",
                                 profile_tsv1=CWD + "/art_profiles/longshot_ART_profile.R1.txt",
                                 profile_tsv2=CWD + "/art_profiles/longshot_ART_profile.R2.txt",
                                 quality_shift = -2,
                                 fold_cover = 9,
                                 frag_mean = 136,
                                 frag_std = 94),
                dict(profname="longshot",
                                 profile_tsv1=CWD + "/art_profiles/longshot_ART_profile.R1.txt",
                                 profile_tsv2=CWD + "/art_profiles/longshot_ART_profile.R2.txt",
                                 quality_shift = -3,
                                 fold_cover = 9,
                                 frag_mean = 136,
                                 frag_std = 94),
                            dict(profname="longshot",
                                 profile_tsv1=CWD + "/art_profiles/longshot_ART_profile.R1.txt",
                                 profile_tsv2=CWD + "/art_profiles/longshot_ART_profile.R2.txt",
                                 quality_shift = 0,
                                 fold_cover = 6,
                                 frag_mean = 300,
                                 frag_std = 50),
                            dict(profname="default",
                                 profile_tsv1 = CWD + "/simulations/bin/art/art_3.11.14/Illumina_profiles/EmpMiSeq250R1.txt",
                                 profile_tsv2 = CWD + "/simulations/bin/art/art_3.11.14/Illumina_profiles/EmpMiSeq250R2.txt",
                                 quality_shift = 0,
                                 fold_cover = 6,
                                 frag_mean = 300,
                                 frag_std = 50),
                            dict(profname="default",
                                 profile_tsv1 = CWD + "/simulations/bin/art/art_3.11.14/Illumina_profiles/EmpMiSeq250R1.txt",
                                 profile_tsv2 = CWD + "/simulations/bin/art/art_3.11.14/Illumina_profiles/EmpMiSeq250R2.txt",
                                 quality_shift = 0,
                                 fold_cover = 9,
                                 frag_mean = 136,
                                 frag_std = 94)]:


                profname = args["profname"]
                profile_tsv1 = args["profile_tsv1"]
                profile_tsv2 = args["profile_tsv2"]
                quality_shift = args["quality_shift"]
                fold_cover = args["fold_cover"]
                frag_mean = args["frag_mean"]
                frag_std = args["frag_std"]

                profile_read_dir =  reads_dir + os.sep + profname
                if not os.path.exists(profile_read_dir):
                    os.makedirs(profile_read_dir)

                art_filename_prefix =  "{}_indel{}_qs{}_cover{}_fragmean{}_fragstd{}_seed{}".format(profname, INDEL_RATE, quality_shift,
                                                                                                    fold_cover, frag_mean, frag_std, seed)
                art_output_prefix = profile_read_dir + os.sep + art_filename_prefix
                create_art_reads(art_exe=ADAPTER_ENABLED_ART_EXE,
                                 ref_fasta=popn_fasta,
                                 profile_tsv1=profile_tsv1, profile_tsv2=profile_tsv2,
                                 ir=INDEL_RATE, ir2=INDEL_RATE, dr=INDEL_RATE, dr2=INDEL_RATE, qs=quality_shift, qs2=quality_shift,
                                 fold_cover=fold_cover,
                                 frag_mean=frag_mean, frag_std=frag_std, output_prefix=art_output_prefix, seed=seed)


    for dirpath, dirname, filenames in os.walk(reads_dir):
        for sam in fnmatch.filter(filenames, "*.sam"):
            if sam.find("errFree") >= 0:
                continue
            profile_sams.extend([dirpath + os.sep + sam])


    cmp_err_rates(outcsv=CWD + "/simulations/data/benchmark_art_profile/seq_error_rate.csv", sams=profile_sams)