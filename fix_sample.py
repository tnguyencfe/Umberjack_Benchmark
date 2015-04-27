"""
Rename small.cov*.indiv*.window*... to just small.cov*.indiv*
"""
import os
import fnmatch
import config.settings as settings
import logging
import shutil

settings.setup_logging()

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = 1


SIM_DATA_DIR =  os.path.dirname(os.path.realpath(__file__)) + os.sep + "simulations/data"
SIM_OUT_DIR =   os.path.dirname(os.path.realpath(__file__)) + os.sep +"simulations/out"

MAPQ_CUTOFF = 20  # alignment quality cutoff
# We set this to 1 because we know that our error rate is good.  We only want breadth to determine sliding windows.
MAX_PROP_N = 1  # maximum proportion of N bases in MSA-aligned sequence
READ_QUAL_CUTOFF = 20   # Phred quality score cutoff [0,40]


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

RENAME_TARGET_DIR = SIM_DATA_DIR + os.sep + "small.cov5.indiv1000.codon500"
TARGET_PATTERN = "small.cov5.indiv1000.codon500.window200.breadth0.6.depth100.0"
RENAME = "small.cov5.indiv1000.codon500"

def recurse_rename(target_dir, pattern, rename):
    """

    :param str  target_dir:  target directory in which all files matching the pattern will be renamed
    :return:
    """

    for dirpath, dirnames, file_basenames in os.walk(target_dir):
        for file_basename in fnmatch.filter(file_basenames, "*" + pattern + "*"):
            old_filename = dirpath + os.sep + file_basename
            new_filename = dirpath + os.sep + file_basename.replace(pattern, rename)
            shutil.move(old_filename, new_filename)
            LOGGER.debug("Renamed " + old_filename + " to " + new_filename )



if __name__ == "__main__":
    recurse_rename(RENAME_TARGET_DIR, TARGET_PATTERN, RENAME)