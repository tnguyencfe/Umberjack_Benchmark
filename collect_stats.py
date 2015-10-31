import Utility
import csv
import glob
import os
import hyphy.hyphy_handler as hyphy_handler
import Bio.SeqIO as SeqIO
import re
import math
import fnmatch
from collections import namedtuple
import config.settings as settings
import logging
import multiprocessing
import subprocess
import Bio.Phylo as Phylo
from test_topology import TestTopology
settings.setup_logging()
import tempfile
import run_sliding_window_tree as simulator

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = 1

PROCS = 6

def get_tree_len_depth(treefilename):
    """
    Returns tuple of (sum of all branch lengths in tree (excluding root branch), deepest root to tip distance)
    :param treefilename:
    :return: total branch length sum of tree (excluding root branch), deepest root to tip distance
    :rtype: (float, float)
    """
    tree = Phylo.read(treefilename, "newick")

    root_branch_length = tree.clade.branch_length if not tree.clade.branch_length is None else 0  # set to 1.0  for some odd reason
    tree_branch_length  = tree.clade.total_branch_length()  # sum of all branch lengths including root branch length
    unroot_tree_len = tree_branch_length  - root_branch_length
    clade_depths = tree.depths()
    longest_depth = 0.0

    for clade, depth in clade_depths.iteritems():
        if clade.is_terminal():
            if longest_depth < depth:
                longest_depth = depth

    return unroot_tree_len, longest_depth

def get_tree_dist(window_treefile, full_popn_treefile, window_fasta, win_nuc_range):
    """
    finds the distance from expected population tree to actual window tree
    :param window_treefile:
    :param full_popn_treefile:
    :return:
    """


    expected_subsample_treefilename = os.path.basename(full_popn_treefile).replace(".nwk",".resample.missing.{}.nwk".format(win_nuc_range))
    expected_treefile_tmp = tempfile.NamedTemporaryFile(mode="w+", prefix=expected_subsample_treefilename, dir=os.path.dirname(window_treefile), delete=False)
    expected_treefile_tmp.close()
    TestTopology.resample_tree_from_fasta(treefile=full_popn_treefile, out_treefile=expected_treefile_tmp.name, fastafile=window_fasta)


    rf_dist = TestTopology.get_rf_dist(expected_treefile_tmp.name, window_treefile)

    os.remove(expected_treefile_tmp.name)
    return rf_dist



def collect_dnds(output_dir, output_csv_filename, full_popn_fasta, comments=None):
    """
    Collects everything related to dnds into 1 table.  Does not do any aggregation of values.  Useful for debugging.
    :return:
    """
    LOGGER.debug("Collect dnds for " + output_csv_filename)
    with open(output_csv_filename, 'w') as fh_out:
        if comments:
            fh_out.write(comments)

        writer = csv.DictWriter(fh_out, fieldnames=["Window_Start", "Window_End",
                                                    "Reads",  # Max read depth for the window (not necessary for the codon site)
                                                    "CodonSite",  # 1-based codon site
                                                    "CodonDepth",  # Total unambiguous codon (depth) at the codon site
                                                    "AADepth", # Total depth of codons that code unambiguously for 1 AA.
                                                    "Conserve",  # Average per-base fraction of conservation across the codon.  Includes N's and gaps.
                                                    "Entropy",  # Average per-base metric entropy across the codon.  Includes N's and gaps.
                                                    "ConserveTrueBase",  # Average per-base fraction of conservation across the codon.  Excludes N's and gaps
                                                    "EntropyTrueBase",  # Average per-base fraction of entropy across the codon.  Excludes N's and gaps
                                                    "N",  # Observed Nonsynonymous substitutions
                                                    "S",  # Observed Nonsynonymous substitutions
                                                    "EN",  # Expected Nonsynonymous substitutions
                                                    "ES",  # Expected Synonymous substitutions
                                                    "dN", "dS",
                                                    "dN_minus_dS", # dN-dS scaled by the total substitutions at the codon site
                                                    "AmbigBase",  # N nucleotide
                                                    "Pad", # left or right pad gap
                                                    "Err",  # Nucleotide errors within the codon
                                                    "Err_N", # nonsynonymous AA change due to sequence error
                                                    "Err_S",  # synonymous AA change due to sequence error
                                                    "Ambig_N",  # Ambiguous base changes the AA.  Should be always 0
                                                    "Ambig_S", # ambigous base does not change the AA.
                                                    "TreeLen",  # Tree length
                                                    "TreeDepth", # deepest tip to root distance
                                                    "TreeDist",  # distance from actual to expected tree
                                                    ]
                                )
        writer.writeheader()
        for slice_fasta_filename in glob.glob(output_dir + os.sep + "*.fasta"):

            # *.{start bp}_{end bp}.fasta filenames use 1-based nucleotide position numbering
            slice_fasta_fileprefix = slice_fasta_filename.split('.fasta')[0]

            win_nuc_range = slice_fasta_fileprefix.split('.')[-1]

            tree_len = None
            tree_depth = None
            tree_dist = None
            tree_filename = slice_fasta_filename.replace(".fasta", ".tree")
            if os.path.exists(tree_filename):
                # NB:  FastTree tree length in nucleotide substitutions / site.  HyPhy converts trees to codon substitution/site to count codon substitutions along phylogeny
                # Parse the HyPhy dnds tsv to get dN, dS,
                tree_len, tree_depth = get_tree_len_depth(tree_filename)
            else:
                tree_filename = slice_fasta_filename.replace(".fasta", ".nwk")
                if os.path.exists(tree_filename):
                    tree_len, tree_depth = get_tree_len_depth(tree_filename)






            # Window ends at this 1-based nucleotide position with respect to the reference
            if win_nuc_range.find("_") <= 0:  # the full genome msa.fasta file won't have a window range
                continue
            win_start_nuc_pos_1based_wrt_ref, win_end_nuc_pos_1based_wrt_ref = [int(x) for x in win_nuc_range.split('_')]
            # Window starts at this 1-based codon position with respect to the reference
            win_start_codon_1based_wrt_ref = win_start_nuc_pos_1based_wrt_ref/Utility.NUC_PER_CODON + 1

            # TODO:  Hack to remove windows that are less than the expected size.  Remove this once ubmerjack fixed
            actual_win_size = win_end_nuc_pos_1based_wrt_ref - win_start_nuc_pos_1based_wrt_ref + 1

            win_size_match = re.search(pattern=r"window(\d+)", string=os.path.basename(os.path.abspath(os.path.dirname(slice_fasta_filename))))
            if win_size_match:
                desired_win_size = int(win_size_match.group(1))
            else:
                raise ValueError("Window Slice fasta does not obey execpted naming convention for simulated data " + slice_fasta_filename)
            if desired_win_size != actual_win_size:
                LOGGER.warn("Slice fasta is not desired window size=" + str(desired_win_size) + " slicefasta=" + slice_fasta_filename)
                continue

            reads = Utility.get_total_seq_from_fasta(slice_fasta_filename)
            consensus = Utility.Consensus()
            consensus.parse(slice_fasta_filename)

            total_codons = consensus.get_alignment_len()/Utility.NUC_PER_CODON # if the last codon doesn't have enuf chars, then hyphy ignores it

            ns, pad, seq_err, err_aa_change, err_aa_nochange, ambig_aa_change, ambig_aa_nochange = error_by_codonpos(slice_fasta_filename, win_start_nuc_pos_1based_wrt_ref, full_popn_fasta)





            outrow = dict()
            outrow["Window_Start"] = win_start_nuc_pos_1based_wrt_ref
            outrow["Window_End"] = win_end_nuc_pos_1based_wrt_ref
            outrow["Reads"] = reads

            dnds_tsv_filename = slice_fasta_filename.replace(".fasta", ".dnds.tsv")
            if os.path.exists(dnds_tsv_filename):
                full_popn_treefile = full_popn_fasta.replace(".fasta", ".tree")
                if not os.path.exists(full_popn_treefile):
                    full_popn_treefile = full_popn_fasta.replace(".fasta", ".nwk")
                    if not os.path.exists(full_popn_treefile):
                        raise ValueError("Can't find full population tree (with .tree or .nwk suffix)" + full_popn_treefile)

                tree_dist = get_tree_dist(window_treefile=tree_filename, full_popn_treefile=full_popn_treefile, window_fasta=slice_fasta_filename,
                                  win_nuc_range=win_nuc_range)

                with open(dnds_tsv_filename, 'rU') as fh_dnds_tsv:
                    reader = csv.DictReader(fh_dnds_tsv, delimiter='\t')
                    for codonoffset_0based, codon_row in enumerate(reader):    # Every codon site is a row in the *.dnds.tsv file
                        nucoffset_0based = codonoffset_0based*Utility.NUC_PER_CODON
                        outrow["CodonSite"] = win_start_codon_1based_wrt_ref + codonoffset_0based
                        outrow["CodonDepth"] = consensus.get_codon_depth(codon_pos_0based=codonoffset_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                        outrow["AADepth"] = consensus.get_unambig_codon2aa_depth(codon_pos_0based=codonoffset_0based)
                        outrow["Conserve"] = consensus.get_ave_conserve(nucoffset_0based, nucoffset_0based + Utility.NUC_PER_CODON, is_count_ambig=True, is_count_gaps=True, is_count_pad=True)
                        outrow["Entropy"] = consensus.get_ave_metric_entropy(nucoffset_0based, nucoffset_0based + Utility.NUC_PER_CODON, is_count_ambig=True, is_count_gaps=True, is_count_pad=True)
                        outrow["ConserveTrueBase"] = consensus.get_ave_conserve(nucoffset_0based, nucoffset_0based + Utility.NUC_PER_CODON, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                        outrow["EntropyTrueBase"] = consensus.get_ave_metric_entropy(nucoffset_0based, nucoffset_0based + Utility.NUC_PER_CODON, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                        outrow["N"] = float(codon_row[hyphy_handler.HYPHY_TSV_N_COL])
                        outrow["S"] = float(codon_row[hyphy_handler.HYPHY_TSV_S_COL])
                        outrow["ES"] = float(codon_row[hyphy_handler.HYPHY_TSV_EXP_S_COL])
                        outrow["EN"] = float(codon_row[hyphy_handler.HYPHY_TSV_EXP_N_COL])
                        outrow["dN"] = float(codon_row[hyphy_handler.HYPHY_TSV_DN_COL])
                        outrow["dS"] = float(codon_row[hyphy_handler.HYPHY_TSV_DS_COL])
                        outrow["dN_minus_dS"] = float(codon_row[hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL])

                        outrow["AmbigBase"] = ns[codonoffset_0based]
                        outrow["Pad"] = pad[codonoffset_0based]
                        outrow["Err"] = seq_err[codonoffset_0based]
                        outrow["Err_N"] = err_aa_change[codonoffset_0based]
                        outrow["Err_S"] = err_aa_nochange[codonoffset_0based]
                        outrow["Ambig_N"] = ambig_aa_change[codonoffset_0based]
                        outrow["Ambig_S"] = ambig_aa_nochange[codonoffset_0based]
                        outrow["TreeLen"] = tree_len
                        outrow["TreeDepth"] = tree_depth
                        outrow["TreeDist"] = tree_dist
                        writer.writerow(outrow)
            else:
                for codonoffset_0based in range(total_codons):
                    nucoffset_0based = codonoffset_0based*Utility.NUC_PER_CODON
                    outrow["CodonSite"] = win_start_codon_1based_wrt_ref + codonoffset_0based
                    outrow["CodonDepth"] = consensus.get_codon_depth(codon_pos_0based=codonoffset_0based, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                    outrow["AADepth"] = consensus.get_unambig_codon2aa_depth(codon_pos_0based=codonoffset_0based)
                    outrow["Conserve"] = consensus.get_ave_conserve(nucoffset_0based, nucoffset_0based + Utility.NUC_PER_CODON, is_count_ambig=True, is_count_gaps=True, is_count_pad=True)
                    outrow["Entropy"] = consensus.get_ave_metric_entropy(nucoffset_0based, nucoffset_0based + Utility.NUC_PER_CODON, is_count_ambig=True, is_count_gaps=True, is_count_pad=True)
                    outrow["ConserveTrueBase"] = consensus.get_ave_conserve(nucoffset_0based, nucoffset_0based + Utility.NUC_PER_CODON, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                    outrow["EntropyTrueBase"] = consensus.get_ave_metric_entropy(nucoffset_0based, nucoffset_0based + Utility.NUC_PER_CODON, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                    outrow["AmbigBase"] = ns[codonoffset_0based]
                    outrow["Pad"] = pad[codonoffset_0based]
                    outrow["Err"] = seq_err[codonoffset_0based]
                    outrow["Err_N"] = err_aa_change[codonoffset_0based]
                    outrow["Err_S"] = err_aa_nochange[codonoffset_0based]
                    outrow["Ambig_N"] = ambig_aa_change[codonoffset_0based]
                    outrow["Ambig_S"] = ambig_aa_nochange[codonoffset_0based]
                    outrow["TreeLen"] = tree_len
                    outrow["TreeDepth"] = tree_depth
                    writer.writerow(outrow)


def error_by_codonpos(slice_msa_fasta, slice_start_wrt_ref_1based, full_popn_fasta):
    # Collect read error stats per window - codon

    #full_popn_recdict = SeqIO.to_dict(SeqIO.parse("/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov2.indiv1k.codon400.bwa.rand/mixed/small.cov2.indiv1k.codon400.bwa.rand.mixed.fasta", "fasta"))
    full_popn_recdict = SeqIO.to_dict(SeqIO.parse(full_popn_fasta, "fasta"))


    longest_seq = Utility.get_longest_seq_size_from_fasta(slice_msa_fasta)
    total_codons = longest_seq/Utility.NUC_PER_CODON  # don't include the last codon if it isn't fully 3 chars long
    ns = [0] * total_codons
    pad = [0] * total_codons
    seq_err = [0] * total_codons
    err_aa_change = [0] * total_codons
    err_aa_nochange = [0] * total_codons
    ambig_aa_change = [0] * total_codons  # aa change due to ambiguous base.  Should be 0
    ambig_aa_nochange = [0] * total_codons  # no aa change even though ambiguous base


    for record in SeqIO.parse(slice_msa_fasta, "fasta"):
        template, read = record.id.split("_")

        # The read is a subset of the slice
        # first nongap char 0-based nucleotide position wrt slice
        read_start_wrt_slice_0based = re.search(r"[^\-]", str(record.seq)).start()
        # last nongap char 0-based nucleotide position wrt slice
        read_end_wrt_slice_0based = re.search(r"[^\-][\-]*$", str(record.seq)).start()

        # The 0-based codon position with respect to the slice for the first codon that isn't left padded with gap chars
        read_codon_start_wrt_slice_0based = int(math.ceil(float(read_start_wrt_slice_0based)/Utility.NUC_PER_CODON))
        # The 0-based codon position with respect to the slice for the last codon that isn't right padded with gap chars
        read_codon_end_wrt_slice_0based = ((read_end_wrt_slice_0based+1)/Utility.NUC_PER_CODON) - 1

        # NB:  the nonpadded portion of the read might not start on a codon, but the MSA should always start on a codon
        for nuc_pos_wrt_slice_0based in range(0, longest_seq, Utility.NUC_PER_CODON):

            # Ignore codons that aren't fully 3 characters long - this is what hyphy does
            if longest_seq - nuc_pos_wrt_slice_0based < Utility.NUC_PER_CODON:
                continue
            read_codon = str(record.seq[nuc_pos_wrt_slice_0based:nuc_pos_wrt_slice_0based + Utility.NUC_PER_CODON])


            nuc_pos_wrt_ref_0based = slice_start_wrt_ref_1based-1 + nuc_pos_wrt_slice_0based
            template_codon = str(full_popn_recdict[template].seq[nuc_pos_wrt_ref_0based:nuc_pos_wrt_ref_0based + Utility.NUC_PER_CODON])

            codon_pos_wrt_slice_0based = nuc_pos_wrt_slice_0based/Utility.NUC_PER_CODON
            is_codon_has_err = False
            is_codon_has_ambig = False
            for i in range(0, Utility.NUC_PER_CODON):
                if read_codon[i] == "N":
                    is_codon_has_ambig = True
                    ns[codon_pos_wrt_slice_0based] += 1
                elif (read_codon[i] == "-" and (codon_pos_wrt_slice_0based < read_codon_start_wrt_slice_0based or
                                                        codon_pos_wrt_slice_0based > read_codon_end_wrt_slice_0based)):
                    pad[codon_pos_wrt_slice_0based] += 1
                    is_codon_has_ambig = True
                elif read_codon[i] != template_codon[i]:  # sequence error
                    is_codon_has_err = True
                    seq_err[codon_pos_wrt_slice_0based] += 1

            # aa change or no aa-change due to sequence error
            if is_codon_has_err and Utility.CODON2AA.get(read_codon) and Utility.CODON2AA.get(template_codon):
                if Utility.CODON2AA.get(read_codon) == Utility.CODON2AA.get(template_codon):
                    err_aa_nochange[codon_pos_wrt_slice_0based] += 1
                else:
                    err_aa_change[codon_pos_wrt_slice_0based] += 1
            elif is_codon_has_ambig and Utility.CODON2AA.get(read_codon) and Utility.CODON2AA.get(template_codon):
                if Utility.CODON2AA.get(read_codon) == Utility.CODON2AA.get(template_codon):
                    ambig_aa_nochange[codon_pos_wrt_slice_0based] += 1
                else:
                    ambig_aa_change[codon_pos_wrt_slice_0based] += 1


    return ns, pad, seq_err, err_aa_change, err_aa_nochange, ambig_aa_change, ambig_aa_nochange



def make1csv(inferred_dnds_dir, sim_data_dir, output_csv_filename, inferred_collate_dnds_csvs):
    """
    Puts all the collate_dnds, full population csv, expected dnds info into 1 csv for checking what causes inaccurate
    inferred dn/ds.
    Careful - there is about 514MB worth of collatednds csv data
    :return:
    """
    FullPopnCons = namedtuple("FullPopnCons", ["Conservation", "Entropy"])
    FullPopnDnDs = namedtuple("FullPopnDnDs", ["N", "S", "EN", "ES", "dNdS", "dN_minus_dS"])

    LOGGER.debug("Writing all collated inferred, expected dnds to " + output_csv_filename)
    with open(output_csv_filename, 'w') as fh_out:
        writer = csv.DictWriter(fh_out, fieldnames=["Window_Start",
                                                    "Window_End",
                                                    "CodonSite",
                                                    "File",
                                                    "Reads.Act", # max read depth for entire slice
                                                    "UnambigCodonRate.Act", # Total unambiguous codon (depth) at the codon site / max read depth for entire slice
                                                    "AADepth.Act",  # Total codons that code for only 1 amino acid at the codon site
                                                    "PopSize.Act",  # Population size
                                                    "ConserveTrueBase.Act",
                                                    "EntropyTrueBase.Act",  # Average per-base fraction of entropy across the codon site.  Excludes N's and gaps
                                                    "AmbigPadBaseRate.Act",  # Average N or gaps per codon at this site
                                                    "ErrBaseRate.Act",  # Average erroneous bases per codon at this site
                                                    "N.Act", "S.Act",
                                                    #"EN.Act", "ES.Act",
                                                    "dNdS.Act",
                                                    "dN_minus_dS.Act",
                                                    "TreeLen.Act",  # length of window tree in nucleotide subs/site
                                                    "TreeDepth.Act",  # depth of longest branch in nucleotide subs/site
                                                    "TreeDist.Act", # distance from actual to expected tree in Robinson Foulds
                                                    "ConserveTrueBase.Exp",
                                                    "EntropyTrueBase.Exp",
                                                    "N.Exp", "S.Exp",
                                                    #"EN.Exp", "ES.Exp",
                                                    "dNdS.Exp",
                                                    "dN_minus_dS.Exp"
                                                    ])
        writer.writeheader()
        if inferred_collate_dnds_csvs:

            #for dirpath, dirnames, filenames in  os.walk(inferred_dnds_dir):
            # if dirpath.find("window") < 0:
            #     continue
            # for filename in fnmatch.filter(filenames, "collate_dnds.csv"):
            for inferred_collate_dnds_csv in inferred_collate_dnds_csvs:
                dirpath = os.path.dirname(inferred_collate_dnds_csv)
                filename = os.path.basename(inferred_collate_dnds_csv)

                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/small.cov2.indiv3000.codon500/consensus/window200.breadth0.75.depth300.0
                sim_popn_name = os.path.basename(os.path.abspath(dirpath + os.sep + os.pardir + os.sep + os.pardir))
                window_traits = os.path.basename(dirpath)

                LOGGER.debug("sim_popn_name=" + sim_popn_name)


                breadth = None
                depth = None
                indiv = None
                indiv_match = re.search(pattern=r"\.indiv(\d+)\.", string=sim_popn_name)
                if indiv_match:
                    indiv = int(indiv_match.group(1))
                else:
                    raise ValueError("Sim name doesn't obey convention " + sim_popn_name)


                window_matches = re.search(pattern=r"\.breadth(0\.\d+)\.depth(\d+)", string=window_traits)
                if window_matches:
                    breadth = float(window_matches.group(1))
                    depth = float(window_matches.group(2))/indiv
                else:
                    raise ValueError("Window traits dir name doesn't obey convention " + window_traits)

                inferred_collate_dnds_csv = dirpath + os.sep + filename
                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/small.cov1.indiv1000.codon500.window200.breadth0.6.depth100.0/mixed/small.cov1.indiv1000.codon500.window200.breadth0.6.depth100.0.mixed.dnds.tsv
                full_popn_dnds_tsv = sim_data_dir + os.sep + sim_popn_name + os.sep + "mixed" + os.sep + sim_popn_name + ".mixed.dnds.tsv"
                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/small.cov1.indiv1000.codon500.window200.breadth0.6.depth100.0/mixed/small.cov1.indiv1000.codon500.window200.breadth0.6.depth100.0.mixed.conserve.csv
                full_popn_conserve_csv = sim_data_dir + os.sep + sim_popn_name + os.sep + "mixed" + os.sep + sim_popn_name + ".mixed.conserve.csv"
                full_popn_fasta = sim_data_dir + os.sep + sim_popn_name + os.sep + "mixed" + os.sep + sim_popn_name + ".mixed.fasta"
                LOGGER.debug("Merge Sim Inferred collated dnds=" + inferred_collate_dnds_csv)
                LOGGER.debug("Merge Sim expected dnds tsv = " + full_popn_dnds_tsv)
                LOGGER.debug("Merge Sim expected conservation csv = " + full_popn_conserve_csv)
                LOGGER.debug("Merge Sim full popn fasta = " + full_popn_fasta)


                total_indiv = Utility.get_total_seq_from_fasta(full_popn_fasta)

                #NucSite	Conserve	Entropy	NucDepth	CodonDepth

                codonsite_2_full_cons= dict()
                with open(full_popn_conserve_csv, 'rU') as fh_full_cons:
                    full_cons_reader = csv.DictReader(fh_full_cons)

                    total_cons = 0.0
                    total_ent = 0.0
                    for row_idx, row in enumerate(full_cons_reader):
                        nucsite = int(row["NucSite"])
                        total_cons += float(row["Conserve"])
                        total_ent += float(row["Entropy"])
                        if nucsite-1 != row_idx:
                            raise ValueError("Error in nuc site indexing in " + full_popn_conserve_csv)
                        if nucsite % 3 == 0:
                            codonsite = nucsite/3
                            codon_ave_cons = total_cons/3.0
                            codon_ave_ent = total_ent/3.0 / math.log(total_indiv)  # hack because full population conserve.csv are shannon entropy not metric entropy
                            full_popn_cons = FullPopnCons(Conservation=codon_ave_cons, Entropy=codon_ave_ent)
                            codonsite_2_full_cons[codonsite] = full_popn_cons
                            total_cons = 0.0
                            total_ent = 0.0

                # Observed S Changes	Observed NS Changes	E[S Sites]	E[NS Sites]	Observed S. Prop.	P{S}	dS	dN	dN-dS	P{S leq. observed}	P{S geq. observed}	Scaled dN-dS
                codonsite_2_full_dnds= dict()
                with open(full_popn_dnds_tsv, 'rU') as fh_full_dnds:
                    full_dnds_reader = csv.DictReader(fh_full_dnds, delimiter="\t")
                    for row_idx, row in enumerate(full_dnds_reader):
                        codonsite = row_idx+1
                        if float(row[hyphy_handler.HYPHY_TSV_DS_COL]):
                            dnds = float(row[hyphy_handler.HYPHY_TSV_DN_COL])/float(row[hyphy_handler.HYPHY_TSV_DS_COL])
                        else:
                            dnds = None
                        dN_minus_dS = row[hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL]
                        full_popn_dnds = FullPopnDnDs(N=row[hyphy_handler.HYPHY_TSV_N_COL], S=row[hyphy_handler.HYPHY_TSV_S_COL],
                                                      EN=row[hyphy_handler.HYPHY_TSV_EXP_N_COL], ES=row[hyphy_handler.HYPHY_TSV_EXP_S_COL],
                                                      dNdS=dnds, dN_minus_dS=dN_minus_dS)
                        codonsite_2_full_dnds[codonsite] = full_popn_dnds

                if len(codonsite_2_full_dnds.keys()) != len(codonsite_2_full_cons.keys()):
                    raise ValueError("full population dnds does not have same number of codon sites as conservation:",
                                     full_popn_dnds_tsv, ", ", full_popn_conserve_csv)



                #Window_Start, Window_End, Reads, CodonSite, CodonDepth, ConserveTrueBase, EntropyTrueBase, N, S, dN, dS, AmbigBase, Pad, Err
                with open(inferred_collate_dnds_csv, 'rU') as fh_inf_dnds:
                    inf_dnds_reader = csv.DictReader(fh_inf_dnds)
                    for row_idx, row in enumerate(inf_dnds_reader):
                        codonsite = int(row["CodonSite"])
                        outrow = dict()
                        outrow["Window_Start"] = row["Window_Start"]
                        outrow["Window_End"] = row["Window_End"]
                        outrow["CodonSite"] = codonsite
                        outrow["File"] = inferred_collate_dnds_csv
                        outrow["Reads.Act"] = row["Reads"]
                        outrow["UnambigCodonRate.Act"] = float(row["CodonDepth"])/float(row["Reads"])
                        outrow["AADepth.Act"]  = row["AADepth"] if  row.get("AADepth") else None  # TODO: hack since some earlier simulations don't have this value
                        outrow["PopSize.Act"] = total_indiv
                        outrow["ConserveTrueBase.Act"] = row["ConserveTrueBase"]
                        outrow["EntropyTrueBase.Act"] = row["EntropyTrueBase"]
                        outrow["AmbigPadBaseRate.Act"] = (int(row["AmbigBase"]) + int(row["Pad"]))/float(row["Reads"])
                        outrow["ErrBaseRate.Act"] = int(row["Err"])/float(row["Reads"])
                        # If it never made it past FastTree into hyphy, then the substitutions will be empty string
                        if row["N"] != "" and row["S"] != "":
                            outrow["N.Act"] = float(row["N"])
                            outrow["S.Act"] = float(row["S"])
                            #outrow["EN.Act"] = float(row["EN"])
                            #outrow["ES.Act"] = float(row["ES"])
                            if row["dS"] and float(row["dS"]) != 0:
                                outrow["dNdS.Act"] = float(row["dN"])/float(row["dS"])

                            outrow["dN_minus_dS.Act"] = row["dN_minus_dS"]
                        outrow["TreeLen.Act"] = row["TreeLen"]
                        outrow["TreeDepth.Act"] = row["TreeDepth"]
                        outrow["TreeDist.Act"] = row["TreeDist"]


                        if not codonsite_2_full_cons.get(codonsite):
                            raise ValueError("Missing codon site" + str(codonsite) + " in " + full_popn_conserve_csv)
                        outrow["ConserveTrueBase.Exp"] = codonsite_2_full_cons[codonsite].Conservation
                        outrow["EntropyTrueBase.Exp"] = codonsite_2_full_cons[codonsite].Entropy

                        if not codonsite_2_full_dnds.get(codonsite):
                            raise ValueError("Missing codon site" + str(codonsite) + " in " + inferred_collate_dnds_csv)
                        outrow["N.Exp"] = codonsite_2_full_dnds[codonsite].N
                        outrow["S.Exp"] = codonsite_2_full_dnds[codonsite].S
                        #outrow["EN.Exp"] = codonsite_2_full_dnds[codonsite].EN
                        #outrow["ES.Exp"] = codonsite_2_full_dnds[codonsite].ES
                        outrow["dNdS.Exp"] = codonsite_2_full_dnds[codonsite].dNdS
                        outrow["dN_minus_dS.Exp"] = codonsite_2_full_dnds[codonsite].dN_minus_dS

                        writer.writerow(outrow)


def collect_dnds_helper(args):
    """
    Calls collect_dnds() using keyword args
    :param args:
    :return:
    """
    ret_val = -1
    try:
        ret_val = collect_dnds(**args)
    except Exception, e:
        LOGGER.exception("Failure in collect_dnds with args=" + str(args))
        raise e
    return ret_val

def recollect_dnds_all(all_inferred_dnds_dir, sim_data_dir):
    """
    Recollects the collate_dnds.csv info for simulated data.
    Only collects the simulated data specified in the sim_args.tsv file.

    :param all_inferred_dnds_dir:
    :param sim_data_dir:
    :return:
    """
    popn_groups, umberjack_group_to_args = simulator.parse_sim_args_tsv()

    pool = multiprocessing.Pool(PROCS)
    args_itr = []

    i = 0
    for umberjackgroup, popn_groups_per_ugroup in umberjack_group_to_args.iteritems():
        for popn_group in popn_groups_per_ugroup:
            sample_ref_outdir = simulator.get_sample_ref_outdir(umberjackgroup, popn_group)
            #for filename in fnmatch.filter(filenames, "collate_dnds.csv"):
            for filename in fnmatch.filter(filenames, "actual_dnds_by_site.csv"):
                #/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/smTall.cov5.indiv1000.codon500.window350.breadth0.6.depth100.0/consensus/window350/collate_dnds.csv
                sim_name = os.path.basename(os.path.abspath(dirpath + os.sep + os.pardir + os.sep + os.pardir))


                LOGGER.debug("Recollecting sim_name=" + sim_name)
                inferred_collate_dnds_csv = dirpath + os.sep + "collate_dnds.csv"
                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/small.cov1.indiv1000.codon500.window200.breadth0.6.depth100.0/mixed/small.cov1.indiv1000.codon500.window200.breadth0.6.depth100.0.mixed.fasta
                full_popn_fasta = sim_data_dir + os.sep + sim_name + os.sep + "mixed" + os.sep + sim_name + ".mixed.fasta"

                LOGGER.debug("Recollating Inferred collated dnds=" + inferred_collate_dnds_csv)
                LOGGER.debug("Recollating Sim Full Popn Fasta = " + full_popn_fasta)

                # # TODO:  remove me
                # if dirpath.find("/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/small.cov5.indiv1000.codon140.scale50/consensus/window200.breadth0.9.depth10.0.errFree") < 0:
                #     continue
                args_itr.append(dict(output_dir=dirpath, output_csv_filename=inferred_collate_dnds_csv, full_popn_fasta=full_popn_fasta))
                inferred_collate_dnds_csvs.extend([inferred_collate_dnds_csv])
                if i >= 10:
                    break

                i += 1

            if i >= 10:
                break


        for result in pool.imap_unordered(collect_dnds_helper, args_itr, 1):
            pass

        pool.terminate()
        return inferred_collate_dnds_csvs


if __name__ == "__main__":
    SIM_OUT_DIR = os.path.dirname(os.path.realpath(__file__)) + "/simulations/out"
    SIM_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + "/simulations/data"
    OUTPUT_INF_EXP_COLLATE_CSV = SIM_OUT_DIR + os.sep + "collate_all.treedist.csv"
    inferred_collate_dnds_csvs = recollect_dnds_all(all_inferred_dnds_dir=SIM_OUT_DIR, sim_data_dir=SIM_DATA_DIR)
    make1csv(inferred_dnds_dir=SIM_OUT_DIR, sim_data_dir=SIM_DATA_DIR, output_csv_filename=OUTPUT_INF_EXP_COLLATE_CSV, inferred_collate_dnds_csvs=inferred_collate_dnds_csvs)

    # LOGGER.debug("About to generate rhtml")
    # Rscript_wdir =  os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "R")
    # subprocess.check_call(["Rscript", "-e",
    #                        ("library(knitr); " +
    #                         "setwd('{}'); ".format(Rscript_wdir) +
    #                         "spin('find_covar_accuracy_multisample.R', knit=FALSE); " +
    #                         "knit2html('./find_covar_accuracy_multisample.Rmd', stylesheet='./markdown_bigwidth.css')")],
    #                       shell=False, env=os.environ)
    # LOGGER.debug("Done generate rhtml")