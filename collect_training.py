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
from test.test_topology import TestTopology
settings.setup_logging()
import tempfile
import run_sliding_window_tree as simulator
from argparse import ArgumentParser
from test.simulations.SimData import SimData
import scipy.stats


LOGGER = logging.getLogger(__name__)
LOGGER.propagate = 1

PROCS = 10

def get_tree_len_depth(treefilename, polytomy_brlen_thresh=0.0):
    """
    Returns tuple of (sum of all branch lengths in tree (excluding root branch), deepest root to tip distance)
    :param treefilename:
    :param float polytomy_brlen_thresh:  min branch length to be considered no polytomy
    :return: total branch length sum of tree (excluding root branch), deepest root to tip distance, total polytomies
    :rtype: (float, float, int)
    """
    tree = Phylo.read(treefilename, "newick")

    root_branch_length = tree.clade.branch_length if not tree.clade.branch_length is None else 0  # set to 1.0  for some odd reason
    tree_branch_length  = tree.clade.total_branch_length()  # sum of all branch lengths including root branch length
    unroot_tree_len = tree_branch_length  - root_branch_length
    clade_depths = tree.depths()
    longest_depth = 0.0

    total_polytomies = 0
    for clade, depth in clade_depths.iteritems():
        if clade.is_terminal():
            if longest_depth < depth:
                longest_depth = depth
        if not clade.branch_length or clade.branch_length < polytomy_brlen_thresh:
            if clade.is_terminal():
                total_polytomies += 1
            # else:
            #     total_polytomies += len(clade.clades)


    return unroot_tree_len, longest_depth, total_polytomies


def get_break_ratio(sim_data, win_start, win_end):
    """
    :param SimData sim_data: simulated data instance
    :param int win_start:  1-based window nucleotide start
    :param int win_end:  1based window nucleotide end
    :return float: sum of (min bases on one side of breakpoint / bases on other side of breakpoint) across all breakpoint in window
    """
    break_ratio = 0.0
    for (recomb_start, recomb_end) in sim_data.get_recombo_breaks():
        if win_start <= recomb_start <= win_end:  # breakpoint within window
            # breakpoint is the 1-based nucleotide position of start of the strand switch
            left = recomb_start - win_start
            right = win_end - recomb_start + 1.0
            break_ratio += min(left, right)/max(left, right)
    return break_ratio



def count_resolved(subs_tsv):
    """
    Counts the substitutions that were resolved vs observed
    :return dict:
    {int 0based site:  (float total resolved Nonsyn subs, float total resolved syn subs,
                float total observed nonsyn subs, float total observed syn subs)}
    """
    site_to_counts = dict()
    with open(subs_tsv, 'rU') as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        IDX_TOTAL_RES_NS = 0
        IDX_TOTAL_RES_S = 1
        IDX_TOTAL_OBS_NS = 2
        IDX_TOTAL_OBS_S = 3
        for row in reader:
            site = int(row["Site"])
            nonsyn = float(row["ObservedNonsynSubst"])
            syn = float(row["ObservedSynSubst"])
            child_codon = row["ChildCodon"]
            is_ambig = child_codon.find(",") >= 0

            if site not in site_to_counts:
                site_to_counts[site] = [0.0, 0.0, 0.0, 0.0]

            if is_ambig:
                site_to_counts[site][IDX_TOTAL_RES_NS] += nonsyn
            else:
                site_to_counts[site][IDX_TOTAL_OBS_NS] += nonsyn
            if is_ambig:
                site_to_counts[site][IDX_TOTAL_RES_S] += syn
            else:
                site_to_counts[site][IDX_TOTAL_OBS_S] += syn

    return site_to_counts



def collect_dnds(output_dir, output_csv_filename, sim_data_config, comments=None):
    """
    Collects everything related to dnds into 1 table.  Does not do any aggregation of values.  Useful for debugging.
    :return:
    """
    LOGGER.debug("Collect dnds for " + output_csv_filename)
    with open(output_csv_filename, 'w') as fh_out:

        sim_data = SimData(sim_data_config)
        full_popn_fasta = sim_data.get_fasta()
        full_popn_breaks = sim_data.get_recombo_breaks()

        if comments:
            fh_out.write(comments)

        writer = csv.DictWriter(fh_out,
                                fieldnames=[
                                            "Window_Start", "Window_End",
                                            "Reads",  # Max read depth for the window (not necessary for the codon site)
                                            "CodonSite",  # 1-based codon site
                                            "CodonDepth",  # Total unambiguous codon (depth) at the codon site
                                            "AADepth", # Total depth of codons that code unambiguously for 1 AA.
                                            # "ConserveAllCodon",  # Average per-base fraction of conservation across the codon.  Includes N's and gaps.
                                            # "EntropyAllCodon",  # Average per-base metric entropy across the codon.  Includes N's and gaps.
                                            "ConserveCodon",  # Average per-base fraction of conservation across the codon.  Excludes N's and gaps
                                            "EntropyCodon",  # Average per-base fraction of entropy across the codon.  Excludes N's and gaps
                                            "N",  # Observed Nonsynonymous substitutions
                                            "S",  # Observed Nonsynonymous substitutions
                                            "EN",  # Expected Nonsynonymous substitutions
                                            "ES",  # Expected Synonymous substitutions
                                            "dN", "dS",
                                            "dN_minus_dS", # dN-dS scaled by the tree length
                                            "unscaled_dN_minus_dS", # dN-dS
                                            "Ambig",  # N nucleotide
                                            "Pad", # left or right pad gap
                                            "Gap",   # internal gap between true bases on both sides
                                            "Err",  # Nucleotide errors within the codon
                                            "Err_N", # nonsynonymous AA change due to sequence error
                                            "Err_S",  # synonymous AA change due to sequence error
                                            "Ambig_N",  # Ambiguous base changes the AA.  Should be always 0
                                            "Ambig_S", # ambigous base does not change the AA.
                                            "TreeLen",  # Tree length
                                            "TreeDepth", # deepest tip to root distance
                                            "TreeDist",  # distance from actual to expected tree
                                            "Is_Break",  # Whether a strand switch starts on this codon site
                                            "BreakRatio",  # sum across window breakpoints (ratio of bases on either side of breakpoint)
                                            "Polytomy",  # total polytomies in tree
                                            "P_SameCodonFreq",  # log10 probability that sliced codon frequency distro is same as full population distro
                                            # Total substitutions that were resolved vs observed
                                            "ResolvedPerSub"
                                ]
        )
        writer.writeheader()
        for slice_fasta_filename in glob.glob(output_dir + os.sep + "*.*_*.fasta"):

            # don't use hyphy ancestral fasta  or fullgene msa fasta or expected files
            if (slice_fasta_filename.endswith(".anc.fasta") or
                    slice_fasta_filename.endswith(".msa.fasta") or
                        slice_fasta_filename.find("expected") >= 0):
                continue

            # *.{start bp}_{end bp}.fasta filenames use 1-based nucleotide position numbering
            slice_fasta_fileprefix = slice_fasta_filename.split('.fasta')[0]

            win_nuc_range = slice_fasta_fileprefix.split('.')[-1]
            # # Window ends at this 1-based nucleotide position with respect to the reference
            if win_nuc_range.find("_") <= 0:  # the full genome msa.fasta file won't have a window range
                continue
            win_start_nuc_pos_1based_wrt_ref, win_end_nuc_pos_1based_wrt_ref = [int(x) for x in win_nuc_range.split('_')]
            # Window starts at this 1-based codon position with respect to the reference
            win_start_codon_1based_wrt_ref = win_start_nuc_pos_1based_wrt_ref/Utility.NUC_PER_CODON + 1

            break_ratio = get_break_ratio(sim_data=sim_data, win_start=win_start_nuc_pos_1based_wrt_ref, win_end=win_end_nuc_pos_1based_wrt_ref)

            slice_aln = Utility.Consensus()
            slice_aln.parse(slice_fasta_filename)
            codon_width = slice_aln.get_alignment_len()/Utility.NUC_PER_CODON # if the last codon doesn't have enuf chars, then hyphy ignores it

            tree_len = None
            tree_depth = None
            tree_dist = None
            total_polytomies = None
            slice_tree_filename = slice_fasta_fileprefix + ".nwk"
            if os.path.exists(slice_tree_filename):
                # NB:  FastTree tree length in nucleotide substitutions / site.
                # HyPhy converts trees to codon substitution/site to count codon substitutions along phylogeny
                # Parse the HyPhy dnds tsv to get dN, dS,
                polytomy_brlen_thresh = 1.0/(3 * codon_width)  # branch length treshold below which node is considered polytomy
                tree_len, tree_depth, total_polytomies = get_tree_len_depth(slice_tree_filename, polytomy_brlen_thresh=polytomy_brlen_thresh)

                # If there is recombination, there may be multiple trees.
                # Use the full population tree corresponding to slice portion of the genome.

                tree_dist = TestTopology.calc_window_tree_dist(sim_data=sim_data,
                                                               window_fasta=slice_fasta_filename,
                                                               window_treefile=slice_tree_filename,
                                                               win_start=win_start_nuc_pos_1based_wrt_ref,
                                                               win_end=win_end_nuc_pos_1based_wrt_ref)


            (seq_err, err_aa_change, err_aa_nochange,
             ambig_aa_change, ambig_aa_nochange) = error_by_codonpos(slice_fasta_filename,
                                                                     win_start_nuc_pos_1based_wrt_ref,
                                                                     full_popn_fasta)

            full_popn_aln = Utility.Consensus()
            full_popn_aln.parse(full_popn_fasta)

            dnds_tsv_filename = slice_fasta_fileprefix + ".dnds.tsv"
            subs_tsv_filename = slice_fasta_fileprefix + ".subst.tsv"



            fh_dnds_tsv = None
            reader = None
            site_to_subcounts = dict()
            try:
                if os.path.exists(dnds_tsv_filename) and os.path.getsize(dnds_tsv_filename):
                    fh_dnds_tsv = open(dnds_tsv_filename, 'rU')
                    site_to_subcounts = count_resolved(subs_tsv_filename)
                    reader = csv.DictReader(fh_dnds_tsv, delimiter='\t')

                for codonoffset_0based in xrange(codon_width):
                    nucoffset_0based = codonoffset_0based*Utility.NUC_PER_CODON
                    outrow = dict()
                    outrow["Window_Start"] = win_start_nuc_pos_1based_wrt_ref
                    outrow["Window_End"] = win_end_nuc_pos_1based_wrt_ref
                    outrow["Reads"] = slice_aln.get_total_seqs()
                    outrow["CodonSite"] = win_start_codon_1based_wrt_ref + codonoffset_0based
                    outrow["CodonDepth"] = slice_aln.get_codon_depth(codon_pos_0based=codonoffset_0based,
                                                                     is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                    outrow["AADepth"] = slice_aln.get_unambig_codon2aa_depth(codon_pos_0based=codonoffset_0based)
                    outrow["ConserveCodon"] = slice_aln.get_codon_conserve(codonoffset_0based,
                                                                              is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                    outrow["EntropyCodon"] = slice_aln.get_codon_shannon_entropy(codonoffset_0based,
                                                                                        is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                    outrow["Ambig"] = (slice_aln.get_ambig_count(pos_0based=nucoffset_0based) +
                                               slice_aln.get_ambig_count(pos_0based=nucoffset_0based+1) +
                                               slice_aln.get_ambig_count(pos_0based=nucoffset_0based+2))
                    outrow["Pad"] = (slice_aln.get_pad_count(pos_0based=nucoffset_0based) +
                                     slice_aln.get_pad_count(pos_0based=nucoffset_0based+1) +
                                     slice_aln.get_pad_count(pos_0based=nucoffset_0based+2))
                    outrow["Gap"] = (slice_aln.get_gap_count(pos_0based=nucoffset_0based) +
                                     slice_aln.get_gap_count(pos_0based=nucoffset_0based+1) +
                                     slice_aln.get_gap_count(pos_0based=nucoffset_0based+2))
                    outrow["Err"] = seq_err[codonoffset_0based]
                    outrow["Err_N"] = err_aa_change[codonoffset_0based]
                    outrow["Err_S"] = err_aa_nochange[codonoffset_0based]
                    outrow["Ambig_N"] = ambig_aa_change[codonoffset_0based]
                    outrow["Ambig_S"] = ambig_aa_nochange[codonoffset_0based]
                    outrow["TreeLen"] = tree_len
                    outrow["TreeDepth"] = tree_depth
                    outrow["TreeDist"] = tree_dist


                    outrow["Is_Break"] = 0
                    for nuc_strand_start_wrt_ref_base1, nuc_strand_end_wrt_ref_base1 in full_popn_breaks:
                        nuc_pos_wrt_ref_base1 = win_start_nuc_pos_1based_wrt_ref + nucoffset_0based
                        # If there are no recombination breaks, full_popn_breaks still contains the full genome as a contiguous section
                        # Don't consider first position as breakpoint
                        if len(full_popn_breaks) > 1 and nuc_pos_wrt_ref_base1 == nuc_strand_start_wrt_ref_base1 > 1:
                            outrow["Is_Break"] = 1

                    outrow["BreakRatio"] = break_ratio
                    outrow["Polytomy"] = total_polytomies

                    #  log-likelihood ratio test that codon count distributions are similar between window and full population
                    full_popn_codon_freq = full_popn_aln.get_codon_freq(codon_pos_0based=win_start_codon_1based_wrt_ref + codonoffset_0based -1,
                                                                        is_count_pad=False, is_count_gaps=False, is_count_ambig=False)
                    slice_codon_freq = slice_aln.get_codon_freq(codon_pos_0based=codonoffset_0based,
                                             is_count_pad=False, is_count_gaps=False, is_count_ambig=False)


                    pval_same = cmp_freq_distro(full_popn_codon_freq, slice_codon_freq)
                    outrow["P_SameCodonFreq"] = pval_same

                    if reader:
                        resolved_ns, resolved_s, obs_ns, obs_s = site_to_subcounts[codonoffset_0based]
                        total_subs = resolved_ns +  resolved_s +  obs_ns +  obs_s
                        if total_subs:
                            outrow["ResolvedPerSub"] = (resolved_ns + resolved_s) / float(total_subs)
                        else:
                            outrow["ResolvedPerSub"] = 0

                        dnds_info = reader.next()  # Every codon site is a row in the *.dnds.tsv file
                        if codonoffset_0based != int(dnds_info["Site"]):
                            # dnds tsv specified the codon site in 0-based coordinates in Site field wrt Slice
                            raise ValueError("Inconsistent site numbering " + str(codonoffset_0based) + " in " + dnds_tsv_filename)

                        outrow["N"] = dnds_info[hyphy_handler.HYPHY_TSV_N_COL]
                        outrow["S"] = dnds_info[hyphy_handler.HYPHY_TSV_S_COL]
                        outrow["ES"] = dnds_info[hyphy_handler.HYPHY_TSV_EXP_S_COL]
                        outrow["EN"] = dnds_info[hyphy_handler.HYPHY_TSV_EXP_N_COL]
                        outrow["dN"] = dnds_info[hyphy_handler.HYPHY_TSV_DN_COL]
                        outrow["dS"] = dnds_info[hyphy_handler.HYPHY_TSV_DS_COL]
                        outrow["dN_minus_dS"] = dnds_info[hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL]
                        outrow["unscaled_dN_minus_dS"] = dnds_info[hyphy_handler.HYPHY_TSV_DN_MINUS_DS_COL]

                        if  abs(total_subs -
                                (float(dnds_info[hyphy_handler.HYPHY_TSV_N_COL]) + float(dnds_info[hyphy_handler.HYPHY_TSV_S_COL]))) > 1e-2:
                            raise ValueError("Inconsitent total subs at 0-based site " +
                                             str(codonoffset_0based) + " wrt  " + subs_tsv_filename + " and " + dnds_tsv_filename +
                                             " " + str(total_subs) + " " + str(float(dnds_info[hyphy_handler.HYPHY_TSV_N_COL]) + float(dnds_info[hyphy_handler.HYPHY_TSV_S_COL])))

                    writer.writerow(outrow)

                if reader:
                    try:
                        dnds_info = reader.next()
                        if dnds_info and len(dnds_info) > 0:
                            raise ValueError("dnds TSV has more codons than expected " + dnds_tsv_filename)
                    except StopIteration:  # We want the reader to have no more rows
                        pass

            finally:
                if fh_dnds_tsv and not fh_dnds_tsv.closed:
                    fh_dnds_tsv.close()


def cmp_freq_distro(obs_freq_distro, exp_freq_distro, is_scale=True):
    """
    Does a G-test to check if the frequency distributions are the same.
    See http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.stats.chi2_contingency.html

    If one distribution is missing symbols form the other distribution, then sets the frequency for that symbol to zero.

    If one of the distros contains all zero counts, then returns a pvalue of  None

    :param dict obs_freq_distro:  {str symbol:  int count}
    :param dict exp_freq_distro:   {str symbol:  int count}
    :param bool is_scale:  whether to scale the observed to the expected frequency counts
    :return float:  probability that the frequency distributions are the same
    """
    symbols = set(obs_freq_distro.keys() + exp_freq_distro.keys())
    obs_freq_distro_arr = [obs_freq_distro[symbol] if symbol in obs_freq_distro else 0 for symbol in symbols]
    exp_freq_distro_arr = [exp_freq_distro[symbol] if symbol in exp_freq_distro else 0 for symbol in symbols]

    if is_scale:
        scale = sum(exp_freq_distro_arr) / float(sum(obs_freq_distro_arr))
        for i, obs_freq in enumerate(obs_freq_distro_arr):
            obs_freq_distro_arr[i] = obs_freq * scale

    if obs_freq_distro_arr.count(0) == len(obs_freq_distro_arr) or exp_freq_distro_arr.count(0) == len (exp_freq_distro_arr):
        pval = None
    else:
        _, pval, _, _ = scipy.stats.chi2_contingency(observed=[obs_freq_distro_arr, exp_freq_distro_arr], lambda_="log-likelihood")

    return pval


def error_by_codonpos(slice_msa_fasta, slice_start_wrt_ref_1based, full_popn_fasta):
    """
    Collect read error stats per window - codon
    If the codon translates to multiple AA, then that is considered an amino acid change
    :param slice_msa_fasta:
    :param slice_start_wrt_ref_1based:
    :param full_popn_fasta:
    :return:
    """


    full_popn_recdict = SeqIO.to_dict(SeqIO.parse(full_popn_fasta, "fasta"))


    longest_seq = Utility.get_longest_seq_size_from_fasta(slice_msa_fasta)
    total_codons = longest_seq/Utility.NUC_PER_CODON  # don't include the last codon if it isn't fully 3 chars long
    seq_err = [0] * total_codons
    err_aa_change = [0] * total_codons
    err_aa_nochange = [0] * total_codons
    ambig_aa_change = [0] * total_codons  # aa change due to ambiguous base
    ambig_aa_nochange = [0] * total_codons  # no aa change even though ambiguous base


    for record in SeqIO.parse(slice_msa_fasta, "fasta"):
        template, read = record.id.split("_")

        # NB:  the nonpadded portion of the read might not start on a codon, but the MSA should always start on a codon
        for nuc_pos_wrt_slice_0based in range(0, longest_seq, Utility.NUC_PER_CODON):

            # Ignore codons that aren't fully 3 characters long - this is what hyphy does
            if longest_seq - nuc_pos_wrt_slice_0based < Utility.NUC_PER_CODON:
                continue
            read_codon = str(record.seq[nuc_pos_wrt_slice_0based:nuc_pos_wrt_slice_0based + Utility.NUC_PER_CODON])


            nuc_pos_wrt_ref_0based = slice_start_wrt_ref_1based-1 + nuc_pos_wrt_slice_0based
            template_codon = str(full_popn_recdict[template].seq[nuc_pos_wrt_ref_0based:nuc_pos_wrt_ref_0based + Utility.NUC_PER_CODON])

            if not Utility.CODON2AA.get(template_codon):
                raise ValueError("Template Codon should translate unambiguously to an amino acid! templatecodon" + template_codon)

            codon_pos_wrt_slice_0based = nuc_pos_wrt_slice_0based/Utility.NUC_PER_CODON
            is_codon_has_err = False
            is_codon_has_ambig = False
            for i in range(0, Utility.NUC_PER_CODON):
                if read_codon[i] == "N" or read_codon[i] == "-":
                    is_codon_has_ambig = True
                elif read_codon[i] != template_codon[i]:  # sequence error
                    is_codon_has_err = True
                    seq_err[codon_pos_wrt_slice_0based] += 1


            # If there are both errors and ambiguous/pad bases, then find out the effect of each alone.
            if is_codon_has_ambig and is_codon_has_err:
                # Test if ambiguous base alone causes AA change by correcting erroneous bases
                corrected_read_codon = template_codon
                for i, read_base in enumerate(read_codon):
                    if read_base == "N" or read_base == "-":  # keep ambiguous bases
                        corrected_read_codon = corrected_read_codon[0:i]  + "N" + corrected_read_codon[i+1:]

                if Utility.CODON2AA.get(corrected_read_codon) == Utility.CODON2AA.get(template_codon):
                    ambig_aa_nochange[codon_pos_wrt_slice_0based] += 1
                else:
                    ambig_aa_change[codon_pos_wrt_slice_0based] += 1

                # Test if erroneous bases alone causes AA change by correcting ambiguous bases
                corrected_read_codon = template_codon
                for i, read_base in enumerate(read_codon):
                    # keep erroneous bases
                    if read_base != template_codon[i] and read_base != "N" and read_base != "-":
                        corrected_read_codon = corrected_read_codon[0:i]  + read_base + corrected_read_codon[i+1:]

                if Utility.CODON2AA.get(corrected_read_codon) == Utility.CODON2AA.get(template_codon):
                    err_aa_nochange[codon_pos_wrt_slice_0based] += 1
                else:
                    err_aa_change[codon_pos_wrt_slice_0based] += 1

            # aa change or no aa-change due to sequence error
            # CODON2AA doesn't accept codons with gaps.  Convert gaps to Ns.

            elif is_codon_has_err:
                if Utility.CODON2AA.get(read_codon.replace("-", "N")) == Utility.CODON2AA.get(template_codon):
                    err_aa_nochange[codon_pos_wrt_slice_0based] += 1
                else:
                    err_aa_change[codon_pos_wrt_slice_0based] += 1

            # aa change or no aa change due to ambiguous base
            elif is_codon_has_ambig:
                if Utility.CODON2AA.get(read_codon.replace("-", "N")) == Utility.CODON2AA.get(template_codon):
                    ambig_aa_nochange[codon_pos_wrt_slice_0based] += 1
                else:
                    ambig_aa_change[codon_pos_wrt_slice_0based] += 1


    return seq_err, err_aa_change, err_aa_nochange, ambig_aa_change, ambig_aa_nochange


def read_codon_csv(csv_file, codon_site_field="CodonSite", is_base0=False, delimiter=","):
    """
    Read in an entire csv as a list as dict of dicts.
    :param str codon_site_field: the column name that describes the codon site
    :param bool is_base0:  whether the codon sites are 0based or 1based
    :return dict: {int 1-based codonsite: csv row dict}
    """
    site_to_row = {}
    with open(csv_file, 'rU') as fh:
        #EG) CodonSite	ConserveCodon	Entropy	NucDepth	CodonDepth
        reader = csv.DictReader(fh, delimiter=delimiter)
        for row_idx, row in enumerate(reader):
            site = int(row[codon_site_field])
            if is_base0:
                site_1based = site + 1
            else:
                site_1based = site
            site_to_row[site_1based] = row

    return site_to_row


def make1csv(output_csv_filename, sim_args_tsv):
    """
    Puts all the collate_dnds, full population csv, expected dnds info into 1 csv for checking what causes inaccurate
    inferred dn/ds.
    Careful - there is about 514MB worth of collatednds csv data
    :return:
    """
    LOGGER.debug("Writing all collated inferred, expected dnds to " + output_csv_filename)
    with open(output_csv_filename, 'w') as fh_out:
        writer = csv.DictWriter(fh_out, fieldnames=["Window_Start",
                                                    "Window_End",
                                                    "CodonSite",
                                                    "File",
                                                    "Is_Break",  # whether the site is a recombinant breakpoint (start of new strand)
                                                    "BreakRatio.Act",  # sum across breakpoints (ratio of bases on either side of breakpoint)
                                                    "Reads.Act", # max read depth for entire slice
                                                    # Total unambiguous codon (depth) at the codon site / max read depth for entire slice
                                                    # This number will be in [0, 1]
                                                    "UnambigCodonRate.Act",
                                                    "AADepth.Act",  # Total codons that code for only 1 amino acid at the codon site
                                                    "PopSize.Act",  # Population size
                                                    "ConserveCodon.Act",
                                                    "EntropyCodon.Act",  # Excludes codons with N's and gaps
                                                     # Average N, gaps, or pads per codon at this site.
                                                     # This is different from UnambigCodonRate.Act because this tells us whether there are
                                                     # lots of codons that are only partial, indicating poor reads or alignment,
                                                     # versus completely uncovered by reads.
                                                     # This nubmer will be in [0, 3]
                                                    "UnknownPerCodon.Act",
                                                    "ErrPerCodon.Act",  # Average erroneous bases per codon at this site
                                                    "N.Act", "S.Act",
                                                    "EN.Act", "ES.Act",
                                                    "dNdS.Act",
                                                    "dN_minus_dS.Act",
                                                    "TreeLen.Act",  # length of window tree in nucleotide subs/site
                                                    "TreeDepth.Act",  # depth of longest branch in nucleotide subs/site
                                                    "Polytomy.Act",  # total polytomies in tree
                                                    #log10 pvalue of log likelihood ratio test (aka Gtest) that codon count distros same as fullpopn
                                                    "P_SameCodonFreq.Act",
                                                    # Fraction of substitutions that were resolved as opposed to observed
                                                    "ResolvedPerSub.Act",
                                                    # distance from actual to expected tree in Robinson Foulds-branch lengths /reads
                                                    "TreeDistPerRead.Act",
                                                    "ConserveCodon.Exp",
                                                    "EntropyCodon.Exp",
                                                    "N.Exp", "S.Exp",
                                                    "EN.Exp", "ES.Exp",
                                                    "dNdS.Exp",
                                                    "dN_minus_dS.Exp"
                                                    ])

        writer.writeheader()


        popn_groups, umberjack_group_to_args = simulator.parse_sim_args_tsv(sim_args_tsv)
        for umberjackgroup, popn_groups_per_ugroup in umberjack_group_to_args.iteritems():
            for popn_group in popn_groups_per_ugroup:
                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/simdatasetname
                sim_popn_name = popn_group.dataset
                sim_data = SimData(popn_group.config_file)
                sim_data_dir = sim_data.sim_data_dir



                #/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/window350.breadth0.6.depth100.0.qual20/simdatasetname/consensus/collate_dnds.csv
                sam_ref_outdir = simulator.get_umberjack_outdir_from_simargs_tsv(popn_group=popn_group, umberjack_group=umberjackgroup)
                inferred_collate_dnds_csv = sam_ref_outdir + os.sep + "collate_dnds.csv"

                LOGGER.debug("Merge sim_name=" + sim_popn_name + " collatednds=" + inferred_collate_dnds_csv)


                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/simdatasetname/subs/simdatasetname.dnds.tsv
                full_popn_dnds_tsv = sim_data_dir + os.sep + "subs" + os.sep + sim_popn_name + ".dnds.tsv"
                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/simdatasetname/fullpopn/simdatasetname.conserve.csv
                full_popn_conserve_csv = sim_data_dir + os.sep + "fullpopn" + os.sep + sim_popn_name + "_TRUE.conserve.csv"


                total_indiv = popn_group.indiv
                total_codon_sites = popn_group.codonsites

                #CodonSite	ConserveCodon	Entropy	NucDepth	CodonDepth
                codonsite_2_full_cons = read_codon_csv(csv_file=full_popn_conserve_csv, codon_site_field="CodonSite", is_base0=False)

                # File,Window_Start,Window_End,Reads,CodonSite,CodonDepth,AADepth,ConserveAllCodon,EntropyAllCodon,ConserveCodon,EntropyCodon,N,S,EN,ES,dN,dS,dN_minus_dS,Ambig,Pad,Err,Err_N,Err_S,Ambig_N,Ambig_S,TreeLen,T

                # Site	Observed S Changes	Observed NS Changes	E[S Sites]	E[NS Sites]	dS	dN	dN-dS	Scaled dN-dS
                codonsite_2_full_dnds= read_codon_csv(csv_file=full_popn_dnds_tsv, codon_site_field="Site", is_base0=True, delimiter="\t")

                if (len(codonsite_2_full_dnds.keys()) != len(codonsite_2_full_cons.keys()) or
                            len(codonsite_2_full_dnds.keys()) != total_codon_sites or
                            len(codonsite_2_full_cons.keys()) != total_codon_sites):
                    raise ValueError("full population dnds does not have same number of codon sites as conservation:",
                                     full_popn_dnds_tsv, ", ", full_popn_conserve_csv)

                # Don't read entire umberjack collate_dnds.csv into memory cuz most likely huge
                # File,Window_Start,Window_End,Reads,CodonSite,CodonDepth,AADepth,ConserveAllCodon,EntropyAllCodon,ConserveCodon,EntropyCodon,N,S,EN,ES,dN,dS,dN_minus_dS,Ambig,Pad,Err,Err_N,Err_S,Ambig_N,Ambig_S,TreeLen,...
                with open(inferred_collate_dnds_csv, 'rU') as fh_inf_dnds:
                    inf_dnds_reader = csv.DictReader(fh_inf_dnds)
                    for row_idx, row in enumerate(inf_dnds_reader):
                        codonsite = int(row["CodonSite"])
                        outrow = dict()
                        outrow["Window_Start"] = row["Window_Start"]
                        outrow["Window_End"] = row["Window_End"]
                        outrow["CodonSite"] = codonsite
                        outrow["File"] = sim_data.name # inferred_collate_dnds_csv
                        outrow["Reads.Act"] = row["Reads"]
                        outrow["UnambigCodonRate.Act"] = float(row["CodonDepth"])/float(row["Reads"])
                        outrow["AADepth.Act"]  = row["AADepth"] if  row.get("AADepth") else None  # TODO: hack since some earlier simulations don't have this value
                        outrow["PopSize.Act"] = total_indiv
                        outrow["ConserveCodon.Act"] = row["ConserveCodon"]
                        outrow["EntropyCodon.Act"] = row["EntropyCodon"]
                        outrow["UnknownPerCodon.Act"] = (int(row["Ambig"]) + int(row["Pad"]) + int(row["Gap"]))/float(row["Reads"])
                        outrow["ErrPerCodon.Act"] = int(row["Err"])/float(row["Reads"])
                        # If it never made it past FastTree into hyphy, then the substitutions will be empty string
                        if row["N"] != "" and row["S"] != "":
                            outrow["N.Act"] = float(row["N"])
                            outrow["S.Act"] = float(row["S"])
                            outrow["EN.Act"] = float(row["EN"])
                            outrow["ES.Act"] = float(row["ES"])
                            if row["dS"] and float(row["dS"]) != 0:
                                outrow["dNdS.Act"] = float(row["dN"])/float(row["dS"])

                            outrow["dN_minus_dS.Act"] = row["dN_minus_dS"]
                        outrow["TreeLen.Act"] = row["TreeLen"]
                        outrow["TreeDepth.Act"] = row["TreeDepth"]
                        if row["TreeDist"]:
                            outrow["TreeDistPerRead.Act"] = float(row["TreeDist"])/float(row["Reads"])
                        outrow["Is_Break"] = row["Is_Break"]
                        outrow["BreakRatio.Act"] = row["BreakRatio"]
                        outrow["Polytomy.Act"] = row["Polytomy"]
                        outrow["P_SameCodonFreq.Act"] = row["P_SameCodonFreq"]
                        outrow["ResolvedPerSub.Act"] = row["ResolvedPerSub"]


                        if not codonsite_2_full_cons.get(codonsite):
                            raise ValueError("Missing codon site" + str(codonsite) + " in " + full_popn_conserve_csv)
                        outrow["ConserveCodon.Exp"] = codonsite_2_full_cons[codonsite]["ConserveCodon"]
                        outrow["EntropyCodon.Exp"] = codonsite_2_full_cons[codonsite]["EntropyCodon"]

                        if not codonsite_2_full_dnds.get(codonsite):
                            raise ValueError("Missing codon site" + str(codonsite) + " in " + inferred_collate_dnds_csv)

                        outrow["N.Exp"] = codonsite_2_full_dnds[codonsite][hyphy_handler.HYPHY_TSV_N_COL]
                        outrow["S.Exp"] = codonsite_2_full_dnds[codonsite][hyphy_handler.HYPHY_TSV_S_COL]
                        outrow["EN.Exp"] = codonsite_2_full_dnds[codonsite][hyphy_handler.HYPHY_TSV_EXP_N_COL]
                        outrow["ES.Exp"] = codonsite_2_full_dnds[codonsite][hyphy_handler.HYPHY_TSV_EXP_S_COL]

                        if (codonsite_2_full_dnds[codonsite][hyphy_handler.HYPHY_TSV_S_COL] and
                                    float(codonsite_2_full_dnds[codonsite][hyphy_handler.HYPHY_TSV_S_COL]) != 0):
                            outrow["dNdS.Exp"] = (float(codonsite_2_full_dnds[codonsite][hyphy_handler.HYPHY_TSV_DN_COL])/
                                                  float(codonsite_2_full_dnds[codonsite][hyphy_handler.HYPHY_TSV_DS_COL]))

                        outrow["dN_minus_dS.Exp"] = codonsite_2_full_dnds[codonsite][hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL]

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

def recollect_dnds_from_tsv(sim_args_tsv):
    """
    Recollects the collate_dnds.csv info for simulated data.
    Only collects the simulated data specified in the sim_args.tsv file.

    :param str sim_args_tsv:  filepath to the TSV containing the simulation dataset arguments and umberjack arguments
    :return:
    """
    popn_groups, umberjack_group_to_args = simulator.parse_sim_args_tsv(sim_args_tsv)


    args_itr = []
    for umberjackgroup, popn_groups_per_ugroup in umberjack_group_to_args.iteritems():
        for popn_group in popn_groups_per_ugroup:
            # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/simdatasetname
            sim_name = popn_group.dataset
            sim_data_config = popn_group.config_file

            #/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/window350.breadth0.6.depth100.0.qual20/simdatasetname/collate_dnds.csv
            sam_ref_outdir = simulator.get_umberjack_outdir_from_simargs_tsv(umberjack_group=umberjackgroup, popn_group=popn_group)
            inferred_collate_dnds_csv = sam_ref_outdir + os.sep + "collate_dnds.csv"
            LOGGER.debug("Recollecting sim_data_config=" + sim_data_config + " Inferred collated dnds=" + inferred_collate_dnds_csv)

            args_itr.append(dict(output_dir=sam_ref_outdir, output_csv_filename=inferred_collate_dnds_csv,
                                 sim_data_config=sim_data_config))

    total = 0
    thepool = multiprocessing.Pool(PROCS)
    for result in thepool.imap_unordered(collect_dnds_helper, args_itr, 1):
        total += 1

    thepool.terminate()
    LOGGER.debug("Collected " + str(total) + " stats from (simulated dataset, umberjack settings) combinations")


def recollect_dnds(all_inferred_dnds_dir, sim_config_file):
    """
    Recollects the collate_dnds.csv info for simulated data
    given the simulated data directory and umberjack output directory.
    Assumes that the simulated data and umberjack output directory represent the same dataset.

    :param all_inferred_dnds_dir:
    :param str sim_config_file:  filepath to simulated dataset config file passed to sim_pipeline.py to create it
    :return:
    """

    LOGGER.debug("Recollecting dataset=" + sim_config_file + " with umberjack output " + all_inferred_dnds_dir)
    sim_data = SimData(sim_config_file)
    full_popn_fasta = sim_data.get_fasta()

    inferred_collate_dnds_csv = all_inferred_dnds_dir + os.sep + "collate_dnds.csv"
    LOGGER.debug("Recollating Inferred collated dnds=" + inferred_collate_dnds_csv)
    LOGGER.debug("Recollating Sim Full Popn Fasta = " + full_popn_fasta)

    collect_dnds(output_dir=all_inferred_dnds_dir, output_csv_filename=inferred_collate_dnds_csv, sim_data_config=sim_config_file)

    return inferred_collate_dnds_csv


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("-c", help="simulated data config")
    parser.add_argument("-o", help="simulated data umberjack output directory")
    parser.add_argument("-t", help="simulated data tsv file")
    parser.add_argument("-f", help="output CSV concatenated from dnds, conservation stats from individual datasets and their umberjack outputs")


    args = parser.parse_args()

    SIM_OUT_DIR = os.path.dirname(os.path.realpath(__file__)) + "/simulations/out"
    SIM_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + "/simulations/data"
    OUTPUT_INF_EXP_COLLATE_CSV = SIM_OUT_DIR + os.sep + "collate_all.treedist.csv"

    if args.t:
        recollect_dnds_from_tsv(sim_args_tsv=args.t)
    else:
        inferred_collate_dnds_csvs = recollect_dnds(all_inferred_dnds_dir=args.o, sim_config_file=args.c)

    if args.f:
        output_concat_stat_csv = args.f
    else:
        output_concat_stat_csv = OUTPUT_INF_EXP_COLLATE_CSV
    make1csv(output_csv_filename=output_concat_stat_csv, sim_args_tsv=args.t)

    # LOGGER.debug("About to generate rhtml")
    # Rscript_wdir =  os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + "R")
    # subprocess.check_call(["Rscript", "-e",
    #                        ("library(knitr); " +
    #                         "setwd('{}'); ".format(Rscript_wdir) +
    #                         "spin('find_covar_accuracy_multisample.R', knit=FALSE); " +
    #                         "knit2html('./find_covar_accuracy_multisample.Rmd', stylesheet='./markdown_bigwidth.css')")],
    #                       shell=False, env=os.environ)
    # LOGGER.debug("Done generate rhtml")