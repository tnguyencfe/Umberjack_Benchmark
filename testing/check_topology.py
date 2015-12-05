"""
one-off script to check why increasing ave robinson foulds distance between window tree and true tree
increases the accuracy of umberjack  in find_covar_accuracy_multisample.R
"""


import run_sliding_window_tree
import os
from test.test_topology import TestTopology
from test.simulations.SimData import  SimData
import fasttree.fasttree_handler
import Utility
import Bio.Phylo
import csv
import hyphy.hyphy_handler
import logging
import config.settings
import collect_training
import copy
import subprocess
import shutil
import sam.sam_handler

config_dict = copy.deepcopy(config.settings.DEFAULT_LOG_CONFIG_DICT)
config_dict["root"]["level"] = "ERROR"
config.settings.setup_logging(config_file="", config_dict=config_dict)

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.ERROR)

OUTPUT_COLLATE_SAMEFULL_CSV = run_sliding_window_tree.SIM_OUT_DIR + os.sep + "collate_all.treedist.test.recombo.same.full.csv"
OUTPUT_COLLATE_UMBERJACK_CSV = run_sliding_window_tree.SIM_OUT_DIR + os.sep + "collate_all.treedist.test.recombo.umberjack.csv"
OUTPUT_COLLATE_ERRFREE_UMBERJACK_CSV = run_sliding_window_tree.SIM_OUT_DIR + os.sep + "collate_all.treedist.test.recombo.errfree.umberjack.csv"
OUTPUT_COLLATE_CSV = run_sliding_window_tree.SIM_OUT_DIR + os.sep + "collate_all.treedist.test.recombo.csv"
SIM_ARGS_TSV = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) + os.sep + os.pardir  + os.sep +
                               "sim_config" + os.sep + "sim_args.test.recombo.tsv")

ART_PROFILE = "../../../art_profiles/longshot_ART_profile."
INDIV = 1000
NUCSITES = 300
COVERAGE = 2
TREELENS = [20]
FRAG_AVE = 375
FRAG_STD = 50
SELECTION_RATE = 0.01
GENERATIONS = 5000
PROCS = 6
DATASET_OUT_DIR = run_sliding_window_tree.SIM_OUT_DIR + os.sep + "window300.breadth0.875.depth10.qual20"

def ave_window_umberjack_dnds_sq_diff(sim_data, window_dnds_tsv, window_start):
    """
    :param test.simulations.SimData.SimData sim_data: instance
    :param str window_dnds_tsv: filepath to window dnds csv
    :param int window_start:  1-based nucleotide window start
    :return float:  average (square distance from umberjack dn-ds and true dn-ds)
    """
    expected_dnds_tsv = sim_data.get_dnds()  # 0based codons
    total_diff = 0.0
    total_sites = 0.0
    with open(window_dnds_tsv, 'rU') as fh_actual, open(expected_dnds_tsv, 'rU') as fh_exp:
        reader_exp = csv.DictReader(fh_exp, delimiter="\t")
        reader_act = csv.DictReader(fh_actual, delimiter="\t")
        for row_idx, row_exp in enumerate(reader_exp):
            exp_site_base0 = int(row_exp["Site"])

            row_act = reader_act.next()
            act_site_base0 = int(row_act["Site"]) + window_start - 1
            if exp_site_base0 != act_site_base0:
                raise ValueError("Expected and actual out of sync")

            dnds_exp = float(row_exp[hyphy.hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL])
            dnds_act = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL])

            total_diff += (dnds_exp - dnds_act) ** 2
            total_sites += 1

    return total_diff / total_sites


def ave_window_umberjack_dnds_abs_diff(sim_data, window_dnds_tsv, window_start):
    """
    :param test.simulations.SimData.SimData sim_data: instance
    :param str window_dnds_tsv: filepath to window dnds csv
    :param int window_start:  1-based nucleotide window start
    :return float:  average (square distance from umberjack dn-ds and true dn-ds)
    """
    expected_dnds_tsv = sim_data.get_dnds()  # 0based codons
    total_diff = 0.0
    total_sites = 0.0
    with open(window_dnds_tsv, 'rU') as fh_actual, open(expected_dnds_tsv, 'rU') as fh_exp:
        reader_exp = csv.DictReader(fh_exp, delimiter="\t")
        reader_act = csv.DictReader(fh_actual, delimiter="\t")
        for row_idx, row_exp in enumerate(reader_exp):
            exp_site_base0 = int(row_exp["Site"])

            row_act = reader_act.next()
            act_site_base0 = int(row_act["Site"]) + window_start - 1
            if exp_site_base0 != act_site_base0:
                raise ValueError("Expected and actual out of sync")

            dnds_exp = float(row_exp[hyphy.hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL])
            dnds_act = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL])

            total_diff += abs(dnds_exp - dnds_act)
            total_sites += 1

    return total_diff / total_sites


def spit_results(dataset, window_fasta, window_treefile, window_dnds_tsv, window_start, window_end):
    """
    :return float, float:  robinson fould distance, ave squared difference between actual and expected dnds
    """
    conf = run_sliding_window_tree.SIM_DATA_DIR + os.sep + dataset + os.sep +  dataset + ".config"
    sim_data = SimData(config_file=conf)
    sim_data_treefiles = sim_data.get_recombo_trees()

    diff = TestTopology.calc_window_tree_dist(sim_data=sim_data,
                                       window_fasta=window_fasta,
                                       window_treefile=window_treefile,
                                       win_start=window_start, win_end=window_end)
    print "dataset=" + dataset + " rf=" + str(diff)
    print "\twindow tree=" + window_treefile

    tree = Bio.Phylo.read(window_treefile, "newick")
    print "\t window tree length = " + str(tree.total_branch_length())

    for recombo_treefile in sim_data_treefiles:
        tree = Bio.Phylo.read(recombo_treefile, "newick")
        print "\t recombo tree " + recombo_treefile + " length= " + str(tree.total_branch_length())

    sq_acc = ave_window_umberjack_dnds_sq_diff(sim_data=sim_data, window_dnds_tsv=window_dnds_tsv, window_start=window_start)
    print "\t ave sq difference from true dnds =" + str(sq_acc)

    abs_acc = ave_window_umberjack_dnds_abs_diff(sim_data=sim_data, window_dnds_tsv=window_dnds_tsv, window_start=window_start)
    print "\t ave abs difference from true dnds =" + str(abs_acc)

    return diff, sq_acc, abs_acc


def get_sim_data(dataset):
    conf = run_sliding_window_tree.SIM_DATA_DIR + os.sep + dataset + os.sep +  dataset + ".config"
    sim_data = SimData(config_file=conf)
    return sim_data



def test_window_same_full(dataset):
    """
    Window is same size as full population sequence.
    Window contains exactly same sequences as full population.
    :return:
    """
    sim_data = get_sim_data(dataset)

    # Instead of using umberjack window fasta, we use the full population fasta
    window_fasta = sim_data.get_fasta()
    # Instead of using umberjack window tree, we use a fasttree made from full population fasta
    # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/window300.breadth0.875.depth10.qual20/TestRecomboNone
    umberjack_output_prefix = DATASET_OUT_DIR + os.sep + dataset + os.sep + dataset + ".repro"
    window_treefile = umberjack_output_prefix + ".nwk"
    fasttree.fasttree_handler.make_tree(fasta_fname=window_fasta,
                                        out_tree_fname=window_treefile, threads=PROCS, debug=True)


    hyphy.hyphy_handler.calc_dnds(codon_fasta_filename=window_fasta,
                                  tree_filename=window_treefile,
                                  hyphy_filename_prefix=umberjack_output_prefix,
                                  threads=PROCS, debug=True)
    window_dnds_tsv = umberjack_output_prefix + ".dnds.tsv"
    window_start = 1
    genome_size = Utility.get_len_1st_seq(window_fasta)
    spit_results(dataset=dataset, window_fasta=window_fasta, window_treefile=window_treefile,
                 window_dnds_tsv=window_dnds_tsv,
                 window_start=window_start, window_end=genome_size)


def test_window_reads(dataset, window_start, window_end):
    """
    Windows taken directly from umberjack output on reads.
    :return:
    """

    umberjack_output_prefix = DATASET_OUT_DIR + os.sep + dataset + os.sep + dataset + ".{}_{}".format(window_start, window_end)
    window_fasta = umberjack_output_prefix + ".fasta"
    window_treefile = umberjack_output_prefix + ".nwk"
    window_dnds_tsv = umberjack_output_prefix + ".dnds.tsv"

    spit_results(dataset=dataset, window_fasta=window_fasta, window_treefile=window_treefile,
                 window_dnds_tsv=window_dnds_tsv,
                 window_start=window_start, window_end=window_end)



def test_errFree_window_reads(dataset, window_start, window_end):
    """
    Windows artificially created from error-free ART reads.  Tests the effect of windowing compared to no windowing
    (ie tests effect of missing data).
    Removes effect of inaccurate alignment and sequencing error.
    :return:
    """
    umberjack_output_prefix = DATASET_OUT_DIR + os.sep + dataset + os.sep + "errfree" +  os.sep + dataset + ".repro.errfree.{}_{}".format(window_start, window_end)
    window_fasta = umberjack_output_prefix + ".fasta"
    window_treefile = umberjack_output_prefix + ".nwk"
    window_dnds_tsv = umberjack_output_prefix + ".dnds.tsv"

    if not os.path.exists(os.path.dirname(window_fasta)):
        os.makedirs(os.path.dirname(window_fasta))

    # Instead of using umberjack window fasta, slice the error free ART generated reads and then artificially reproduce umberjack.
    # We can use the error free MSA generated the by sim_pipeline.py and feed that into umberjack, but we don't so that
    # we can see the effect of missing data without effect of any unknown umberjack bugs.

    # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/TestRecomboNone/reads/TestRecomboNone.reads.errFree.sort.query.sam
    sim_data = get_sim_data(dataset)
    art_errfree_sam = sim_data.get_art_sam(is_errfree=True)
    sam.sam_handler.create_msa_slice_from_sam(sam_filename=art_errfree_sam, ref=None,
                                              out_fasta_filename=window_fasta,
                                              mapping_cutoff=0,  # doesn't matter, since err free sam
                                              read_qual_cutoff=0,  # doesn't matter since err free sam
                                              max_prop_N=1.0, # doesn't matter since err free sam
                                              # TODO:  take these from the sim_args_tsv file instead of hardcoding
                              breadth_thresh=0.875, start_pos=window_start, end_pos=window_end,
                              do_insert_wrt_ref=False, do_mask_stop_codon=True,
                              do_remove_dup=False, ref_len=0)

    fasttree.fasttree_handler.make_tree(fasta_fname=window_fasta,
                                        out_tree_fname=window_treefile, threads=PROCS, debug=True)


    hyphy.hyphy_handler.calc_dnds(codon_fasta_filename=window_fasta,
                                  tree_filename=window_treefile,
                                  hyphy_filename_prefix=umberjack_output_prefix,
                                  threads=PROCS, debug=True)


    spit_results(dataset=dataset, window_fasta=window_fasta, window_treefile=window_treefile,
                 window_dnds_tsv=window_dnds_tsv,
                 window_start=window_start, window_end=window_end)


def make1csv_same_full(output_csv_filename, sim_args_tsv):
    """
    Puts all the collate_dnds, full population csv, expected dnds info into 1 csv for checking what causes inaccurate
    inferred dn/ds.
    Instead of putting read sequencing umberjack results in the csv, instead it puts umberjack results
    from the full population.  Ie the windows are made on the sequences from the full population instead of the ART reads.

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
                                                    "UnambigCodonRate.Act", # Total unambiguous codon (depth) at the codon site / max read depth for entire slice
                                                    "AADepth.Act",  # Total codons that code for only 1 amino acid at the codon site
                                                    "PopSize.Act",  # Population size
                                                    "ConserveCodon.Act",
                                                    "EntropyCodon.Act",  # Excludes codons with N's and gaps
                                                    "UnknownPerCodon.Act",  # Average N or gaps per codon at this site
                                                    "ErrPerCodon.Act",  # Average erroneous bases per codon at this site
                                                    "N.Act", "S.Act",
                                                    "EN.Act", "ES.Act",
                                                    "dNdS.Act",
                                                    "dN_minus_dS.Act",
                                                    "TreeLen.Act",  # length of window tree in nucleotide subs/site
                                                    "TreeDepth.Act",  # depth of longest branch in nucleotide subs/site
                                                    "Polytomy.Act",
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


        popn_groups, umberjack_group_to_args = run_sliding_window_tree.parse_sim_args_tsv(sim_args_tsv)
        for umberjackgroup, popn_groups_per_ugroup in umberjack_group_to_args.iteritems():
            for popn_group in popn_groups_per_ugroup:
                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/simdatasetname
                sim_popn_name = popn_group.dataset
                sim_data = SimData(popn_group.config_file)
                sim_data_dir = sim_data.sim_data_dir


                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/simdatasetname/subs/simdatasetname.dnds.tsv
                full_popn_dnds_tsv = sim_data_dir + os.sep + "subs" + os.sep + sim_popn_name + ".dnds.tsv"
                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/simdatasetname/fullpopn/simdatasetname.conserve.csv
                full_popn_conserve_csv = sim_data_dir + os.sep + "fullpopn" + os.sep + sim_popn_name + "_TRUE.conserve.csv"


                # Instead of using umberjack window fasta made from reads, we use window made from full population fasta
                window_fasta = sim_data.get_fasta()
                umberjack_output_prefix = DATASET_OUT_DIR + os.sep + popn_group.dataset + os.sep + popn_group.dataset
                window_treefile = umberjack_output_prefix + ".repro.nwk"
                window_dnds_tsv = umberjack_output_prefix + ".repro.dnds.tsv"
                window_start = 1
                window_end = Utility.get_len_1st_seq(window_fasta)

                LOGGER.debug("Merge sim_name=" + sim_popn_name + " full popn window dnds tsv =" + window_dnds_tsv)

                total_indiv = popn_group.indiv
                total_codon_sites = popn_group.codonsites

                #CodonSite	ConserveCodon	Entropy	NucDepth	CodonDepth
                codonsite_2_full_cons = collect_training.read_codon_csv(csv_file=full_popn_conserve_csv, codon_site_field="CodonSite", is_base0=False)

                # File,Window_Start,Window_End,Reads,CodonSite,CodonDepth,AADepth,ConserveAllCodon,EntropyAllCodon,ConserveCodon,EntropyCodon,N,S,EN,ES,dN,dS,dN_minus_dS,Ambig,Pad,Err,Err_N,Err_S,Ambig_N,Ambig_S,TreeLen,T

                # Site	Observed S Changes	Observed NS Changes	E[S Sites]	E[NS Sites]	dS	dN	dN-dS	Scaled dN-dS
                codonsite_2_full_dnds= collect_training.read_codon_csv(csv_file=full_popn_dnds_tsv, codon_site_field="Site", is_base0=True, delimiter="\t")

                if (len(codonsite_2_full_dnds.keys()) != len(codonsite_2_full_cons.keys()) or
                            len(codonsite_2_full_dnds.keys()) != total_codon_sites or
                            len(codonsite_2_full_cons.keys()) != total_codon_sites):
                    raise ValueError("full population dnds does not have same number of codon sites as conservation:",
                                     full_popn_dnds_tsv, ", ", full_popn_conserve_csv)

                aln = Utility.Consensus()
                aln.parse(msa_fasta_filename=window_fasta)
                window_reads = aln.get_total_seqs()
                window_tree_dist = TestTopology.calc_window_tree_dist(sim_data=sim_data,
                                       window_fasta=window_fasta,
                                       window_treefile=window_treefile,
                                       win_start=window_start, win_end=window_end)
                full_popn_breaks = sim_data.get_recombo_breaks()

                break_ratio = collect_training.get_break_ratio(sim_data=sim_data, win_start=window_start, win_end=window_end)

                polytomy_brlen_thresh = 1.0/(3 * total_codon_sites)  # branch length treshold below which node is considered polytomy
                window_treelen, window_treedepth, total_polytomies = collect_training.get_tree_len_depth(window_fasta, polytomy_brlen_thresh=polytomy_brlen_thresh)
                with open(window_dnds_tsv, 'rU') as fh_actual:
                    reader_act = csv.DictReader(fh_actual, delimiter="\t")
                    for row_idx, row_act in enumerate(reader_act):
                        act_codonsite_offset_base0 = int(row_act["Site"])
                        act_codonsite_base0 = act_codonsite_offset_base0 + window_start - 1


                        act_nucsite_offset_base0 = act_codonsite_base0 * 3
                        codonsite_base1 = act_codonsite_base0 + 1
                        unambig_codon_depth = aln.get_codon_depth(codon_pos_0based=act_codonsite_offset_base0, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                        outrow = dict()
                        outrow["Window_Start"] = window_start
                        outrow["Window_End"] = window_end
                        outrow["CodonSite"] = codonsite_base1
                        outrow["File"] = "Full_" + sim_data.name
                        outrow["Reads.Act"] = window_reads
                        outrow["UnambigCodonRate.Act"] = float(unambig_codon_depth)/window_reads
                        outrow["AADepth.Act"]  = aln.get_unambig_codon2aa_depth(codon_pos_0based=act_codonsite_offset_base0)
                        outrow["PopSize.Act"] = total_indiv
                        outrow["ConserveCodon.Act"] = aln.get_codon_conserve(codon_pos_0based=act_codonsite_offset_base0,
                                                                             is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                        outrow["EntropyCodon.Act"] = aln.get_codon_shannon_entropy(codon_pos_0based=act_codonsite_offset_base0,
                                                                             is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                        outrow["UnknownPerCodon.Act"] = float(aln.get_gap_count(pos_0based=act_nucsite_offset_base0) +
                                                              aln.get_ambig_count(pos_0based=act_nucsite_offset_base0) +
                                                              aln.get_pad_count(pos_0based=act_nucsite_offset_base0)) / window_reads
                        outrow["ErrPerCodon.Act"] = 0
                        # If it never made it past FastTree into hyphy, then the substitutions will be empty string
                        if row_act[hyphy.hyphy_handler.HYPHY_TSV_N_COL] and row_act[hyphy.hyphy_handler.HYPHY_TSV_S_COL]:
                            outrow["N.Act"] = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_N_COL])
                            outrow["S.Act"] = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_S_COL])
                            outrow["EN.Act"] = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_EXP_N_COL])
                            outrow["ES.Act"] = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_EXP_S_COL])
                            if row_act["dS"] and float(row_act[hyphy.hyphy_handler.HYPHY_TSV_S_COL]) != 0:
                                outrow["dNdS.Act"] = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_DN_COL])/float(row_act[hyphy.hyphy_handler.HYPHY_TSV_DS_COL])
                            outrow["dN_minus_dS.Act"] = row_act[hyphy.hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL]
                        outrow["TreeLen.Act"] = window_treelen
                        outrow["TreeDepth.Act"] = window_treedepth
                        outrow["TreeDistPerRead.Act"] = float(window_tree_dist)/window_reads

                        outrow["Is_Break"] = 0
                        for nuc_strand_start_wrt_ref_base1, nuc_strand_end_wrt_ref_base1 in full_popn_breaks:
                            nuc_pos_wrt_ref_base1 = window_start + act_nucsite_offset_base0
                            # If there are no recombination breaks, full_popn_breaks still contains the full genome as a contiguous section
                            # Don't consider first position as breakpoint
                            if len(full_popn_breaks) > 1 and nuc_pos_wrt_ref_base1 == nuc_strand_start_wrt_ref_base1 > 1:
                                outrow["Is_Break"] = 1

                        outrow["BreakRatio.Act"] = break_ratio
                        outrow["Polytomy.Act"] = total_polytomies

                        if not codonsite_2_full_cons.get(codonsite_base1):
                            raise ValueError("Missing codon site" + str(codonsite_base1) + " in " + full_popn_conserve_csv)
                        outrow["ConserveCodon.Exp"] = codonsite_2_full_cons[codonsite_base1]["ConserveCodon"]
                        outrow["EntropyCodon.Exp"] = codonsite_2_full_cons[codonsite_base1]["EntropyCodon"]

                        if not codonsite_2_full_dnds.get(codonsite_base1):
                            raise ValueError("Missing codon site" + str(codonsite_base1) + " in " + window_dnds_tsv)

                        outrow["N.Exp"] = codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_N_COL]
                        outrow["S.Exp"] = codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_S_COL]
                        outrow["EN.Exp"] = codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_EXP_N_COL]
                        outrow["ES.Exp"] = codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_EXP_S_COL]

                        if (codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_S_COL] and
                                    float(codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_S_COL]) != 0):
                            outrow["dNdS.Exp"] = (float(codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_DN_COL])/
                                                  float(codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_DS_COL]))

                        outrow["dN_minus_dS.Exp"] = codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL]

                        writer.writerow(outrow)


def make1csv_errfree(output_csv_filename, sim_args_tsv):
    """
    Puts all the collate_dnds, full population csv, expected dnds info into 1 csv for checking what causes inaccurate
    inferred dn/ds.
    Instead of putting read sequencing umberjack results in the csv, instead it puts umberjack results
    from the full population.  Ie the windows are made on the sequences from the full population instead of the ART reads.

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
                                                    "UnambigCodonRate.Act", # Total unambiguous codon (depth) at the codon site / max read depth for entire slice
                                                    "AADepth.Act",  # Total codons that code for only 1 amino acid at the codon site
                                                    "PopSize.Act",  # Population size
                                                    "ConserveCodon.Act",
                                                    "EntropyCodon.Act",  # Excludes codons with N's and gaps
                                                    "UnknownPerCodon.Act",  # Average N or gaps per codon at this site
                                                    "ErrPerCodon.Act",  # Average erroneous bases per codon at this site
                                                    "N.Act", "S.Act",
                                                    "EN.Act", "ES.Act",
                                                    "dNdS.Act",
                                                    "dN_minus_dS.Act",
                                                    "TreeLen.Act",  # length of window tree in nucleotide subs/site
                                                    "TreeDepth.Act",  # depth of longest branch in nucleotide subs/site
                                                    "Polytomy.Act",
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


        popn_groups, umberjack_group_to_args = run_sliding_window_tree.parse_sim_args_tsv(sim_args_tsv)
        for umberjackgroup, popn_groups_per_ugroup in umberjack_group_to_args.iteritems():
            for popn_group in popn_groups_per_ugroup:
                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/simdatasetname
                sim_popn_name = popn_group.dataset
                sim_data = SimData(popn_group.config_file)
                sim_data_dir = sim_data.sim_data_dir


                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/simdatasetname/subs/simdatasetname.dnds.tsv
                full_popn_dnds_tsv = sim_data_dir + os.sep + "subs" + os.sep + sim_popn_name + ".dnds.tsv"
                # /home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/simdatasetname/fullpopn/simdatasetname.conserve.csv
                full_popn_conserve_csv = sim_data_dir + os.sep + "fullpopn" + os.sep + sim_popn_name + "_TRUE.conserve.csv"


                # Instead of using umberjack window fasta made from reads, we use window made from error free ART reads with perfect alignment
                window_start = 1
                window_end = NUCSITES
                umberjack_output_prefix = DATASET_OUT_DIR + os.sep + sim_popn_name + os.sep + "errfree" + os.sep + sim_popn_name + ".repro.errfree.{}_{}".format(window_start, window_end)
                window_fasta = umberjack_output_prefix + ".fasta"
                window_treefile = umberjack_output_prefix + ".nwk"
                window_dnds_tsv = umberjack_output_prefix + ".dnds.tsv"



                LOGGER.debug("Merge sim_name=" + sim_popn_name + " full popn window dnds tsv =" + window_dnds_tsv)

                total_indiv = popn_group.indiv
                total_codon_sites = popn_group.codonsites

                #CodonSite	ConserveCodon	Entropy	NucDepth	CodonDepth
                codonsite_2_full_cons = collect_training.read_codon_csv(csv_file=full_popn_conserve_csv, codon_site_field="CodonSite", is_base0=False)

                # File,Window_Start,Window_End,Reads,CodonSite,CodonDepth,AADepth,ConserveAllCodon,EntropyAllCodon,ConserveCodon,EntropyCodon,N,S,EN,ES,dN,dS,dN_minus_dS,Ambig,Pad,Err,Err_N,Err_S,Ambig_N,Ambig_S,TreeLen,T

                # Site	Observed S Changes	Observed NS Changes	E[S Sites]	E[NS Sites]	dS	dN	dN-dS	Scaled dN-dS
                codonsite_2_full_dnds= collect_training.read_codon_csv(csv_file=full_popn_dnds_tsv, codon_site_field="Site", is_base0=True, delimiter="\t")

                if (len(codonsite_2_full_dnds.keys()) != len(codonsite_2_full_cons.keys()) or
                            len(codonsite_2_full_dnds.keys()) != total_codon_sites or
                            len(codonsite_2_full_cons.keys()) != total_codon_sites):
                    raise ValueError("full population dnds does not have same number of codon sites as conservation:",
                                     full_popn_dnds_tsv, ", ", full_popn_conserve_csv)

                aln = Utility.Consensus()
                aln.parse(msa_fasta_filename=window_fasta)
                window_reads = aln.get_total_seqs()
                window_tree_dist = TestTopology.calc_window_tree_dist(sim_data=sim_data,
                                       window_fasta=window_fasta,
                                       window_treefile=window_treefile,
                                       win_start=window_start, win_end=window_end)
                full_popn_breaks = sim_data.get_recombo_breaks()

                break_ratio = collect_training.get_break_ratio(sim_data=sim_data, win_start=window_start, win_end=window_end)

                polytomy_brlen_thresh = 1.0/(3 * total_codon_sites)  # branch length treshold below which node is considered polytomy
                window_treelen, window_treedepth, total_polytomies = collect_training.get_tree_len_depth(window_fasta, polytomy_brlen_thresh=polytomy_brlen_thresh)

                with open(window_dnds_tsv, 'rU') as fh_actual:
                    reader_act = csv.DictReader(fh_actual, delimiter="\t")
                    for row_idx, row_act in enumerate(reader_act):
                        act_codonsite_offset_base0 = int(row_act["Site"])
                        act_codonsite_base0 = act_codonsite_offset_base0 + window_start - 1


                        act_nucsite_offset_base0 = act_codonsite_base0 * 3
                        codonsite_base1 = act_codonsite_base0 + 1
                        unambig_codon_depth = aln.get_codon_depth(codon_pos_0based=act_codonsite_offset_base0, is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                        outrow = dict()
                        outrow["Window_Start"] = window_start
                        outrow["Window_End"] = window_end
                        outrow["CodonSite"] = codonsite_base1
                        outrow["File"] = "ErrFree_" + sim_data.name
                        outrow["Reads.Act"] = window_reads
                        outrow["UnambigCodonRate.Act"] = float(unambig_codon_depth)/window_reads
                        outrow["AADepth.Act"]  = aln.get_unambig_codon2aa_depth(codon_pos_0based=act_codonsite_offset_base0)
                        outrow["PopSize.Act"] = total_indiv
                        outrow["ConserveCodon.Act"] = aln.get_codon_conserve(codon_pos_0based=act_codonsite_offset_base0,
                                                                             is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                        outrow["EntropyCodon.Act"] = aln.get_codon_shannon_entropy(codon_pos_0based=act_codonsite_offset_base0,
                                                                             is_count_ambig=False, is_count_gaps=False, is_count_pad=False)
                        outrow["UnknownPerCodon.Act"] = float(aln.get_gap_count(pos_0based=act_nucsite_offset_base0) +
                                                              aln.get_ambig_count(pos_0based=act_nucsite_offset_base0) +
                                                              aln.get_pad_count(pos_0based=act_nucsite_offset_base0)) / window_reads
                        outrow["ErrPerCodon.Act"] = 0
                        # If it never made it past FastTree into hyphy, then the substitutions will be empty string
                        if row_act[hyphy.hyphy_handler.HYPHY_TSV_N_COL] and row_act[hyphy.hyphy_handler.HYPHY_TSV_S_COL]:
                            outrow["N.Act"] = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_N_COL])
                            outrow["S.Act"] = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_S_COL])
                            outrow["EN.Act"] = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_EXP_N_COL])
                            outrow["ES.Act"] = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_EXP_S_COL])
                            if row_act["dS"] and float(row_act[hyphy.hyphy_handler.HYPHY_TSV_S_COL]) != 0:
                                outrow["dNdS.Act"] = float(row_act[hyphy.hyphy_handler.HYPHY_TSV_DN_COL])/float(row_act[hyphy.hyphy_handler.HYPHY_TSV_DS_COL])
                            outrow["dN_minus_dS.Act"] = row_act[hyphy.hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL]
                        outrow["TreeLen.Act"] = window_treelen
                        outrow["TreeDepth.Act"] = window_treedepth
                        outrow["TreeDistPerRead.Act"] = float(window_tree_dist)/window_reads

                        outrow["Is_Break"] = 0
                        for nuc_strand_start_wrt_ref_base1, nuc_strand_end_wrt_ref_base1 in full_popn_breaks:
                            nuc_pos_wrt_ref_base1 = window_start + act_nucsite_offset_base0
                            # If there are no recombination breaks, full_popn_breaks still contains the full genome as a contiguous section
                            # Don't consider first position as breakpoint
                            if len(full_popn_breaks) > 1 and nuc_pos_wrt_ref_base1 == nuc_strand_start_wrt_ref_base1 > 1:
                                outrow["Is_Break"] = 1

                        outrow["BreakRatio.Act"] = break_ratio
                        outrow["Polytomy.Act"] = total_polytomies

                        if not codonsite_2_full_cons.get(codonsite_base1):
                            raise ValueError("Missing codon site" + str(codonsite_base1) + " in " + full_popn_conserve_csv)
                        outrow["ConserveCodon.Exp"] = codonsite_2_full_cons[codonsite_base1]["ConserveCodon"]
                        outrow["EntropyCodon.Exp"] = codonsite_2_full_cons[codonsite_base1]["EntropyCodon"]

                        if not codonsite_2_full_dnds.get(codonsite_base1):
                            raise ValueError("Missing codon site" + str(codonsite_base1) + " in " + window_dnds_tsv)

                        outrow["N.Exp"] = codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_N_COL]
                        outrow["S.Exp"] = codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_S_COL]
                        outrow["EN.Exp"] = codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_EXP_N_COL]
                        outrow["ES.Exp"] = codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_EXP_S_COL]

                        if (codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_S_COL] and
                                    float(codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_S_COL]) != 0):
                            outrow["dNdS.Exp"] = (float(codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_DN_COL])/
                                                  float(codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_DS_COL]))

                        outrow["dN_minus_dS.Exp"] = codonsite_2_full_dnds[codonsite_base1][hyphy.hyphy_handler.HYPHY_TSV_SCALED_DN_MINUS_DS_COL]

                        writer.writerow(outrow)

def concat_collate_csv(output_csv, append_csv):
    """
    Append training dnds csv  to existing training dnds csv
    :return:
    """
    with open(output_csv, 'a') as fh_out:
        # Strip out the header and append to file
        subprocess.check_call(('tail', '-n', '+2', append_csv), stdout=fh_out)

# TestCase:  make a window where all the sequences are exactly like the full population.  No recombination.
# (not reads, the full sequence)
#   - do umberjack on window
# Expected:  we expect that the window tree won't have the exact same topology due to greedy fasttree.
# But we expect the robinson foulds distance to be minimal.  This testcase should have the best concordance.
test_window_same_full(dataset="TestRecomboNone")

# TestCase:  make a window where all the sequences are exactly like the full population.  1 Recombination in middle.
#   - do umberjack on window
# Expected:  3rd best concordance.
test_window_same_full(dataset="TestRecomboMid")

# TestCase:  make a window where all the sequences are exactly like the full population.  Lots Recombination in middle.
#   - do umberjack on window
# Expected:  3rd best concordance.
test_window_same_full(dataset="TestRecomboLots")

# TestCase:  use the umberjack window from the dataset reads. But window is entire size of genome.  No recombination
# Expected:  best concordance out of windows made from dataset reads.
# No expecations on where concordance ranks wrt windows made from full population sequences.
test_window_reads(dataset="TestRecomboNone", window_start=1, window_end=NUCSITES)

# TestCase:  use the umberjack window from the dataset reads. But window is entire size of genome.  1 Recombination in middle.
# Expected:  2nd best concordance out of windows made from dataset reads.
# No expecations on where concordance ranks wrt windows made from full population sequences.
test_window_reads(dataset="TestRecomboMid", window_start=1, window_end=NUCSITES)

# TestCase:  use the umberjack window from the dataset reads. But window is entire size of genome.  Lots recombination in middle.
# Expected:  3rd best concordance out of windows made from dataset reads.
# No expecations on where concordance ranks wrt windows made from full population sequences.
test_window_reads(dataset="TestRecomboLots", window_start=1, window_end=NUCSITES)


# TestCase:  use the umberjack window from the error free dataset reads with the perfect read alignment (from ART).
# Window is entire size of genome.
# Expected:  concordance is better than umberjack run on imperfectly aligned reads (via BWA).
#
test_errFree_window_reads(dataset="TestRecomboNone", window_start=1, window_end=NUCSITES)
test_errFree_window_reads(dataset="TestRecomboMid", window_start=1, window_end=NUCSITES)
test_errFree_window_reads(dataset="TestRecomboLots", window_start=1, window_end=NUCSITES)

make1csv_same_full(output_csv_filename=OUTPUT_COLLATE_SAMEFULL_CSV, sim_args_tsv=SIM_ARGS_TSV)
make1csv_errfree(output_csv_filename=OUTPUT_COLLATE_ERRFREE_UMBERJACK_CSV, sim_args_tsv=SIM_ARGS_TSV)


subprocess.check_call(["python", "collect_training.py",
                       "-t", SIM_ARGS_TSV,
                      "-f", OUTPUT_COLLATE_UMBERJACK_CSV])


# concatenate the csv files
#shutil.copy(OUTPUT_COLLATE_UMBERJACK_CSV, OUTPUT_COLLATE_CSV)
#shutil.copy(OUTPUT_COLLATE_SAMEFULL_CSV, OUTPUT_COLLATE_CSV)
#concat_collate_csv(output_csv=OUTPUT_COLLATE_CSV, append_csv=OUTPUT_COLLATE_SAMEFULL_CSV)
#concat_collate_csv(output_csv=OUTPUT_COLLATE_CSV, append_csv=OUTPUT_COLLATE_ERRFREE_UMBERJACK_CSV)
