# Find the error rate of ART simulated reads by parsing the cigar strings in the SAM file produced by Art


#tu3527-read2   147     otu3527 4646    99      17=1X93=1X18=1X119=     =
import csv
import re
import subprocess

CIGAR_RE = re.compile('[0-9]+[MIDNSHPX=]')

class AlignType:
    ALN_HIT = "M"
    INSERTION = "I"
    DELETION = "D"
    SKIPPED = "N"
    SOFTCLIP = "S"
    HARDCLIP = "H"
    PAD = "P"
    SEQ_MATCH = "="
    SEQ_MISMATCH = "X"

def apply_cigar (cigar, seq, qual):
    """
    Parse SAM CIGAR and apply to the SAM nucleotide sequence.
    Remove soft-clipped sequences.  They may be valid polymorphisms but there is no alignment information for clipped
    sequences so we have no way to know how they line up with other sequences.
    Bases with low quality are not removed here - that can be done with merge_pairs()
    Left-pads and right-pads sequence so that it lines up with the reference.

    :return:  tuple [left and right padded sequence, left and right padded quality]
    :rtype : tuple [str, str]
    :param str cigar: SAM cigar field
    :param str seq: SAM sequence field
    :param str qual: SAM quality field
    """
    newseq = ''
    newqual = ''
    tokens = CIGAR_RE.findall(cigar)
    if len(tokens) == 0:
        return None, None

    left = 0

    for token in tokens:
        length = int(token[:-1])

        if token[-1] == 'S':
            left += length

        # Matching sequence: carry it over
        elif token[-1] == 'M':
            newseq += seq[left:(left+length)]
            newqual += qual[left:(left+length)]
            left += length

        # Deletion relative to reference: pad with gaps
        elif token[-1] == 'D':
            newseq += '-'*length
            newqual += '!'*length 		# Assign fake placeholder score (Q=-1)

        # Insertion relative to reference:
        elif token[-1] == 'I':
            # newseq += seq[left:(left+length)]
            # newqual += qual[left:(left+length)]
            left += length
            continue

        else:
            raise Exception("Unable to handle CIGAR token: {} - quitting".format(token))

    return newseq, newqual


def get_unclipped_seq_start (cigar, pos):
    """
    Parse SAM CIGAR.
    Returns the start position of the unclipped sequence.

    :param str cigar: SAM cigar field
    :param int pos: SAM pos field indicating start position of clipped sequence.
    """
    tokens = CIGAR_RE.findall(cigar)
    if len(tokens) == 0:
        return None

    clip_len  = 0
    first_token = tokens[0]
    length = int(first_token[:-1])
    if first_token[-1] == 'S':
        clip_len  = length

    return pos - clip_len


def get_unclipped_start_end (cigar, pos):
    """
    Parse SAM CIGAR.
    Returns the (start, end) position of the unclipped sequence.  1based positiions

    :param str cigar: SAM cigar field
    :param int pos: SAM pos field indicating start position of clipped sequence.
    """
    tokens = CIGAR_RE.findall(cigar)
    if len(tokens) == 0:
        return None

    clip_len  = 0
    first_token = tokens[0]
    length = int(first_token[:-1])
    if first_token[-1] == 'S' or first_token[-1] == "H":
        clip_len  = length

    return pos - clip_len


def count_mismatch(samfile):
    with open(samfile, 'rU') as fh:
        for line in fh:  # skip the headers
            if line.startswith("@PG"):
                break

        total_mismatch = 0
        total_insert = 0
        total_deletion = 0
        total_hardclip = 0
        total_softclip = 0
        total_reads = 0
        total_bases_aligned = 0
        total_bases = 0
        total_match = 0
        reader = csv.DictReader(fh, delimiter="\t", fieldnames=["Qname",
                                            "Flag",
                                            "Rname",
                                            "Pos",
                                            "Mapq",
                                            "Cigar",
                                            "Rnext",
                                            "Pnext",
                                            "Tlen",
                                            "Seq",
                                            "Qual"])


        for row in reader:
            total_reads += 1
            cigar = row["Cigar"]
            tokens = CIGAR_RE.findall(cigar)
            for token in tokens:
                # last letter in token indicates the type of alignment.  The number prior indicate the length of the alignment
                align_type = token[-1]
                align_len = int(token[:-1])

                if align_type != AlignType.DELETION and align_type != AlignType.PAD:
                    total_bases += align_len

                if align_type  != AlignType.HARDCLIP and align_type != AlignType.SOFTCLIP and align_type != AlignType.DELETION and align_type != AlignType.PAD:
                    total_bases_aligned += align_len

                if align_type == AlignType.SEQ_MISMATCH:
                    total_mismatch += align_len
                if align_type == AlignType.SEQ_MATCH or align_type == AlignType.ALN_HIT:
                    total_match += align_len
                elif align_type == AlignType.DELETION:
                    total_deletion += align_len
                elif align_type == AlignType.INSERTION:
                    total_insert += align_len
                elif align_type == AlignType.HARDCLIP:
                    total_hardclip += align_len
                elif align_type == AlignType.SOFTCLIP:
                    total_softclip += align_len

        print "\nResults for " + samfile
        print "========================="
        print "Total Reads=" + str(total_reads)
        print "Total Mismatch=" + str(total_mismatch) + " Fraction of Bases=" + str(float(total_mismatch)/total_bases)
        print "Total Insert=" + str(total_insert) + " Fraction of Bases=" + str(float(total_insert)/total_bases)
        print "Total Deletion=" + str(total_deletion) + " Fraction of Bases=" + str(float(total_deletion)/total_bases)
        print "Total Hard Clip=" + str(total_hardclip) + " Fraction of Bases=" + str(float(total_hardclip)/total_bases)
        print "Total Soft Clip=" + str(total_softclip) + " Fraction of Bases=" + str(float(total_softclip)/total_bases)
        print "Total Seq Bases Aligned=" + str(total_bases_aligned) + " Fraction of Bases=" + str(float(total_bases_aligned)/total_bases)
        print "Total Bases Hit (M= cigar, no indel)=" + str(total_match) + " Fraction of Bases=" + str(float(total_match)/total_bases)
        print "Total Bases=" + str(total_bases)



def get_next_record(sam_handle):
    for line in sam_handle:
        if not line.startswith("@"):
            qname, flag_str, rname, pos_str, mapq_str, cigar  = line.split("\t")[:6]
            flag = int(flag_str)
            pos = int(pos_str)
            mapq = int(mapq_str)
            return qname, flag, rname, pos, mapq, cigar

    return None

def is_record_same_ahead(query1, flag1, query2, flag2):
    """
    Returns if record1 is at same or ahead of record 2in terms of query sorted.
    :param query1:
    :param flag1:
    :param query2:
    :param flag2:
    :return:
    """
    IS_FIRST =                 0x040
    IS_SECOND =                0x080
    IS_SECONDARY_ALIGNMENT =   0x100
    IS_CHIMERIC_ALIGNMENT =    0x800
    if ((query1 > query2) or
            (query1 == query2 and flag1 & IS_FIRST and flag2 & IS_FIRST) or
            (query1 == query2 and flag1 & IS_SECOND and flag2 & IS_SECOND) or
            (query1 == query2 and flag1 & IS_SECOND and flag2 & IS_FIRST)):
        return True
    return False

def is_record_ahead(query1, flag1, query2, flag2):
    """
    Returns if record1 is at same or ahead of record 2in terms of query sorted.
    :param query1:
    :param flag1:
    :param query2:
    :param flag2:
    :return:
    """
    IS_FIRST =                 0x040
    IS_SECOND =                0x080
    IS_SECONDARY_ALIGNMENT =   0x100
    IS_CHIMERIC_ALIGNMENT =    0x800
    if ((query1 > query2) or
            (query1 == query2 and flag1 & IS_FIRST and flag2 & IS_SECOND)):
        return True
    return False

def is_record_same(query1, flag1, query2, flag2):
    """
    Returns if record1 is at same or ahead of record 2in terms of query sorted.
    :param query1:
    :param flag1:
    :param query2:
    :param flag2:
    :return:
    """
    IS_FIRST =                 0x040
    IS_SECOND =                0x080
    IS_SECONDARY_ALIGNMENT =   0x100
    IS_CHIMERIC_ALIGNMENT =    0x800

    if ((query1 == query2 and flag1 & IS_FIRST and flag2 & IS_FIRST) or
            (query1 == query2 and flag1 & IS_SECOND and flag2 & IS_SECOND)):
        return True
    return False



def cmp_art_aln_sam(art_samfile, aln_samfile):
    """
    Compares the true read placement as generated by the ART read simulator vs aligned read placement as determined by an aligner.
    ASSUMES the art sam file and the aln_samefile have been queryname sorted.
    ASSUMES that the art sam file will contain every read.
    ASSUMES that the aln sam file can have >=0 alignments for every read.
    :param art_samfile:
    :param aln_samfile:
    :return:
    """
    total_wrong = 0
    total_art_mate_alns = 0.0
    print "\nComparing Alignments for " + art_samfile + " and " + aln_samfile
    print "========================="
    with open(art_samfile, 'rU') as fh_art, open(aln_samfile, 'rU') as fh_aln:

        is_more_art = True
        last_qname_art = None
        last_qname_aln = None
        last_flag_art = None
        last_flag_aln = None
        last_pos_art = None

        qname_art = None
        flag_art = 0
        rname_art = None
        mapq_art = None
        pos_art = None
        cigar_art = None
        art_rec = qname_art, flag_art, rname_art, pos_art, mapq_art, cigar_art

        qname_aln = None
        flag_aln = 0
        rname_aln = None
        mapq_aln = None
        pos_aln = None
        cigar_aln = None
        aln_rec = qname_aln, flag_aln, rname_aln, pos_aln, mapq_aln, cigar_aln

        while is_more_art:
            if not last_qname_art or is_record_same_ahead(last_qname_aln, last_flag_aln, last_qname_art, last_flag_art):  # advance art if the aln is further ahead or same
                art_rec = get_next_record(fh_art)
                if art_rec:
                    total_art_mate_alns += 1
                    qname_art, flag_art, rname_art, pos_art, mapq_art, cigar_art  = art_rec
                else:
                    is_more_art = False



            if not last_qname_aln or is_record_same_ahead(last_qname_art, last_flag_art, last_qname_aln, last_flag_aln):  # advance aln if art is further ahead or same
                aln_rec = get_next_record(fh_aln)
                if aln_rec:
                    qname_aln, flag_aln, rname_aln, pos_aln, mapq_aln, cigar_aln  = aln_rec

            if is_record_ahead(qname_aln, flag_aln, qname_art, flag_art):
                print("Missing alignment for " + qname_art + " flag=" + str(flag_art) + " current alnqname=" + qname_aln)
                total_wrong += 1
            elif is_record_same(qname_aln, flag_aln, qname_art, flag_art):
                if pos_aln != pos_art:
                    unclipped_pos_aln = get_unclipped_seq_start(cigar=cigar_aln, pos=pos_aln)
                    if unclipped_pos_aln != pos_art:
                        print "Wrong position for " + qname_art + " artpos=" + str(pos_art) + " alnpos(unclipped)=" + str(unclipped_pos_aln)
                        total_wrong += 1
            elif is_record_same(qname_aln, flag_aln, last_qname_art, last_flag_art):
                if pos_aln != last_pos_art:
                    unclipped_pos_aln = get_unclipped_seq_start(cigar=cigar_aln, pos=pos_aln)
                    if unclipped_pos_aln != last_pos_art:
                        print "Wrong position for " + last_qname_art + " artpos=" + str(last_pos_art) + " alnpos(unclipped)=" + str(unclipped_pos_aln)
                        total_wrong += 1

            last_qname_art = qname_art
            last_qname_aln = qname_aln
            last_flag_art = flag_art
            last_flag_aln = flag_aln
            last_pos_art = pos_art





    print "\nComparison results for " + art_samfile + " and " + aln_samfile
    print ("Total Mates wrong = " + str(total_wrong) + "\n Fraction Mates wrong = " + str(total_wrong/total_art_mate_alns))










if __name__ == "__main__":
    # ART_SAMFILE = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/indelible/cov1x/sample_genomes.100.1x.sam"
    # count_mismatch(ART_SAMFILE)
    #
    # BOWTIE_CONSENSUS_ART_SAMFILE = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/indelible/cov1x/sample_genomes.100.1x.consensus.sam"
    # count_mismatch(BOWTIE_CONSENSUS_ART_SAMFILE)
    #
    # CFE_SAMFILE = "/home/thuy/projects/miseq_metrics/simulated_reads/sample_genomes.100.1x.sam"
    # count_mismatch(CFE_SAMFILE)

    # BOWTIE_CONSENSUS_CFE_SAMFILE = "/home/thuy/projects/miseq_metrics/simulated_reads/alignedVsConsensus/sample_genomes.100.1x.consensus.sam"
    # count_mismatch(BOWTIE_CONSENSUS_CFE_SAMFILE)
    #
    # ART_SMALL_SAMFILE = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small/mixed/reads/small.mixed.reads.sam"
    # count_mismatch(ART_SMALL_SAMFILE)
    #
    # BOWTIE_CONSENSUS_ART_SMALL_SAMFILE = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small/mixed/aln/small.mixed.reads.consensus.bowtie.sort.sam"
    # count_mismatch(BOWTIE_CONSENSUS_ART_SMALL_SAMFILE)
    #
    # ART_SMALL_COV5_SAMFILE = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov5.indiv200.codon400/mixed/reads/small.cov5.indiv200.codon400.mixed.reads.sam"
    # count_mismatch(ART_SMALL_COV5_SAMFILE)
    #
    # ART_SMALL_COV2_SAMFILE = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov2.indiv200.codon400/mixed/reads/small.cov2.indiv200.codon400.mixed.reads.sam"
    # count_mismatch(ART_SMALL_COV2_SAMFILE)
    #
    # ART_SMALL_COV10_SAMFILE = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov10.indiv200.codon400/mixed/reads/small.cov10.indiv200.codon400.mixed.reads.sam"
    # count_mismatch(ART_SMALL_COV10_SAMFILE)


    # ART_SMALL_COV2_INDIV1k_SAMFILE = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov2.indiv1k.codon400.bwa.rand/mixed/reads/small.cov2.indiv1k.codon400.bwa.rand.mixed.reads.sort.query.sam"
    # count_mismatch(ART_SMALL_COV2_INDIV1k_SAMFILE)
    #
    #
    # BWA_SMALL_COV2_INDIV1k_SAMFILE = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov2.indiv1k.codon400.bwa.rand/mixed/aln/small.cov2.indiv1k.codon400.bwa.rand.mixed.reads.consensus.bwa.sort.query.sam"
    # count_mismatch(BWA_SMALL_COV2_INDIV1k_SAMFILE)
    # # Quantify the alignment errors introduced by BWA for typical reads
    # cmp_art_aln_sam(ART_SMALL_COV2_INDIV1k_SAMFILE, BWA_SMALL_COV2_INDIV1k_SAMFILE)
    #
    #
    # ART_SMALL_COV2_INDIV1k_ERRFREE_SAMFILE = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov2.indiv1k.codon400.bwa.rand/mixed/reads/small.cov2.indiv1k.codon400.bwa.rand.mixed.reads.errFree.sort.query.sam"
    # count_mismatch(ART_SMALL_COV2_INDIV1k_ERRFREE_SAMFILE)
    #
    # BWA_SMALL_COV2_INDIV1k_ERRFREE_SAMFILE = "/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov2.indiv1k.codon400.bwa.rand/mixed/aln/small.cov2.indiv1k.codon400.bwa.rand.mixed.reads.errFree.consensus.bwa.sort.query.sam"
    # count_mismatch(BWA_SMALL_COV2_INDIV1k_ERRFREE_SAMFILE)
    #
    # # Quantify the alignment errors introduced by BWA for error free reads
    # cmp_art_aln_sam(ART_SMALL_COV2_INDIV1k_ERRFREE_SAMFILE, BWA_SMALL_COV2_INDIV1k_ERRFREE_SAMFILE)
    # BWA_HCV_SEQ = "/home/thuy/gitrepo/HCV_genotype/out2/isolate_gtype_NS5B/genotype/6/HCV_AnyRegion_Genotype6.pairaln.bwa.sam"
    # count_mismatch(BWA_HCV_SEQ)


    # For Art's fake miseq reads
    HXB2_ENV_FAKE_MISEQ = "/home/thuy/projects/fake_miseq/hxb2_env.sam"
    count_mismatch(HXB2_ENV_FAKE_MISEQ)


