# exports to pickled form the interval of exons in transcript reads for selected gene using CIGAR

import pandas as pd
import re
import pickle

def getReadIntervals(samfile):
    """
    Yields the intervals of reads from SAM file, computed from CIGAR.
    If read contains "N" or "D" CIGAR chars, yields all intervals: 10M5N50M10D => yield [(start1, end1), (start2, end2)]
    """
    re_cigar = re.compile("(\d+)([MIDNSHPX=])", re.I)
    reads = pd.read_table(samfile, sep="\t")

    for index, row in reads.iterrows():
        start = row["POS"]
        end = row["POS"]
        intervals = []
        cigar = re.findall(re_cigar, row["CIGAR"])
        for c in cigar:
            if c[1] == "M":
                end += int(c[0]) - 1
            elif c[1] == "N" or c[1] == "D":
                intervals.append((start, end))
                end += int(c[0]) - 1
                start = end
            elif c[1] == "I":
                end -= int(c[0]) - 1
        intervals.append((start, end))
        if len(intervals) > 1:
            yield intervals

gene_id = "uc021pzf.2"
gene_file = "data/refGene-chr10-{}-exons.gtf".format(gene_id)
transcript_file = "transcript"
transcript_file_suffix = ".sam"
transcript_filtered_file = "transcript-{}".format(gene_id)

# get first exon starting position
gene = pd.read_table(gene_file, sep="\t")
min_pos = int(gene.loc[[0]]["start"])

# filter reads - include only reads starting at position >= first exon starting position
reads = pd.read_table("data/" + transcript_file + transcript_file_suffix, sep="\t")
reads = reads[(reads["POS"] >= min_pos)]
reads.to_csv(transcript_filtered_file + transcript_file_suffix, index=False, sep="\t")

# intervals serialization ("pickling")
with open("data/{}-cigarIntervals.pickle".format(transcript_filtered_file), mode='bw') as f:
    pickle.dump(list(getReadIntervals(transcript_filtered_file + transcript_file_suffix)), f)