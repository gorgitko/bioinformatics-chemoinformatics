import pandas as pd
import pickle

data_dir = "data/"
gene_id = "uc021pzf.2"
gene_file_suffix = ".gtf"
gene_file = "refGene-chr10-{}-exons".format(gene_id)
cigar_pickle_file = "transcript-{}-cigarIntervals.pickle".format(gene_id)
cyto_file = gene_file + ".cyto"

# load read intervals from serialized file
with open(data_dir + cigar_pickle_file, mode='br') as f:
    intervals = pickle.load(f)

# load exons from .GTF
exons = {}
data = pd.read_table(data_dir + gene_file + gene_file_suffix, sep="\t")
first = max(data.index)

# create dictionary with keys as exon names and values with starting and ending position of exon
# gene is on "-" strand so first exon is with highest starting position
# {"exon_XX": (start, end)}
for index, row in data.iterrows():
    exons["exon_{:02d}".format(first)] = (int(row["start"]), int(row["end"]))
    first -= 1

joined_exons = {}

# iterate through read intervals and add connected exons as keys and their counts to dictionary
# keys are exon tuples: e.g. (exon2, exon5)
for i in intervals:
    lastExon = -1
    for ii in i:
        for exon in sorted(exons.keys()):
            if (ii[0] <= exons[exon][0] < ii[1]) or (ii[0] < exons[exon][1] <= ii[1]) or \
                (exons[exon][0] <= ii[0] <= ii[1] <= exons[exon][1]) or (ii[0] <= exons[exon][0] <= exons[exon][1] <= ii[1]):
                if (lastExon, exon) not in joined_exons.keys():
                    if lastExon != exon and lastExon != -1:
                        joined_exons[(lastExon, exon)] = 1
                else:
                    joined_exons[(lastExon, exon)] += 1
                lastExon = exon

# write table for cytoscape
with open(cyto_file, mode="w") as f:
    f.write("NODE_1\tNODE_2\tINTERACTION\tSUM_LINKS\n")
    for exons in joined_exons.keys():
        f.write("{}\t{}\tpd\t{}\n".format(exons[0], exons[1], joined_exons[exons]))