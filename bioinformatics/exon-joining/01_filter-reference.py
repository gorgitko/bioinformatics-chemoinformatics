# exports the gene (resp. the rows for one gene) from reference genome

import pandas as pd

gene_id = "uc021pzf.2"

reference_file = "refGene-chr10"
reference_file_suffix = ".gtf"
output_file = "{}-{}-exons{}".format(reference_file, gene_id, reference_file_suffix)

# only positions in this interval
min_pos = 110000000
max_pos = 115000000

# filter genes (starting/ending position, only exons, only one gene) + add GTF header
genes = pd.read_table(reference_file + reference_file_suffix, sep="\t", header=None)
genes.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
genes = genes[(genes["start"] >= min_pos) & (genes["end"] <= max_pos) & (genes["feature"] == "exon") & (genes["attribute"].str.contains(gene_id))]
genes.to_csv(output_file, index=False, sep="\t")