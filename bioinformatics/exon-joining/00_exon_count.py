# prints the number of exons for each gene in reference genome
# I made this to choose the gene with optimal count of exons

import pandas as pd

min_pos = 110000000
max_pos = 115000000

genes = pd.read_table("data/refGene-chr10.gtf", sep="\t")
genes.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
genes = genes[(genes["start"] >= min_pos) & (genes["end"] <= max_pos) & (genes["feature"] == "exon")]
grouped_count = pd.DataFrame(genes.groupby("attribute").size())
grouped_count.columns = ["count"]

print(grouped_count.sort_values("count"))