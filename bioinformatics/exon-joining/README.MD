# Exon joining

Requires the [`pandas`](http://pandas.pydata.org/) module.

We have some reference genome (I thinks its mouse's) and transcript reads from
RNA-seq, which is already processed. Our task is to find the exon joining sites
for one chosen gene and visualize them (I am using the [Cytoscape](http://www.cytoscape.org/).
By counting the transcription variants we can get the relative amounts of gene splice variants.
**To simplify this problem we assume that any overlap between read and exon site is ok.**

So for example we have a gene (taken from reference genome) containing exons (copy this to some text editor if it is wrapped):

A  T  C  G  G  T  C  A  A  C  C  G  T  A  C  T  G  G  T  C  G  G  C  T  G  C  C  G  T  G  T  G  C  G  T  A  C
<---intron---><---exon1---><--intron---><---exon2---><--------coding---------><---------exon3--------------->
1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 25 26 27 28 29 30 31 32 33 34 35 36

Some reads map to position 2-34: the aren't any exons = any spliced sites

Some reads map to position 4-12 and 20-34: there was probably spliced **exon2** in this trancript so
the transcription variant is exon1-exon3.

Some reads map to position 16-35: there was probably spliced **exon1** in this transcript so
the transcription variant is exon2-exon3

**reference genome file:** *refGene-chr10.gtf*

**transcript file:** *transcript.sam*

First we get some gene from reference genome (we can use the script *00_exon_count.py*
which tells us how many exons does each gene have). I have chosen the gene called
*uc021pzf.2*.

Then we filter the related information of this gene from reference genome with *01_filter-reference.py*.

Now we can get interval of exons for our gene from transcript reads with *02_cigar_intervals.py*.

Finally we can compute the exon joins with *03_exons.py*. The resulting output for
Cytoscape is:

NODE_1	NODE_2	INTERACTION	SUM_LINKS

exon_06	exon_05	pd	41

exon_12	exon_11	pd	157

exon_04	exon_03	pd	40

exon_19	exon_18	pd	215

...

and visualized (just a part):

![](https://ctrlv.cz/shots/2016/03/31/3s36.png)
