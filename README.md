# Qiime2PhyloSeq
# Language: Python
# Input: TXT
# Output: PREFIX
# Tested with: PluMA 1.1, Python 3.6

PluMA plugin to take OTU and mapping TXT files from Qiime (Caporaso et al, 2010)
and convert into PhyloSeq (McMurdie and Holmes, 2013) format.

The plugin accepts as input a TXT file of keyword-value pairs, tab-delimited.

Keywords:
otufile: Qiime OTU TXT file
mapping: Qiime mapping TXT file

PhyloSeq otu_table, tax_table and metadata will be generated with the user-specified output prefix.

