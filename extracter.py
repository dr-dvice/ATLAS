# Author: David Bellini
# Date: 4/8/22
# A simple script to extract a gene product description from a .gff file
import regex as re

filename = "genome_assemblies_genome_gff"
outfile = 'gffKey2.tsv'
lol = {}
with open(filename) as gffFile:
    for line in gffFile:
        if re.search(r"product=[^;]+", line):
            gene = re.search(r"gene=[^;=]+", line)
            gene = gene.group(0)
            if gene.endswith('product'):
                gene = gene[5:-7]
            else:
                gene = gene[5:]
            product = re.search(r"product=[^;]+", line)
            product = product.group(0)
            product = product.replace("product=", "")
            product = product.rstrip()
            if gene not in lol.keys():
                lol[gene] = product
with open(outfile, 'w') as output:
    for gen,pro in lol.items():
        output.writelines(gen + '\t' + pro + '\n')