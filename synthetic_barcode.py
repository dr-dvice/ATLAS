# Author: David Bellini
# Date: 4/14/22
# A script to apply synthetic barcodes for pseudo single cell rna seq analysis
# Takes 2 arguments: 1. The keyfile path and 2. The directory with the fastq gzipped files. 3. The desired output dir.
# Script can be modified to accept more arguments: UMIsize, MAX_THREADS, etc.
import sys
from os import listdir
from os.path import isfile, join
import gzip
from multiprocessing.pool import ThreadPool

#keyFile = "res.csv"
#fastqgz_directory = "ff15/"
#filecodes = {}
#UMIsize = 8

keyFile = sys.argv[1]
fastqgz_directory = sys.argv[2]
outdir = sys.argv[3]
filecodes = {}
UMIsize = 8
MAX_THREADS = 8

# Retrieve all files in specified directory
fastq_files = [f for f in listdir(fastqgz_directory) if isfile(join(fastqgz_directory, f))]

print("Reading key file...", flush=True)
with open(keyFile) as cellkeys:
    for line in cellkeys:
        vals = line.split(",")
        filecodes[vals[0]] = vals[1]
print("Key file read", flush=True)


def Processfile(fq, outfilename):
    print("Processing " + fq + " ...", flush=True)
    fq_path = fastqgz_directory + fq
    blockbuilder = ""
    idLine = ""
    dnaLine = ""
    umi = ""
    barcode = filecodes[fq]
    # read content from original fastq file and prep for processing
    with gzip.open(fq_path, 'rt') as f:
        with gzip.open(outfilename, 'wt') as outfile:
            file_content = f.read()
            file_content = file_content.split("\n")
            nextLineDNA = False
            # loop through each read and save each piece of info into respective string tag
            for line in file_content:
                if line == '':
                    outfile.write(blockbuilder)
                    break
                elif nextLineDNA:
                    umi = line[0:UMIsize]
                    dnaLine = line[UMIsize:]
                    nextLineDNA = False
                elif line[0] == '@':
                    idLine = line
                    nextLineDNA = True
                elif line[0] == '+':
                    # do all the processing here
                    # here, we add the appropriate barcode and umi to the string
                    idBits = idLine.split(' ')
                    idBits[0] += "_" + barcode.rstrip() + "_" + umi
                    idBits = " ".join(idBits)
                    blockbuilder += (idBits + "\n")
                    # then we add in the stripped DNA
                    blockbuilder += (dnaLine + "\n")
                    blockbuilder += ('+' + '\n')

                    # write and reset
                    outfile.write(blockbuilder)
                    blockbuilder = idline = dnaLine = umi = ""

                else:  # This is the read quality line
                    # We need to strip the quality strings that matched the UMI. This information will be lost.
                    blockbuilder = (line[UMIsize:] + '\n')
    print(fq + " processing complete.", flush=True)

# initialize multithreading pool
pool = ThreadPool(processes=MAX_THREADS)

# loop through all the .fastq.gz files in the specified directory
for fq in (fastq_files):
    outfilename_gen = sys.argv[3] + "processed_" + fq
    # add the files to the multithreading handler
    pool.apply_async(Processfile, args=(fq, outfilename_gen))
pool.close()
pool.join()

print("All Files Processed!", flush=True)