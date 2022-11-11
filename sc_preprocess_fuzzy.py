# Author: David Bellini
# Date: 11/10/22
# A script to conduct all necessary preprocessing on raw single cell data for the cell atlas.
import sys
from os import listdir
from os.path import isfile, join
import gzip
#import re
import regex
from multiprocessing.pool import ThreadPool

#fastqgz_directory = sys.argv[1]
#outdir = sys.argv[2]
#sampleBCfile = sys.argv[3]
fastqgz_directory = "/home/davidbellini/files/atlas/data/"
sampleBCfileName = "samplebarcodes.csv"
cellBCfileName = "cellbarcodes.txt"


BClen = 8 #expected length of the barcodes
MAX_HANDLE_POS = 17 #Handle is not expected to the farther than 16 bases from the start
HANDLE_TO_READ_LEN = 38 #Expected distance from the end of the handle to the actual read.
handle = "GAGAGCGAGACAGAAG"
GAT21 = "GATGGTTGAGGATGTGTGGAG"
sampleBC_codes = {}
cellBC_codes = []
cellBC_codes2 = "--".join(cellBC_codes)

print("Reading sample barcode file...", flush=True)
with open(sampleBCfileName) as sampleBCfile:
    for line in sampleBCfile:
        vals = line.split(",")
        sampleBC_codes[vals[0]] = vals[1]
print("Sample file read", flush=True)
print("Reading Cell BC file...", flush=True)
with open(cellBCfileName) as cellBCfile:
    cellBC_codes = [line.rstrip('\n') for line in cellBCfile]
print("Cell BC file read", flush=True)

#fastq_files = [f for f in listdir(fastqgz_directory) if isfile(join(fastqgz_directory, f))]

def Processfile(fq, outfilename):
    print("Processing " + fq + " ...", flush=True)
    fq_path = fastqgz_directory + fq
    blockbuilder = ""
    #sampleBC = sampleBC_codes[fq]
    falseHandles = 0
    noHandles = 0

    # debug stuff
    noHandleList = []

    # read content from original fastq file and prep for processing
    with gzip.open(fq_path, 'rt') as f:
        with gzip.open(outfilename, 'wt') as outfile:
            read2 = False
            curBC = ""
            umi = ""
            rna_to_write = ""
            trimmed_QC = ""
            idLine = ""
            nextLineDNA = False
            end_of_file = False
            # loop through each read and save each piece of info into respective string tag
            for line in f:
                if read2:
                    sampleBC = line[:8]
                    read2 = False
                    # do all the processing here
                    # we add the appropriate barcode and umi to the string
                    idBits = idLine.split(' ')
                    idBits[0] += "_" + curBC.rstrip() + "_" + sampleBC + "_" + umi
                    idBits = " ".join(idBits)
                    blockbuilder += (idBits)
                    # then we add in the stripped DNA
                    blockbuilder += (rna_to_write)
                    # the plus
                    blockbuilder += ('+' + '\n')
                    # and the stripped QC scores
                    blockbuilder += (trimmed_QC)
                    # write and reset
                    outfile.write(blockbuilder)
                    blockbuilder = curBC = rna_to_write = umi = ""
                    # Skip the next two lines in the file cause we don't care about this read aside from the UMI
                    next(f)
                    next(f)
                    continue
                #We have reached the end of the file, write and break.
                if line == '':
                    outfile.write(blockbuilder)
                    break
                if nextLineDNA:
                    nextLineDNA = False
                    curBC = rna_to_write = trimmed_QC = ""
                    # TODO: Fuzzy handling
                    pattern = "(" + handle + ")"
                    pattern_fuzzy = "(" + handle + "){e<=1}"
                    # without bestmatch, it will match other nucleotides even if the handle is perfect
                    match = regex.search(pattern_fuzzy, line, regex.BESTMATCH)
                    if match == None:
                        noHandles+=1
                        noHandleList.append(line)
                        #skip to start of next read 1 and ignore rest
                        try:
                            for i in range(6):
                                next(f)
                        except StopIteration as e:
                            print(e)
                            end_of_file = True
                        idLine = ""
                        continue
                    handle_start = match.regs[0][0]
                    handle_end = match.regs[0][1]
                    if handle_start > MAX_HANDLE_POS:
                        falseHandles+=1
                        try:
                            for i in range(6):
                                next(f)
                        except StopIteration as e:
                            print(e)
                            end_of_file = True
                        idLine = ""
                        continue
                    BC1 = line[(handle_start-BClen):handle_start]
                    BC2 = line[handle_end:(handle_end + BClen)]
                    curBC = BC1 + "_" + BC2
                    umi = line[:(handle_start-BClen)]
                    rna_to_write = line[(handle_end+HANDLE_TO_READ_LEN):]
                    if end_of_file: break
                    next(f) # I think we can just skip the + line entirely
                elif line[0] == '@':
                    idLine = line
                    nextLineDNA = True
                else:  # This is the read quality line
                    #Skip next ID line and signal to loop to strip the barcode from the read 2, that is
                    #coming up afterwards
                    trimmed_QC = line[(handle_end+HANDLE_TO_READ_LEN):]
                    next(f)
                    read2 = True
    with open("noHandleReads.txt", 'w') as noHandle_file:
        for read in noHandleList:
            noHandle_file.write(read)
    print(fq + " processing complete.", flush=True)
    print(str(falseHandles) + " False Handles detected.")
    print(str(noHandles) + " reads had no valid handle.")

def CheckCellBarcode(bc):
    match = ""
    if bc in cellBC_codes:
        return bc
    else:
        pattern = "(" + bc + "){e<=1}"
        match = regex.match(pattern, cellBC_codes2)
        if match is None: return None
        else: return match.string





Processfile("OL456_D1_merged.fq.gz", "OL456_D1_processed.fq.gz")
#Processfile("test_merged.fq.gz", "test_merged_processed.fq.gz")