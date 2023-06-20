# Author: David Bellini
# Date: 6.19.23
# A script to conduct all necessary preprocessing on raw single cell data for the cell atlas, total RNA format.

import sys
from os import listdir
from os.path import isfile, join
import gzip
import regex
from multiprocessing.pool import ThreadPool

#fastqgz_directory = sys.argv[1]
#outdir = sys.argv[2]
#sampleBCfileName = sys.argv[2]

MAX_THREADS = 8
fastqgz_directory = "/home/davidbellini/files/workspace/june_run/testbed/"
#fastqgz_directory = "C:/Users/david_n6jk0pe/OneDrive - Baylor College of Medicine/gabbiani/ATLAS/"
#sampleBCfileName = "samplebarcodes456.csv"
cellBCfileName = "/home/davidbellini/files/workspace/june_run/testbed/MATQ_Drop_384_barcodes.txt"
outputFolder = "/home/davidbellini/files/workspace/june_run/testbed/results/"

BClen = 10 #expected length of the barcodes
UMIlen = 8
MAX_HANDLE_POS = 15 #Handle is not expected to the farther than 16 bases from the start
#HANDLE_TO_READ_LEN = 38 #Expected distance from the end of the handle to the actual read.
handle = "TGGTAGGTGGTAGAGA"
handlePattern = "(" + handle + "){e<=3}"
#GAT21 = "GATGGTTGAGGATGTGTGGAG"
sampleBC_to_ID = {}
cellBC_codes = []

print("Reading Cell BC file...", flush=True)
with open(cellBCfileName) as cellBCfile:
    cellBC_codes = [line.rstrip('\n') for line in cellBCfile]
print("Cell BC file read", flush=True)

#fastq_files = [f for f in listdir(fastqgz_directory) if isfile(join(fastqgz_directory, f))]

def Processfile(fq):
    print("Processing " + fq + " ...", flush=True)
    fq_path = fastqgz_directory + fq
    blockbuilder = ""
    allHandles = {}
    falseHandles = {}
    noHandles = 0

    # debug stuff
    noHandleList = []
    faultyBarcodes = 0
    totalReads = 0
    read2short = 0
    open((outputFolder + "badBarcodes.txt"), "w").close() #clears badbarcodes file

    # read content from original fastq file and prep for processing
    with gzip.open(fq_path, 'rt') as f:
        name = fq.split(".")[0]
        outfile = gzip.open((outputFolder + name + "_processed.fq.gz"), 'wt')
        read2 = False
        rna_to_write = ""
        QC_line = ""
        idLine = ""
        nextLineDNA = False
        end_of_file = False
        debugMode = True
        # loop through each read and save each piece of info into respective string tag
        for line in f:
            if read2: #Read 2 is where we do most of our processing. Here we need to extract the barcodes and umi,
                # and verify the general integrity of the read.
                read2 = False
                # In the total RNA version, we find look for the handle and the barcodes in the second read.
                handleLoc = regex.search(handlePattern, line, regex.BESTMATCH)
                # check validity
                if handleLoc is None:
                    noHandles+=1
                    noHandleList.append(line)
                    SkipToNextRead(f)
                    ResetVars(idLine, rna_to_write, QC_line)
                    continue
                else:
                    startpos = handleLoc.regs[0][0]
                    endpos = handleLoc.regs[0][1]
                    if debugMode: allHandles[startpos] = allHandles.get(startpos, 0) + 1
                    if startpos > MAX_HANDLE_POS: #If the handle is way farther down than expected
                        falseHandles[startpos] = falseHandles.get(startpos, 0) + 1
                        SkipToNextRead(f)
                        ResetVars(idLine, rna_to_write, QC_line)
                        continue
                #check to see that the read is actually long enough to have all of the information we needed
                UMIStartpos = endpos + BClen
                if ((startpos - BClen) >= 0) and ((endpos + (BClen + UMIlen)) < len(line)):
                    barcode1 = line[(startpos - BClen):startpos]
                    barcode2 = line[endpos:(endpos + BClen)]
                    if debugMode: cBarcode = DebugBadBarcodes(barcode1, barcode2, line)
                    else: cBarcode = CorrectCellBarcodes(barcode1, barcode2)
                    if cBarcode is None:
                        faultyBarcodes+=1
                        SkipToNextRead(f)
                        ResetVars(idLine, rna_to_write, QC_line)
                        continue
                    else: # Checks successfull, here we actually take the read and the umi and create the new ID line for read 1.
                        newIDline = cBarcode + "_" + line[UMIStartpos:(UMIStartpos + UMIlen)]
                        idBits = idLine.split(' ')
                        idBits[0] += newIDline
                        idBits = " ".join(idBits)
                        # Time to jot everything down
                        blockbuilder = idBits + "\n" + rna_to_write + "\n" + "+" + QC_line
                        outfile.write(blockbuilder)
                else:
                    read2short+=1
                # Skip the next two lines in the file cause we don't care about this read aside from the BC info it gives
                SkipToNextRead(f)
                ResetVars(idLine, rna_to_write, QC_line)
                continue
            # Read 1 contains the RNA line. We will make no edits to this for now
            if nextLineDNA:
                rna_to_write = line
                nextLineDNA = False
                totalReads+=1
                if totalReads % 1000000 == 0:
                    print(name + ": " + str(totalReads) + " reads processed", flush=True)
                if end_of_file: break
                next(f) # I think we can just skip the + line entirely
            elif line[0] == '@':
                idLine = line
                nextLineDNA = True
            else:  # This is the read quality line
                QC_line = line
                read2 = True
                # we are going to skip the read 2 id line
                next(f)

        # close the file writer
        outfile.close()
    if debugMode:
        with open((outputFolder + "noHandleReads.txt"), 'w') as noHandle_file:
            for read in noHandleList:
                noHandle_file.write(read)
        with open((outputFolder + "handle_stats.txt"), 'w') as handlestats_file:
            for key in sorted(allHandles.keys()):
                handlestats_file.write(str(key) + ": " + str(allHandles[key]) + '\n')
        # Write the handles

    print(fq + " processing complete:", flush=True)
    falseHandlesNum = sum(falseHandles.values())
    print("\t-" + str(falseHandlesNum) + " False Handles detected.")
    print("\t-" + str(noHandles) + " reads had no valid handle.")
    print("\t-" + str(read2short) + " reads were too short to read the barcodes")
    print("\t-" + str(faultyBarcodes) + " reads with at least one faulty barcodes")
    print("\t-" + str(totalReads) + " total reads")
    # Reads are discarded in this order and are cumulative, so they are output in this order.

    badHandle_percent = round(((noHandles + falseHandlesNum) / totalReads) * 100, 1)
    read2short_percent = round((read2short / totalReads) * 100, 1) #rounds to nearest tenth
    faultyBarcodes_percent = round((faultyBarcodes / totalReads) * 100, 1)
    print("\t-" + "Bad Handle Percentage: " + str(badHandle_percent) + "%")
    print("\t-" + "Short read percentage: " + str(read2short_percent) + "%")
    print("\t-" + "Bad Cell BC Percentage: " + str(faultyBarcodes_percent) + "%")
    print("\t-" + "Total read percentage discarded: " + str(round((read2short_percent + badHandle_percent + faultyBarcodes_percent), 1 )) + "%")

# A function that writes bad barcodes and reads to a file for debug
# This function significantly increases runtime, when bc are not being debugged please opt to use normal "CorrectCellBarcodes" instead.
def DebugBadBarcodes(bc1, bc2, read):
    with open((outputFolder + "badBarcodes.txt"), 'a') as badBarcode_file:
        BC1_check = CheckCellBarcode(bc1)
        BC2_check = CheckCellBarcode(bc2)
        if BC1_check is None:  # If either barcode is faulty
            badBarcode_file.write("BC1 is bad: " + bc1 + " ")
        if BC2_check is None:
            badBarcode_file.write("BC2 is bad: " + bc2 + " ")
        if None in [BC1_check, BC2_check]:
            badBarcode_file.write("BAD read: " + read + '\n')
            return None
        else: return BC1_check + BC2_check

def CorrectCellBarcodes(bc1, bc2):
    BC1_check = CheckCellBarcode(bc1)
    BC2_check = CheckCellBarcode(bc2)
    if None in [BC1_check, BC2_check]:
        return None
    else: return BC1_check + BC2_check

def CheckCellBarcode(bc):
    if bc in cellBC_codes:
        return bc
    else:
        pattern = "(" + bc + "){e<=2}"
        for cellBC in cellBC_codes:
            match = regex.search(pattern, cellBC, regex.BESTMATCH)
            if match is not None: return cellBC
        return None

def SkipToNextRead(file):
    try:
        for i in range(2): #changed to 2 from 6 of dt
            next(file)
    except StopIteration as e:
        print(e)

def ResetVars(id, rna, QC):
    id = rna = QC = ""


# pool = ThreadPool(processes=MAX_THREADS)
#
# # for fq in fastq_files:
# #     #add the files to the multithreading handler
# #     result = pool.apply_async(Processfile, args=(fq, 1))
# #     #print(result.successful())
#
# pool.close()
# pool.join()

#print("All files processed!", flush=True)

Processfile("millionreads.fq.gz")
#Processfile("ye.fq.gz")