# Author: David Bellini
# Date: 11/10/22
# A script to conduct all necessary preprocessing on raw single cell data for the cell atlas.
# PLEASE HAVE FILE NAMES IN THE FORMAT OF: OL_ID1_ID2_ID3_LANENUM_merged.fq.gz
# For example OL_04_05_06_D1_CKD220028955-1A_HK7KLDSX5_L3_1.fq.gz

import sys
from os import listdir
from os.path import isfile, join
import gzip
import regex
from multiprocessing.pool import ThreadPool

fastqgz_directory = sys.argv[1]
#outdir = sys.argv[2]
sampleBCfileName = sys.argv[2]

MAX_THREADS = 8
#fastqgz_directory = "/home/davidbellini/files/atlas/data/"
#fastqgz_directory = "C:/Users/david_n6jk0pe/OneDrive - Baylor College of Medicine/gabbiani/ATLAS/"
#sampleBCfileName = "samplebarcodes456.csv"
cellBCfileName = "cellbarcodes.txt"
outputFolder = "/home/davidbellini/files/atlas/data/"

BClen = 8 #expected length of the barcodes
MAX_HANDLE_POS = 17 #Handle is not expected to the farther than 16 bases from the start
HANDLE_TO_READ_LEN = 38 #Expected distance from the end of the handle to the actual read.
handle = "GAGAGCGAGACAGAAG"
GAT21 = "GATGGTTGAGGATGTGTGGAG"
sampleBC_to_ID = {}
cellBC_codes = []

print("Reading sample barcode file...", flush=True)
with open(sampleBCfileName) as sampleBCfile:
    for line in sampleBCfile:
        vals = line.split(",")
        #sampleBC_codes[vals[0]] = vals[1]
        sampleBC_to_ID[vals[1].rstrip()] = vals[0]  #Barcode -> File name
print("Sample file read", flush=True)
print("Reading Cell BC file...", flush=True)
with open(cellBCfileName) as cellBCfile:
    cellBC_codes = [line.rstrip('\n') for line in cellBCfile]
print("Cell BC file read", flush=True)

#by joining all the codes into one searchable string we avoid having very long for loops that impact performance
#cellBC_codes2 = " ".join(cellBC_codes)
#had to discontinue this as it was causing matches of 7 nucleotides when I need the whole 8 nucl barcode

fastq_files = [f for f in listdir(fastqgz_directory) if isfile(join(fastqgz_directory, f))]

def Processfile(fq, number):
    print("Processing " + fq + " ...", flush=True)
    fq_path = fastqgz_directory + fq
    blockbuilder = ""
    falseHandles = 0
    noHandles = 0

    # debug stuff
    noHandleList = []
    faultyBarcodes = 0
    totalReads = 0
    faultySampleBC = 0
    savedRead1 = ""
    open("badBarcodes.txt", "w").close() #clears badbarcodes file

    # Create a file writer for each sample present
    filename_list = fq.split("_")
    sample_filewriters = GetFileWriters(filename_list)

    # read content from original fastq file and prep for processing
    with gzip.open(fq_path, 'rt') as f:
        read2 = False
        curBC = ""
        #umi = ""
        rna_to_write = ""
        trimmed_QC = ""
        idLine = ""
        nextLineDNA = False
        end_of_file = False
        # loop through each read and save each piece of info into respective string tag
        for line in f:
            if read2:
                read2 = False
                sampleBC = line[:8]
                sampleBC_check = CheckSampleBarcode(sampleBC)
                if sampleBC_check is None:
                    # Log and skip read
                    faultySampleBC += 1
                    next(f)
                    next(f)
                    continue
                # get correct file writer
                outfile = sample_filewriters[sampleBC_to_ID[sampleBC_check]]
                # do all the processing here. First we add the appropriate barcode and umi to the string
                idBits = idLine.split(' ')
                standardizedUMI = savedRead1[:6] #Sometimes this overlaps with the first letter of the barcode if there is a read error where one of the UMI letters is missing
                idBits[0] += "_" + curBC + "_" + standardizedUMI
                idBits = " ".join(idBits)
                blockbuilder += idBits
                # then we add in the stripped DNA
                blockbuilder += rna_to_write
                # the plus
                blockbuilder += ('+' + '\n')
                # and the stripped QC scores
                blockbuilder += trimmed_QC
                # write and reset
                outfile.write(blockbuilder)
                blockbuilder = curBC = rna_to_write = ""
                # Skip the next two lines in the file cause we don't care about this read aside from the UMI
                next(f)
                next(f)
                continue
            #We have reached the end of the file, write and break.
            if line == '':
                #outfile.write(blockbuilder)
                print("UNHANDLED CODE")
                break
            if nextLineDNA:
                savedRead1 = line
                nextLineDNA = False
                totalReads+=1
                curBC = rna_to_write = trimmed_QC = ""
                pattern_fuzzy = "(" + handle + "){e<=1}" #allow for edit distance of 1
                # without bestmatch, it will match other nucleotides even if the handle is perfect
                match = regex.search(pattern_fuzzy, line, regex.BESTMATCH)
                if match is None:
                    noHandles+=1
                    noHandleList.append(line)
                    #skip to start of next read 1 and ignore rest
                    SkipToNextRead(f)
                    continue
                handle_start = match.regs[0][0]
                handle_end = match.regs[0][1]
                if handle_start > MAX_HANDLE_POS:
                    falseHandles+=1
                    SkipToNextRead(f)
                    continue
                BC1 = line[(handle_start-BClen):handle_start]
                BC2 = line[handle_end:(handle_end + BClen)]
                # Verify Cell BC
                #curBC = DebugBadBarcodes(BC1, BC2, line) #takes longer but outputs bad reads into a file
                curBC = CorrectCellBarcodes(BC1, BC2)
                if curBC is None:
                    faultyBarcodes+=1
                    SkipToNextRead(f)
                    continue
                #umi = line[:(handle_start-BClen)]
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

        # close all the file writers
        for writer in sample_filewriters.values():
            writer.close()

    with open("noHandleReads.txt", 'w') as noHandle_file:
        for read in noHandleList:
            noHandle_file.write(read)

    print(fq + " processing complete:", flush=True)
    print("\t-" + str(falseHandles) + " False Handles detected.")
    print("\t-" + str(noHandles) + " reads had no valid handle.")
    print("\t-" + str(faultyBarcodes) + " reads with at least one faulty barcodes")
    print("\t-" + str(faultySampleBC) + " reads with a bad sample BC")
    print("\t-" + str(totalReads) + " total reads")
    # Reads are discarded in this order and are cumulative, so they are output in this order.
    badHandle_percent = round((noHandles / totalReads) * 100, 1)
    badBC_percent = round((faultyBarcodes / totalReads) * 100, 1) #rounds to nearest tenth
    badSample_percent = round((faultySampleBC / totalReads) * 100, 1)
    print("\t-" + "Bad Handle Percentage: " + str(badHandle_percent) + "%")
    print("\t-" + "Bad Cell BC Percentage: " + str(badBC_percent) + "%")
    print("\t-" + "Bad Sample BC percentage: " + str(badSample_percent) + "%")
    print("\t-" + "Total read percentage discarded: " + str(badBC_percent + badHandle_percent + badSample_percent) + "%")

# A function that writes bad barcodes and reads to a file for debug
# This function significantly increases runtime, when bc are not being debugged please opt to use normal "CorrectCellBarcodes" instead.
def DebugBadBarcodes(bc1, bc2, read):
    with open("badBarcodes.txt", 'a') as badBarcode_file:
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
        pattern = "(" + bc + "){e<=1}"
        for cellBC in cellBC_codes:
            match = regex.search(pattern, cellBC, regex.BESTMATCH)
            if match is not None: return cellBC
        return None

def CheckSampleBarcode(bc):
    if bc in sampleBC_to_ID.keys():
        return bc
    else:
        pattern = "(" + bc + "){e<=1}"
        for sampleBCs in sampleBC_to_ID.keys():
            match = regex.search(pattern, sampleBCs, regex.BESTMATCH)
            if match is not None: return sampleBCs #originally this returned match.group(1) but now i'm thinking that
            # all checks should return the og sampleBC rather than the match.group(1) because in a case like this one of the
            # nucleotides was missing. Now this is putting the other methods that use match.group(1) into doubt for me
            # Ended up changing as the barcode was indeed being set to 7 nucleotides rather than the correct 8
        return None

def SkipToNextRead(file):
    try:
        for i in range(6):
            next(file)
    except StopIteration as e:
        print(e)

# Generates a list of writers, one for each sample present in a file
# Requires file names to be formatted as indicated at top of the script
def GetFileWriters(sampleFileList):
    sampleFileWriterz = {}
    sampleList = []
    sampling = True
    sampling_index = 1
    while sampling:
        if sampleFileList[sampling_index].isnumeric():
            sampleList.append(sampleFileList[sampling_index])
            sampling_index+=1
        else:
            sampling = False
            laneID = sampleFileList[sampling_index]
            for sample in sampleList:
                name = (sample + "_" + laneID)
                sampleFileWriterz[sample] = gzip.open((outputFolder + name + "_processed.fq.gz"), 'wt')
    return sampleFileWriterz


pool = ThreadPool(processes=MAX_THREADS)

for fq in fastq_files:
    #add the files to the multithreading handler
    result = pool.apply_async(Processfile, args=(fq, 1))
    #print(result.successful())

pool.close()
pool.join()

#print("All files processed!", flush=True)

#Processfile("OL_04_05_06_D1_merged.fq.gz", "OL456_D1_processed.fq.gz")
#Processfile("OL_04_05_06_test.fq.gz", "456_test_processed.fq.gz")