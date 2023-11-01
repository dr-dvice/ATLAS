# Author: David Bellini
# Date: 11.1.23
# A script to conduct all necessary preprocessing on raw single cell data for the cell atlas, total RNA format.
import os
import sys
from os import listdir
from os.path import isfile, join
import gzip
import regex
from multiprocessing.pool import ThreadPool
#import traceback
#from line_profiler_pycharm import profile
import argparse

parser = argparse.ArgumentParser(
    description="A demultiplexing script using regex for careful barcode correction of total RNA matq-drop experiments"
)
parser.add_argument("fastqgz_directory", type=str, help="the filepath for the directory containing the fastq.gz files to be" \
                                                       "processed")
parser.add_argument("barcode_file", type=str, help="The full path to the barcode file containing 2 columns, one for barcode one and"
                                                   "one for barcode 2")
parser.add_argument("fastqgz_samplestring", type=str, help="The identifying sample characters of the samples. For ex: If I have 3 samples like so:\n"
                                                           "\tGS_D01_CKDL230034544-1A_22F222LT3_L8.fq.gz\n"
                                                           "\tGS_D02_CKDL230034544-1A_22F222LT3_L8.fq.gz\n"
                                                           "\tGS_D03_CKDL230034544-1A_22F222LT3_L8.fq.gz\n"
                                                           "Then the identifying sample characters would be 'GS_D' - 01,02, and 03 will be selected and used to assign"
                                                           "synthetic sample barcodes in the code")
parser.add_argument("-a", "--handle", type=str, help="The handle that should be selected. This algorithm uses the handle to determine the location"
                                             "of all other RNA features (bc, umi, etc)", default="TGGTAGGTGGTAGAGA")
parser.add_argument("-o", "--output", type=str, help="output directory", default="results/")
parser.add_argument("-t", "--threads", type=int, help="number of threads to be used. default is 1. Recommended to use as many threads"
                                                      "as there are files, if possible", default=1)
parser.add_argument("-b", "--barcodelen", type=int, help="length of the barcodes. default is 10", default=10)
parser.add_argument("-u", "--umilen", type=int, help="length of the umis. default is 8", default=8)
parser.add_argument("-m", "--maxhandlepos", type=int, help="the furthest position at which you would expect to find the barcode. "
                                                           "debugMode can be useful to find this if uknown.", default=45)
parser.add_argument("-d", "--debugmode", action=argparse.BooleanOptionalAction, help="calculates and outputs additional statistics, "
                                                                                     "such as handle position distribution, and bad barcodes")

args = vars(parser.parse_args())

fastqgz_directory = args["fastqgz_directory"]
outputFolder = args["output"]
cellBCfileName = args["barcode_file"]
samplestr = args["fastqgz_samplestring"]



# hardcoded shortcuts
#fastqgz_directory = "/home/davidbellini/files/workspace/sg_ol2_analysis/fastq_files_merged_aatrimmed/"
#cellBCfileName = "/home/davidbellini/files/workspace/sg_ol2_analysis/MATQ_Drop_384x384_barcode.tsv"
#outputFolder = "/home/davidbellini/files/workspace/sg_ol2_analysis/results/"


# Parameter Vars
MAX_THREADS = args["threads"]
BClen = args["barcodelen"] #expected length of the barcodes
DIVIDER = BClen + 2
UMIlen = args["umilen"]
MAX_HANDLE_POS = args["maxhandlepos"] #Handle is not expected to the farther than 16 bases from the start
# debugmode and single column mode require fixes if they are to be used in the future
debugMode = args["debugmode"]
#handle = args["handle"]
handle = "TGGTAGGTGGTAGAGA"
handlePattern = "(" + handle + "){e<=3}"
sampleBC_to_ID = {}
cellBC_codes = []


def Read2ColBarcodes(filename):
    barcodes1 = []
    barcodes2 = []
    with open (filename) as file:
        for line in file:
            data = line.strip().split("\t")
            barcodes1.append(data[0])
            barcodes2.append(data[1])
    return barcodes1, barcodes2

print("Reading Cell BC file...", flush=True)


cellBC_list_col1, cellBC_list_col2 = Read2ColBarcodes(cellBCfileName)
cellBC_str_col1 = "..".join(cellBC_list_col1)
cellBC_str_col2 = "..".join(cellBC_list_col2)

print("Cell BC file read", flush=True)

fastq_files = [f for f in listdir(fastqgz_directory) if isfile(join(fastqgz_directory, f))]

def Processfile(fq, synth_code):
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
    BCmutationstats = [0,0,0] #subsitutions insertions deletions

    # read content from original fastq file and prep for processing
    with gzip.open(fq_path, 'rt') as f:
        name = fq.split(".")[0]
        if not os.path.exists(outputFolder) : os.makedirs(outputFolder)
        outfile = gzip.open((outputFolder + name + "_processed.fq.gz"), 'wt')
        read2 = False
        rna_to_write = ""
        QC_line = ""
        idLine = ""
        nextLineDNA = False
        end_of_file = False
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
                    if debugMode: cBarcode = DebugBadBarcodes(barcode1, barcode2, line, BCmutationstats)
                    else: cBarcode = CorrectCellBarcodes2col(barcode1, barcode2, line, BCmutationstats)
                    if cBarcode is None:
                        faultyBarcodes+=1
                        SkipToNextRead(f)
                        ResetVars(idLine, rna_to_write, QC_line)
                        continue
                    else: # Checks successfull, here we actually take the read and the umi and create the new ID line for read 1.
                        # I use to also have a + "_" between synth code and cBarcode but that just makes thing more complicated furthe down the line
                        newIDline = "_" + synth_code + cBarcode + "_" + line[UMIStartpos:(UMIStartpos + UMIlen)]
                        idBits = idLine.split(' ')
                        idBits[0] += newIDline
                        idBits = " ".join(idBits)
                        # Time to jot everything down
                        blockbuilder += idBits  + rna_to_write + "+"  + "\n" + QC_line

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
                    outfile.write(blockbuilder) #WRITE OUT EVERY MILLION LINES. Hopefully the /n are working
                    blockbuilder=""
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
        outfile.write(blockbuilder) #write out what's remanining in the blockbuilder
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
    print("\t-" + str(sum(BCmutationstats)) + " Barcodes were corrected:")
    print("\t\t-" + str(BCmutationstats[0]) + " subsitutions")
    print("\t\t-" + str(BCmutationstats[1]) + " (Insertions not tested)")
    print("\t\t-" + str(BCmutationstats[2]) + " deletions")
    print("\t-" + str(totalReads) + " total reads")
    # Reads are discarded in this order and are cumulative, so they are output in this order.

    badHandle_percent = round(((noHandles + falseHandlesNum) / totalReads) * 100, 1)
    read2short_percent = round((read2short / totalReads) * 100, 1) #rounds to nearest tenth
    faultyBarcodes_percent = round((faultyBarcodes / totalReads) * 100, 1)
    correctedBarcodes_percent = round((sum(BCmutationstats) / totalReads) * 100, 1)
    print("\t-" + "Bad Handle Percentage: " + str(badHandle_percent) + "%")
    print("\t-" + "Short Read Percentage: " + str(read2short_percent) + "%")
    print("\t-" + "Corrected Barcode percentage: " + str(correctedBarcodes_percent) + "%")
    print("\t-" + "Discarded Cell BC Percentage: " + str(faultyBarcodes_percent) + "%")
    print("\t-" + "Total Read Percentage Discarded: " + str(round((read2short_percent + badHandle_percent + faultyBarcodes_percent), 1 )) + "%")

# A function that writes bad barcodes and reads to a file for debug
# This function significantly increases runtime, when bc are not being debugged please opt to use normal "CorrectCellBarcodes" instead.
def DebugBadBarcodes(bc1, bc2, line, BCmutationstats_list):
    with open((outputFolder + "badBarcodes.txt"), 'a') as badBarcode_file:
        BC1_check = CheckCellBarcode2col(bc1, 1, line, BCmutationstats_list)
        BC2_check = CheckCellBarcode2col(bc2, 2, line, BCmutationstats_list)
        if BC1_check is None:
            badBarcode_file.write("BC1 is bad: " + bc1 + " ")
        if BC2_check is None:
            badBarcode_file.write("BC2 is bad: " + bc2 + " ")
        if None in [BC1_check, BC2_check]:
            badBarcode_file.write("BAD read: " + line + '\n')
            return None
        else: return BC1_check + BC2_check

#@profile
def CorrectCellBarcodes2col(bc1, bc2, line, BCmutationstats_list):
    BC1_check = CheckCellBarcode2col(bc1, 1, line, BCmutationstats_list)
    BC2_check = CheckCellBarcode2col(bc2, 2, line, BCmutationstats_list)
    if None in [BC1_check, BC2_check]:
        return None
    else: return BC1_check + BC2_check

#@profile
def CheckCellBarcode2col(bc, col, line, BCmutationstats_list):
    cellBCs = cellBC_list_col1 if col == 1 else cellBC_list_col2
    if bc in cellBCs:
        return bc
    else:
        pattern = "(" + bc + "){e<=1}"
        cellBC_str = cellBC_str_col1 if col == 1 else cellBC_str_col2
        match = regex.search(pattern, cellBC_str, regex.BESTMATCH) #Fuzzy_counts=(subsitutions, insertions, deletions)
        if match is None: return None
        else:
            # If subsitution match
            matchedBC = match.group()
            if (match.span()[1] - match.span()[0]) == BClen and ("." not in match.group()):
                BCmutationstats_list[0]+=1 #In this situation, has to be a subsitutions
                return matchedBC # In reality, i have to do a lot more checking than this and find the right BC
            else: #Its a deletion and we need to work even more to match the right barcode
                BCmutationstats_list[2]+=1
                #Lets use math. If every barcode covers 12 chars in the string, we should be able to use the modulo
                #operator to extract the right barcode.
                barcodeNum = match.span()[1] // DIVIDER
                barcodeLoc = barcodeNum * DIVIDER
                barcode = cellBC_str[barcodeLoc:(barcodeLoc+BClen)]
                return barcode

def SkipToNextRead(file):
    try:
        for i in range(2): #changed to 2 from 6 of dt
            next(file)
    except StopIteration as e:
        print(e)

def ResetVars(id, rna, QC):
    id = rna = QC = ""

def GenerateSyntheticBarcodes(amount):
    nucleotides = ['A', 'C', 'G', 'T']
    sequences = []
    for i in range(amount):
        sequence = ''
        for j in range(3):
            sequence += nucleotides[i % 4]
            i //= 4
        sequences.append(sequence)

    return sequences

syntheticBarcodes = GenerateSyntheticBarcodes(len(fastq_files))

pool = ThreadPool(processes=MAX_THREADS)

for fq in fastq_files:
    if fq == ".DS_Store": continue
    id_str = samplestr + "(\d+)"
    sample_ID_match = regex.search(id_str, fq)
    sample_ID = sample_ID_match.group(1)
    #add the files to the multithreading handler
    try:
        result = pool.apply_async(Processfile, args=(fq, syntheticBarcodes[int(sample_ID)-1]))
    except Exception as e:
        print(e)

    #print(result.successful())
    print(sample_ID + " :" + "was assigned synthetic code: " + syntheticBarcodes[int(sample_ID)-1])


pool.close()
pool.join()

print("All files processed!", flush=True)

Processfile("test.fq.gz", "AAA")
