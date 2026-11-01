#!/bin/env python

import argparse
import re
def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f','--file',help='your input SAM file',required=True)
    parser.add_argument('-o','--outfile',help='your output SAM file',required=True)
    parser.add_argument('-u','--UMI',help='your known UMI input file',required=True)
    return parser.parse_args()
args = get_args()

def get_UMI(QNAME):
     '''
     This function returns the Unique Molecular Identifiers associated with a read in a SAM file.
     '''
     UMI = QNAME.split(":")[-1]
     return UMI
    
def stradedness(BIT):
    '''
    This function returns whether a read is forward or reverse stranded
    True = forward
    False = Reverse
    '''
    if (BIT & 16) == 16: 
        return False
    else:
        return True  

def true_start(CIGAR,POS,BIT):
    '''
    This function finds the "true" 5' starting position of a read from a SAM file.
    '''
    start = POS
    fwd = stradedness(BIT)
    cig_list = re.findall(r"[0-9]+|[A-Z]+",CIGAR) # splits the cigar string into list ['number','character]
    rvs_5prime = 0
    prime3_clip = False
    if fwd: # if forward strand
          for i in cig_list:
            match = re.search(r"[0-9]+", i) # checks to see if item in list is digit
            if match:
                   length = int(i)
            elif i in "MDN": # ignore if character is a M,D, or N
                break
            else:
                if i == "S": # if it is an S, subtract to get the true 5' start
                      start = start - length
                length = 0
    else: # if reverse strand
        for list_pos,i in enumerate(cig_list):
            match = re.search(r"[0-9]+", i) # checks to see if item in list is digit
            if match:
            	length = int(i)
            elif i == "S" and list_pos == 1: # if there is an S at the beginning of the cigar string
                continue
            elif i in "MDNSH": # if the character is one of these characters, add the length.
                start += length
    return start

seen_set = set()
chrom = ''
seen_count = 0
fwd_seen_count = 0
rvs_seen_count = 0
chrom_count = 0
total_count = 0
wrong_UMI = 0

with open(args.file, "r") as i:
    with open(args.outfile, "w") as o:
        for line in i: # for each line in the SAM file
            if line.startswith("@") or line.split("\t")[0] == "QNAME": # if the line starts with "@" (header line) or "QNAME" (part of my test file), add it to the output file
                 o.write(f"{line}")
                 continue
            QNAME = line.split("\t")[0] #extracts the read title with the UMI
            BIT = int(line.split("\t")[1]) #extracts the bitflag of the read for strandedness
            RNAME = line.split("\t")[2] #extracts the chromosome name
            if chrom != '' and chrom != RNAME: #if the previous and current chromosomes are not equal to each other, empty the set and print how many of the previous chromosome reads were excluded.
                seen_set = set()
                print(f"{chrom}\t{chrom_count}")
                chrom_count = 0 #reset exclusion count
            POS = int(line.split("\t")[3]) #extracts the leftmost mapping position
            CIGAR = line.split("\t")[5] #extracts the cigar 
            five_prime = true_start(CIGAR,POS,BIT) # grabs the 5prime pos.
            UMI = get_UMI(QNAME) #gets the UMI from QNAME
            identifier = (five_prime,get_UMI(QNAME),RNAME, stradedness(BIT)) # adds the unique tags associated with a PCR duplicate into a set
            if "N" in UMI: # if the UMI contains an N, add the wrong UMI, and continue
                wrong_UMI += 1
                continue
            if identifier in seen_set: # If these unique tags have already been seen, don't add to the file
                fwd = stradedness(BIT)
                seen_count += 1
                if fwd: # if forward strand
                    fwd_seen_count +=1
                else: # if reverse strand
                    rvs_seen_count +=1
                continue
            else: # write the line to the output file
                 o.write(f"{line}")
                 total_count += 1
                 chrom_count += 1
            seen_set.add(identifier)
            chrom = RNAME
# print statements for results of deduplicication
print(f"chromosome {chrom} sequences:\n{chrom_count}")
print(f"total reads in file {total_count}")
print(f"total # of reads excluded:\ntotal: {seen_count}\nfwd: {fwd_seen_count}\nrevervse:{rvs_seen_count}") 
print(f"Invalid UMIs: {wrong_UMI}")

# testing cigar string
assert true_start("74M", 76718921, 0) == 76718921
assert true_start("3S71M", 76718924, 0) == 76718921
assert true_start("1S71M", 3123097, 0) == 3123096
assert true_start("1S71M", 3123097, 16) == 3123168
assert true_start("41M53N31M", 3147363, 0) == 3147363
assert true_start("41M53N31M", 3147363, 16) == 3147488
assert true_start("68M4S", 3369137, 16) == 3369209
assert true_start("5S68M4S", 3369137, 16) == 3369209
assert true_start("68M7I4S", 3369137, 16) == 3369209
