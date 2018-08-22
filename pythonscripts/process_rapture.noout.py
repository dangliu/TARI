#!/usr/bin/python3

usage = """
This script is for filtering/rearranging/demultiplexing Rapture raw reads.
It was written by Dang Liu. Last updated: Mar 22 2018.
usage:
python3 process_rapture.noout.py barcodes read1 read2 enz_cutsite

"""

# modules here
import sys, re
import gzip
import io
import time

# if input number is incorrect, print the usage
if len(sys.argv) < 5:
	print(usage)
	sys.exit()


# time the process start here
tStart =time.time()


# restriction enzyme cutsite
enz_cutsite = sys.argv[4]
enz_cutsite_len = len(enz_cutsite)

# input barcode file
# format: barcode"\t"sample_ID
file1 = open(sys.argv[1], 'r')
line1 = file1.readline()
barcode_length = len(line1.split("\t")[0])
sample_l = []
# collect barcode info.
barc_dict = {}
while(line1):
	barcode = line1.split("\t")[0]
	sample_ID = line1.split("\t")[1].replace("\n", "")
	barc_dict[barcode] = sample_ID
	sample_l.append(sample_ID)
	line1 = file1.readline()
# close barcode file
file1.close()

# Read read1 and read2; files can be either gz or not
if (".gz" in sys.argv[2]):
	file2 = gzip.open(sys.argv[2], 'rt')
	#file2 = io.BufferedReader(gz2)
	file3 = gzip.open(sys.argv[3], 'rt')
	#file3 = io.BufferedReader(gz3)
else:
	file2 = open(sys.argv[2], 'r')
	file3 = open(sys.argv[3], 'r')


enz_cutsite_n = 0

header1 = file2.readline()
seq1 = file2.readline()
strand1 = file2.readline()
q1 = file2.readline()
header2 = file3.readline()
seq2 = file3.readline()
strand2 = file3.readline()
q2 = file3.readline()
reads_n = 0
retained_n = 0
double_n = 0
not_found_n = 0
while(header1):
	reads_n += 2
	#print("R1: " + header1.replace("\n", ""))
	#print("R1: " + seq1.replace("\n", ""))
	#print("R1: " + strand1.replace("\n", ""))
	#print("R1: " + q1.replace("\n", ""))
	#print("R2: " + header2.replace("\n", ""))
	#print("R2: " + seq2.replace("\n", ""))
	#print("R2: " + strand2.replace("\n", ""))
	#print("R2: " + q2.replace("\n", ""))
	if (seq1[:barcode_length] in barc_dict and seq2[:barcode_length] not in barc_dict):
		print("Process " + barc_dict[seq1[:barcode_length]] + "......")
		retained_n += 2
		if (seq1[barcode_length:barcode_length+enz_cutsite_len] == enz_cutsite):
			enz_cutsite_n += 2
	elif (seq1[:barcode_length] not in barc_dict and seq2[:barcode_length] in barc_dict):
		retained_n += 2
		print("Process " + barc_dict[seq2[:barcode_length]] + "......")
		if (seq2[barcode_length:barcode_length+enz_cutsite_len] == enz_cutsite):
			enz_cutsite_n += 2		
	elif (seq1[:barcode_length] in barc_dict and seq2[:barcode_length] in barc_dict):
		double_n += 2
		print("Double barcodes in %s!! %s and %s." % (re.split(r'\s+', header1)[0], seq1[:barcode_length], seq2[:barcode_length]))
		if (seq1[barcode_length:barcode_length+enz_cutsite_len] == enz_cutsite or seq2[barcode_length:barcode_length+enz_cutsite_len] == enz_cutsite):
			enz_cutsite_n += 2
	else:
		not_found_n += 2
		# Try to find out why there were no barcodes
		print("Barcodes not found in %s, which begin as %s and %s." % (re.split(r'\s+', header1)[0], seq1[:barcode_length], seq2[:barcode_length]))        
	header1 = file2.readline()
	seq1 = file2.readline()
	strand1 = file2.readline()
	q1 = file2.readline()
	header2 = file3.readline()
	seq2 = file3.readline()
	strand2 = file3.readline()
	q2 = file3.readline()
	

print("All done!")
# Add stats here
# use .format instead of %d (for digital), %s (for string) and %.1f (for float with one number after decimal)
# {number} is for the index in .format(list) 
print("-------Stats--------")
print("Process {0} reads, retained {1} reads ({2:.2f}%; {3:.0f} pairs).".format(reads_n, retained_n, retained_n*100/reads_n, retained_n/2))
print("There are {0:.0f} read pairs without restriction enzyme 1 cutting site ({1:.2f}%).".format((reads_n-enz_cutsite_n)/2, (100 - (enz_cutsite_n*100/reads_n))))
print("Filter out {0} reads ({1:.2f}%) with double barcodes.".format(double_n, double_n*100/reads_n))
print("Filter out {0} reads ({1:.2f}%) with barcodes not found.".format(not_found_n, not_found_n*100/reads_n))
print("--------------------")



# close read files
#if (".gz" in sys.argv[2]):
#    gz2.close()
#    gz3.close()
#else:
file2.close()
file3.close()


# time the process end here
tEnd = time.time()
print("It cost %f sec." % (tEnd - tStart))

# last_v20180530
