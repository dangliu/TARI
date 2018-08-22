#!/usr/bin/python3

usage = """
This script is for renaming sample ID in vcf according to smaple_ID.csv.
It was written by Dang Liu. Last updated: May 28 2018.
usage:
python3 process_rapture.py sample_ID.csv vcf

"""

# modules here
import sys, re

# if input number is incorrect, print the usage
if len(sys.argv) < 3:
	print(usage)
	sys.exit()



# input sample_ID.csv file
# format: barcode_ID","sample_ID
file1 = open(sys.argv[1], 'r')
line1 = file1.readline()
# collect info.
sample_dict = {}
while(line1):
	barcode_ID = line1.split(",")[0]
	sample_ID = line1.split(",")[1].replace("\n", "")
	sample_dict[barcode_ID] = sample_ID
	line1 = file1.readline()
# close barcode file
file1.close()


# Read vcf 
file2 = open(sys.argv[2], 'r')
line2 = file2.readline()
while(line2):
	if (re.search("#CHROM", line2)):
		line2_split = re.split(r'\s+', line2.replace("\n", ""))
		new_line2 = line2_split[0]
		for i in line2_split[1:]:
			if (line2_split.index(i) < 9):
				new_line2 += "\t" + i
			else:
				new_line2 += "\t" + sample_dict[i]
		print(new_line2)
	else:
		print(line2.replace("\n", ""))        
	line2 = file2.readline()

file2.close()

# last_v20180528