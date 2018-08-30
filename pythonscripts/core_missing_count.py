#!/usr/bin/python3

usage = """

A script to count allele missing rate from core.txt (output from vcf2core).
This script is written by Dang Liu. Last updated: Aug. 30 2018. 
usage:
python3 core_missing_count.py core.txt

"""

import sys, re

if len(sys.argv) < 2:
	print(usage)
	sys.exit()


f = open(sys.argv[1], 'rt')

line = f.readline()
total_allele_n = len(line.replace("\n", "").split('\t')) - 2
line = f.readline()
while(line):
	line_s = line.replace("\n", "").split("\t")
	sample = line_s[1]
	missing = line_s.count(" ")
	missing_rate = (missing/total_allele_n)*100
	print(sample + "\t" + str(missing) + "\t{0:.2f}".format(missing_rate))
	line = f.readline()

f.close()


# done

#last_v20180830