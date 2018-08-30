#!/usr/bin/python3

usage = """

A script to extract 'QD', 'FS', 'SOR', 'MQ', 'MQRankSum', 'ReadPosRankSum' from vcf INFO field.
This script is written by Dang Liu. Last updated: Aug. 28 2018. 
usage:
python3 vcf_extract.py vcf

"""

import sys, re
import gzip

if len(sys.argv) < 2:
	print(usage)
	sys.exit()

out_f = open(sys.argv[1] + "_out.txt", 'wt')
# make a list to collect info.
INFO_list = ['QD', 'FS', 'SOR', 'MQ', 'MQRankSum', 'ReadPosRankSum']
out_f.write("CHROM\tPOS\tFILTER\tQD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum\n")
#print("QD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum")

n = 1

# Read vcf by each line

# files can be either gz or not
if (".gz" in sys.argv[1]):
	f = gzip.open(sys.argv[1], 'rt')
else:
	f = open(sys.argv[1], 'rt')
line = f.readline()
while(line):
	# Ignore comments
	if (re.match(r'#', line)):
		pass
	else:
		print("Extracting " + str(n) + "th varinat...")
		# split by space and extract info
		line_split = re.split(r'\s+',line)
		CHROM = line_split[0]
		POS = line_split[1]
		FILTER = line_split[6]
		out_f.write(CHROM + "\t" + POS + "\t" + FILTER + "\t")
		INFO = line_split[7]
		# loop by INFO_list
		for i in INFO_list:
			# build a regex as string including the element from INFO_list
			# MQ=54.97;MQRankSum=-1.645e+00;QD=32.09;ReadPosRankSum=0.00;SOR=3.953
			# () can be used to extract the number in group()
			# Use [-+e.0-9]+ instead of [-+e.0-9]* to avoid empty match written into output
			my_regex = re.escape(i) + r'=([-+e.0-9]+)' 
			#print(my_regex)
			#print(re.search(my_regex, line))
			#print(re.search(my_regex, line).group(1))
			# search my_regex pattern in each line
			if (re.search(my_regex, line)):
				#number=re.search(k + r'=(\d+);', line).group(1)
				#INFO_dict[k].append(number)
				# () can be used to extract the number in group()
				number=re.search(my_regex, line).group(1)				
				out_f.write(number)
				#out += number + "\t"
				# if the element is the last one in list write out \n 
				if (INFO_list.index(i) == 5):
					out_f.write("\n")
				else:
					out_f.write("\t")
			else:
				#INFO_dict[k].append('NA')
				out_f.write("NA")
				#out += "NA\t"
				if (INFO_list.index(i) == 5):
					out_f.write("\n")
				else:
					out_f.write("\t")
		#print(out)	
		n += 1			
	line = f.readline()

f.close()
out_f.close()


# done

# last_v20180829