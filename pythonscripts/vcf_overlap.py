#!/usr/bin/python3

usage = """
This script is for comparing the overlapping of two vcf.
It was written by Dang Liu. Last updated: May 29 2018.
usage:
python3 vcf_overlap.py vcf1 vcf2

"""

# modules here
import sys, re

# if input number is incorrect, print the usage
if len(sys.argv) < 3:
	print(usage)
	sys.exit()



# Read vcf1
vcf1_dit = {}
vcf1 = open(sys.argv[1], 'r')
line1 = vcf1.readline()
while(line1):
	if ("#" not in line1):
		line1_s = re.split(r'\s+', line1)
		chro = line1_s[0]
		pos = line1_s[1]
		ref = line1_s[3]
		if (chro not in vcf1_dit):
			vcf1_dit[chro] = {}
			vcf1_dit[chro][pos] = ref
		else:
			vcf1_dit[chro][pos] = ref
	line1 = vcf1.readline()
# close barcode file
vcf1.close()


# Read vcf2 
overlap = 0
incons = 0
vcf2 = open(sys.argv[2], 'r')
line2 = vcf2.readline()
while(line2):
	if ("#" not in line2):
		line2_s = re.split(r'\s+', line2)
		chro2 = line2_s[0]
		pos2 = line2_s[1]
		ref2 = line2_s[3]
		if (chro2 in vcf1_dit):
			if (pos2 in vcf1_dit[chro2]):
				if (ref2 == vcf1_dit[chro2][pos2]):
					print("Found one overlap!")
					overlap += 1
				else:
					print("There is an inconsistency between REFs! {0}* in vcf1 and {1}* in vcf2.".format(vcf1_dit[chro2][pos2], ref2))
					incons += 1
	line2 = vcf2.readline()
vcf2.close()

print("All done! There are {0:.0f} overlaps and {1:.0f} inconsistencies.".format(overlap, incons))

# last_v20180529