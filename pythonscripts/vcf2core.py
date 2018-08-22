#!/usr/bin/python

usage = """
This script is for transformating vcf file of diploid samples to corehunter input.
It was written by Dang Liu. Last updated: Jan 22 2018. 
usage:
python vcf2core.py vcf

"""

# modules here
import sys, re

# if input number is incorrect, print the usage
if len(sys.argv) < 2:
	print usage
	sys.exit()

# input vcf file
vcf = open(sys.argv[1], 'r')


# collections
geno_dit = {'0/0':'0', '1/0':'1', '0/1':'1', '1/1':'2', './.':' '}
sample_dit = {}
k = 0
# processing here
line = vcf.readline()
while(line):
	# get rid of comments
	if ('##' in line):
		pass
	# record sample id
	elif ('#' in line):
		line_s = re.split(r'\s+', line.replace('#', '').replace('\n+', ''))
		l = len(line_s)
		sample_id = line_s[9:l-1]
		n = 9
		for i in sample_id:
			sample_dit[n] = []
			sample_dit[n].append(i)
			n += 1
		#print sample_id
	# record genotype
	else:
		k += 1
		print 'process ' + str(k) + ' loci...'
		line_s = re.split(r'\s+', line.replace('\n+', ''))
		l = len(line_s)
		genotype = line_s[9:l-1]
		n = 9
		for i in genotype:
			geno = re.sub(r':[-:,.0-9]*', '', i)
			sample_dit[n].append(geno_dit[geno])
			n += 1
	line = vcf.readline()
vcf.close()
print 'loci all collected!'

# output here
out_f = open(sys.argv[1].replace('vcf', 'core.txt'), 'w')
out_f.write('ID\tNAME\t')
for i in range(1, k+1):
	if (i != k):
		out_f.write('mk' + str(i) + '\t')
	if (i == k):
		out_f.write('mk' + str(i) + '\n')
for i in sample_dit:
	out_f.write(str(i-8) + '\t')
	check_k = 0
	for g in sample_dit[i]:
		if (check_k != k):
			out_f.write(g + '\t')
		else:
			out_f.write(g + '\n')
		check_k += 1

print sys.argv[1].replace('vcf', 'core.txt') + ' is output.'
print 'All done!'

#print sample_dit


#last_v20180122
