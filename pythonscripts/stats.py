#!/usr/bin/python

#genome_size [\t] N_of_scf [\t] Largest_scf [\t] N50 [\t] L50 [\t] N90 [\t] L90 [\t] valid_nuc [\t] gap
usage = """
This script is written by Dang Liu. Last updated: Oct 03 2017. 
usage:
python stats.py genomic.fa

"""

import sys, re

if len(sys.argv) < 2:
	print usage
	sys.exit()

# define functions here
def stats(numlist):
	sorted_list = sorted(numlist, reverse=True)
	genome_size = sum(sorted_list)
	N_of_scf = len(sorted_list)
	Largest_scf = sorted_list[0]
	sum_i = 0
	check_N50 = 0
	check_N90 = 0
	for i in sorted_list:
		sum_i += i
		if (sum_i >= genome_size/2 and check_N50 == 0):
			N50 = i
			L50 = sorted_list.index(i) + 1
			check_N50 = 1
		elif (sum_i >= genome_size*9/10 and check_N90 == 0):
			N90 = i
			L90 = sorted_list.index(i) + 1
			check_N90 = 1
		else:
			pass
	final_l = [str(genome_size), str(N_of_scf), str(Largest_scf), str(N50), str(L50), str(N90), str(L90)]			
	return final_l

def killo(n):
	kn = "%.2f" % (float(n)/1000.0)
	return kn

def mega(n):
	mn = "%.2f" % (float(n)/1000000.0)
	return mn

# A dcit to collect sequence length
seq_dit = {}

# Count gaps by counting "N"
gap = 0

f = open(sys.argv[1], 'r')
line = f.readline()

while(line):
        if (">" in line):
                ID = line.replace(">", "").replace("\n", "")
                seq_dit[ID] = 0
        else:
                seq = line.replace("\n", "")
                N = seq.count("N")
                seq_len = len(seq)
                gap += N
                seq_dit[ID] += seq_len
        line = f.readline()

f.close()

# A list for stats calculation
len_list = []

for i in seq_dit:
	print i + "\t" + str(seq_dit[i])
	len_list.append(seq_dit[i])

final_l = stats(len_list)
valid_nuc = int(final_l[0]) - gap

final_l.append(str(valid_nuc))
final_l.append(str(gap))

# output
print "%10s %5s %10s %10s %5s %10s %5s %10s %10s" % ("Genom_size", "N.Scf", "Large_scf", "N50", "L50", "N90", "L90", "Valid_nuc", "Gap")
print "%10s %5s %10s %10s %5s %10s %5s %10s %10s" % (killo(final_l[0]), final_l[1], killo(final_l[2]), killo(final_l[3]), final_l[4], killo(final_l[5]), final_l[6], killo(final_l[7]), killo(final_l[8]))
print "%10s %5s %10s %10s %5s %10s %5s %10s %10s" % (mega(final_l[0]), final_l[1], mega(final_l[2]), mega(final_l[3]), final_l[4], mega(final_l[5]), final_l[6], mega(final_l[7]), mega(final_l[8]))
#last_v20171003

