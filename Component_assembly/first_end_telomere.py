#!/usr/bin/python

#Usage: python first_end_telomere.py inputfile outputfile first_coordinate_telomere_subtelomereboundary 


import sys

inputfile = sys.argv[1]
outputfile =  sys.argv[2]
#first coordinate of the telomere-subtelomere boundary
i = int(sys.argv[3]) 
sys.stdout=open(outputfile,"w")
my_file = open(inputfile)

file_contents = my_file.read()
first_kb = file_contents[0:i]
telomere_count = first_kb.count("CCCTAAA")
print (telomere_count)
