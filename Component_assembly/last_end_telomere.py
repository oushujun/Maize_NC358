#!/usr/bin/python

python first_end_telomere.py inputfile outputfile last_coordinate_telomere_subtelomereboundary
  

import sys

inputfile = sys.argv[1]

outputfile =  sys.argv[2]
#last coordinate of the telomere-subtelomere boundary
i = int(sys.argv[3])
sys.stdout=open(outputfile,"w")

my_file = open(inputfile)
file_contents = my_file.read()
length = int(len(file_contents))
last_kb = file_contents[length-i:length]
telomere_count = last_kb.count("TTTAGGG")
print (telomere_count)
