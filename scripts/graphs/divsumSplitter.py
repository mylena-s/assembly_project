import sys
import csv

def divsumSplitter_improved(file):
	with open(file, 'r') as infile, open(file+".csv", 'w') as outfile:
		within_section = False
		for line in infile:
			if line.startswith('Class'):
				within_section = True
				continue  
			if line.startswith('Coverage'):
				within_section = False
				break
			if within_section:
		    		outfile.write(line)
try:
    file = sys.argv[1]
  
except:
	print("usage: python3 divsumSplitter.py file")

divsumSplitter_improved(file)
