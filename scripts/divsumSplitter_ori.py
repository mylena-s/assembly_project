import sys
import csv

def divsum_splitter(file):    
    file_modif = file.split('.')[0]
    with open(file) as file, open( file_modif + '.csv', 'w') as modif:    
        orig_table = csv.reader(file)
        for row in orig_table:
            for index in row:
                if 'Div' in index:
                    file_csv = csv.writer(modif)
                    file_csv.writerow(row)
                    for row in orig_table:
                        csv.writer(modif)
                        file_csv.writerow(row)
    return file_modif + '.tsv'
    
    
try:
    file = sys.argv[1]
  
except:
	print("usage: python3 divsumSplitter.py file")

divsum_splitter(file)
