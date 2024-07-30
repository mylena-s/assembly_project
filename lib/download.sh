# download orthodb 
wget https://data.orthodb.org/download/odb11v0_levels.tab.gz
wget https://data.orthodb.org/download/odb11v0_all_og_fasta.tab.gz
wget https://data.orthodb.org/download/odb11v0_level2species.tab.gz
gunzip odb*gz

# extract fish proteins using https://github.com/tomasbruna/orthodb-clades script
python3 selectclade.py odb11v0_all_og_fasta.tab odb11v0_levels.tab odb11v0_level2species.tab Actinopterygii > Actinopterygii.fasta

# download pfam
mkdir Pfam & cd Pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

#download eggnog database
wget https://github.com/eggnogdb/eggnog-mapper/blob/master/download_eggnog_data.py
download_eggnog_data.py -y --data_dir .

