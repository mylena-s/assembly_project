# isolate query proteins and removes prot with no projections.
cat prot.fasta | grep "PROT | QUERY" -w -A 1 | grep "^\-\-$" -v | awk '{if ($1 ~ /^>/) printf $1"\t"; else print $0}' | sed 's/-//g' | awk -F "\t" '{if ($2 != "") print $1"\n"$2}' > query_prot.fasta

# separate Intact, probably intact and uncertain loss in a new bed file
grep PROJECTION loss_summ_data.tsv | awk '{if ($3 == "I" || $3 == "PI" || $3 == "UL") print $2}' | sort > I_PI_UL.txt
sort -k4,4 query_annotation.bed -o query_annotation.sorted.bed
join -1 1 -2 4 I_PI_UL.txt query_annotation.sorted.bed -t $'\t'  > query_annotation.I_PI_UL.temp
awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' query_annotation.I_PI_UL.temp > query_annotation.I_PI_UL.bed

# list query isoforms
awk '{print $3}' query_annotation.I_PI_UL.bed > query_annotation.I_PI_UL.list
sort query_annotation.I_PI_UL.list | uniq | wc -l query_annotation.I_PI_UL.uniq.list

# separate query isoforms proteins

seqtk subseq query_prot.fasta query_annotation.I_PI_UL.uniq.list > query_annotation.I_PI_UL.fasta
# as not all query annotation isosofrms of bed file are on prot file:
grep ">" query_annotation.I_PI_UL.fasta | sed 's/>//g' > query_annotation.I_PI_UL.list

#get isoforms-gene file for query prot
grep -Ff query_annotation.I_PI_UL.list query_isoforms_I_PI_UL.tsv > query_isoforms_I_PI_UL_uniq.tsv 
# modify with python to run omark

# bed to gtf
bed2gtf --bed query_annotation.I_PI_UL.bed --isoforms query_isoforms_I_PI_UL.tsv -o query_annotation.I_PI_UL.gtf
# standardize gff
agat_convert_sp_gxf2gxf.pl -g query_annotation.I_PI_UL.gtf -o query_annotation.I_PI_UL.gff

