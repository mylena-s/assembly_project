conda activate variant_calling
samtools mpileup -q 30 -Q 20 -C 50 -u -r Chr1 -f CBFH00195_renamed.fasta CBFH00195_renamed_CBFH00195_trimmed.fq_sorted.RG.bam | bcftools call -c -V indels | bamCaller.py 30 out_mask.bed.gz | gzip -c > out.vcf.gz
export PATH=~/genome_report/software/msmc-tools/:$PATH   
generate_multihetsep.py --mask out_mask.bed.gz --chr Chr1 out.vcf.gz > chr1.multihetsep.txt
~/genome_report/software/msmc2_Linux -t 6 -p 1*2+16*1+1*2 -o out.msmc2 chr1.multihetsep.txt

