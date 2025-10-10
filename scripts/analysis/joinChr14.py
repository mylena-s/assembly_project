from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

path = '/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Tuca/rna/rRNA/ribotin_output0/morphs.fa'
meancov = 15 


C14_I = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Tuca/Chr14/merging/Ctuca_14_I.fa"
c14_II=  "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Tuca/Chr14/merging/Ctuca_14_II.fa"

rnaClust = Seq("")
for record in SeqIO.parse(path, "fasta"):
    cov = int(record.id.split("coverage")[1]) // meancov
    record.seq = record.seq * cov
    rnaClust = rnaClust + record.seq


N = Seq("N")
Ns = N * 1000

partI = SeqIO.read(C14_I, "fasta")
partII = SeqIO.read(c14_II, "fasta")
partII.seq = partII.seq.reverse_complement() 

joined = partI.seq + Ns + rnaClust + Ns + partII.seq 
Chr14 = SeqRecord(joined, id="Ctuca_14", name="Ctuca_14", description = "joined")

with open("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Tuca/Chr14/merging/Chr14_joined.fasta", "w") as w:
    SeqIO.write(Chr14, w, "fasta")