# install inspector
git clone https://github.com/Maggi-Chen/Inspector.git
cd ../bin & ln -s ../software/Inspector/*py .
cd ../software

# install genomescope
git clone https://github.com/tbenavi1/genomescope2.0.git
cd genomescope2.0/
Rscript install.R
# install hifiasm 
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make
cd ../../bin & ../software/hifiasm/hifiasm
# install sambamba
wget https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz
gunzip sambamba-1.0.1-linux-amd64-static.gz 
chmod +x sambamba-1.0.1-linux-amd64-static
cd ../bin & ln -s ../software/sambamba-1.0.1-linux-amf64-static sambamba
# install freebayes
cd ../software
wget https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz 
gunzip freebayes-1.3.6-linux-amd64-static.gz
chmod +x freebayes-1.3.6-linux-amd64-static
cd ../bin/ & ln -s ../software/freebayes-1.3.6-linux-amd64-static freebayes
# install htslib
cd ../software
git clone https://github.com/samtools/htslib.git
cd htslib
.configure --prefix=$PWD
make
make install
cd ../../bin & ln -s ../software/htslib/bgzip .

# install vt
git clone https://github.com/atks/vt.git
#change directory to vt
cd vt
#update submodules
git submodule update --init --recursive
#run make, note that compilers need to support the c++0x standard
make
#you can test the build
make test
cd ../../bin 
ln -s ../software/vt/vt .
cd ../software

# install purge_dups #### not used
git clone https://github.com/dfguan/purge_dups.git
cd purge_dups/src && make
cd ../../bin
ln -s ../software/purge_dups/bin/* .
cd ../software

# download simplifyheaders from braker pipeline
cd ../bin/ & wget http://bioinf.uni-greifswald.de/bioinf/downloads/simplifyFastaHeaders.pl
chmod +x simplifyFastaHeaders.pl 
cd ../software
=======
# mitohifi
singularity pull mitohifi.sif docker://ghcr.io/marcelauliano/mitohifi:master

# satellite finder
git clone https://github.com/lh3/srf
cd srf && make
cd ../../bin
ln -s ../software/srf/srf .
ln -s ../software/srf/srf .
cd ../software

# k8 for srf
wget https://github.com/lh3/dipcall/releases/download/v0.3/dipcall-0.3_x64-linux.tar.bz2
tar -jxf dipcall-0.3_x64-linux.tar.bz2
cd ../bin
ln -s ../software/dipcall.kit/k8 .
cd ../software

# nessie for centromere analysis
git clone https://github.com/mylena-s/nessie.git
cd nessie & make
chmod +x to_wig.py
cd ../../bin 
ln -s ../software/nessie/nessie . 
ln -s ../software/nessie/to_wig.py .
cd ../software

#variantqc
singularity pull discvrseq.sif docker://ghcr.io/bimberlab/discvrseq:latest

# wig to bigwig and bigwig to bedgraph
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig
chmod +x wigToBigWig
cd ../bin/ & ln -s ../software/wigToBigWig . 
cd ../software
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigToBedGraph
chmod +x bigWigToBedGraph 
cd ../bin/
ln -s ../software/bigWigToBedGraph .
cd ../software

# plink2
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_amd_avx2_20240418.zip
unzip plink2_linux_amd_avx2_20240418.zip
cd ../bin
ln -s ../software/plink2 .
cd ../software
# vcftools
git clone https://github.com/vcftools/vcftools.git
cd vcftools
./autogen.sh
./configure --prefix=$PWD
make
make install
cd ../../bin
ln -s ../software/vcftools/bin/* .
cd ../software

# beagle
wget https://faculty.washington.edu/browning/beagle/beagle.01Mar24.d36.jar
chmod +x beagle*.jar
cd ../bin/
ln -s ../software/beagle.*.jar .
cd ../software
# twiss
git clone https://github.com/simonhmartin/genomics_general.git
git clone https://github.com/simonhmartin/twisst.git
cd ../bin & ln -s ../software/twiss/*.R .
ln -s ../software/twiss/*py .
ln -s ../software/genomics_general/*py .
ln -s ../software/genomics_general/phylo/*.py
cd ../software/genomics_general/phylo & cp ../genomics.py .
cd ../..
#
# install TE curation software 
git clone https://github.com/clemgoub/TE-Aid.git
cd TE-Aid
chmod +x blastndotplot.R consensus2genome.R Run-c2g.R reduce.cpp
cd ../../bin
ln -s ../software/TE-Aid/*.sh .
ln -s ../software/TE-Aid/*.R .
ln -s ../software/TE-Aid/*.cpp .
cd ../software
#INSTALL LONGPHASE
wget https://github.com/twolinin/longphase/releases/download/v1.7.2/longphase_linux-x64.tar.xz
tar -xJf longphase_linux-x64.tar.xz
cd ../bin
ln -s ../software/longphase_linux-x64 .
cd ../software
#download TE manual annotation
git clone https://github.com/annaprotasio/TE_ManAnnot.git
cd ../bin
ln -s ../software/TE_ManAnnot/bin/* .
	cd ../software
# install bwa
git clone https://github.com/lh3/bwa
cd bwa && make
cd ../../bin/
ln -s ../software/bwa/bwa
cd ../software/
# admixture
wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar -xzvf admixture_linux-1.3.0.tar.gz
cd ../bin/
ln -s ../software/dist/admixture_linux-1.3.0/admixture
# chromopainter
cd assembly_project/software/
wget https://people.maths.bris.ac.uk/~madjl/finestructure/plink2chromopainter.pl.zip
unzip plink2chromopainter.pl.zip
wget https://people.maths.bris.ac.uk/~madjl/finestructure/makeuniformrecfile.pl.zip
unzip makeuniformrecfile.pl.zip 
wget https://people.maths.bris.ac.uk/~madjl/finestructure-old/chromopainter-0.0.4.tar.gz
tar -xzvf chromopainter-0.0.4.tar.gz 
cd chromopainter-0.0.4/
./configure --prefix=$PWD
make
make install
cd ../bin
ln -s ../software/chromopainter-0.0.4/chromopainter
ln -s ../software/makeuniformrecfile.pl
ln -s ../software/plink2chromopainter.pl

#mvftools
cd ../software
git clone https://www.github.com/jbpease/mvftools



#deepvariant
wget https://github.com/dnanexus-rnd/GLnexus/releases/download/v1.4.1/glnexus_cli
chmod +x glnexus_cli
ln -s glnexus_cli ../bin/
