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
#install qualimap
cd ../software
wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip
unzip qualimap_v2.3.zip
cd ../bin & ln -s ../software/qualimap_v2.3/qualimap .
# install htslib
cd ../software
git clone https://github.com/samtools/htslib.git
cd htslib
.configure --prefix=$PWD
make
make install
cd ../../bin & ln -s ../software/htslib/bgzip .