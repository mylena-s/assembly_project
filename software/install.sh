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

# mitohifi
singularity pull mitohifi.sif docker://ghcr.io/marcelauliano/mitohifi:master

