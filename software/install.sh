# install inspector
git clone https://github.com/Maggi-Chen/Inspector.git
cd ../bin & ln -s ../software/Inspector/*py .
cd ../software

# install genomescope
git clone https://github.com/tbenavi1/genomescope2.0.git
cd genomescope2.0/
Rscript install.R
