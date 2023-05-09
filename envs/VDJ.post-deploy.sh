#!/usr/bin/bash
log=$CONDA_PREFIX/deploy.log
srcdir=src
mkdir -p $srcdir

##################################
### baizer
##################################

# R -e "install.packages('baizer', repos='http://cran.r-project.org')" >> $log 2>&1
R -e "devtools::install_github('william-swl/baizer')" >> $log 2>&1
echo ">>> baizer installed" >> $log 2>&1

##################################
### genogamesh
##################################

R -e "devtools::install_github('william-swl/genogamesh')" >> $log 2>&1
echo ">>> genogamesh installed" >> $log 2>&1
