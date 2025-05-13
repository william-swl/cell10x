#!/usr/bin/bash
log=$CONDA_PREFIX/deploy.log

##################################
### plutor
##################################

R -e "devtools::install_github('william-swl/plutor')" >> $log 2>&1
echo ">>> plutor installed" >> $log 2>&1
