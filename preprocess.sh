#!/bin/bash
basedir=${PWD}
retdir=${basedir}/scripts/retrieval
# Create the results folder
cd $1
mkdir Results
# Integrate datasets (6 samples of RNAseq, tiling from Nicolas 2012 and BacillusRegNet from Misirli)
cp ${retdir}/files/* ./
cp ${retdir}/rgife/configuration.conf Results/configuration.conf
Rscript ${retdir}/merge_and_normalize.R $2
rm tiling_array_expressions.txt -f
rm tiling_array_regulations.txt -f
rm gene_regulations.gb -f
