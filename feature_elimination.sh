#!/bin/bash
basedir=${PWD}
retdir=${basedir}/scripts/retrieval

cd $1
# Convert to .arff
Rscript ${retdir}/To_arff.R

# change working directory
cd Results
mkdir Experiment
# Run the RGIFE several times
n=$2
i=1
while [ $i -le $n ]
do
    python ${retdir}/rgife/rgife.py configuration.conf data.arff
    mkdir Experiment/run$i
    mv iterations.tar.gz Experiment/run$i/iterations.tar.gz
    mv summary.txt Experiment/run$i/summary.txt
    mv BestIteration Experiment/run$i/BestIteration
    mv ReferenceIteration Experiment/run$i/ReferenceIteration
    Rscript ${retdir}/plot_regulation.R ${basedir}/${1} Experiment/run$i/summary.txt
    i=`expr $i + 1`
done

# Look for the min max and union biomarker set
python ${retdir}/rgife/polices.py Experiment $n

# Tidy
mkdir polices_result
mv Experimentmax_model.txt polices_result/Experimentmax_model.txt
mv Experimentmin_model.txt polices_result/Experimentmin_model.txt
mv Experimentunion_model.txt polices_result/Experimentunion_model.txt

cd ..
# Analytic/brute-force optimisation of the classifiers
python ${retdir}/feature_elimination.py SVM
python ${retdir}/feature_elimination.py RF

#Generate the gene regulatory table to import into cytoscape
#Rscript ${retdir}/Cytoscape/RGIFE_to_cytoscape.R
