#!/bin/bash
# $1-Base directory (may contain CountTables folder)
# $2-Number of RGIFE iterations
# $3-Correctly formatted locations of stress samples


bash preprocess.sh $1 $3

bash feature_elimination.sh $1 $2
