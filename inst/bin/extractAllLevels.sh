#!/bin/bash

#for reactome, need to sub out , and ;  .

awk '{gsub(/,/," ");print}' Ensembl2Reactome_All_Levels.txt | awk '{gsub(/;/," ");print}' | awk -F "\t" '{print $1"\t"$2"\t"$6"\t"$4}'  > reactomeAllLeveles.txt


#read into csv and print out
