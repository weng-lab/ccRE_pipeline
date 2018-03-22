#!/usr/bin/env python3

import os
import sys
import argparse
import tempfile
import shutil

from helpers.utils import Utils, printt, numLines

#Jill E. Moore
#Weng Lab
#UMass Medical School
#ENCODE Encyclopedia Version 5

class CreateRDHSs(objects):
    def __init__(self, args):
        self.args = args

    def run(self):
        genome=hg38
        dir=~/Lab/ENCODE/Encyclopedia/V5/$genome-DNase/
        scriptDir=/home/moorej3/Projects/ENCODE/Encyclopedia/Version5/cRE-Pipeline

        cd $dir/Processed-DHSs

        echo -e "Combining DHSs..."
        cat output.* > tmp
        mv tmp $dir/$genome-DHS-All.bed

        echo -e "Filtering DHSs..."
        cd $dir
        #awk '{if ($9 > 3 && $3-$2 < 5000 && $7 > 0.075) print $0}' $genome-DHS-All.bed > $genome-DHS-FDR3-5K.bed
        #awk '{if ($9 > 3 ) print $0}' $genome-DHS-All.bed > $genome-DHS-FDR3-5K.bed
        awk '{if ($3-$2 > 150 && $5 > -0.75 ) print $0}' $genome-DHS-All.bed > $genome-DHS-Filtered.bed

        mkdir scratch
        cp $genome-DHS-Filtered.bed scratch/tmp.bed
        cd scratch

        echo -e "Sorting DHSs..."
        sort -k1,1 -k2,2n tmp.bed > sorted
        rm -f rPeaks
        num=$(wc -l sorted | awk '{print $1}')

        echo -e "Merging DHSs..."
        while [ $num -gt 0 ]
        do
            echo -e "\t" $num
            bedtools merge -i sorted -c 4,5 -o collapse,collapse > merge
            python $scriptDir/pick.best.peak.py merge > peak-list
            awk 'FNR==NR {x[$1];next} ($4 in x)' peak-list sorted >> rPeaks
            bedtools intersect -v -a sorted -b rPeaks > remaining
            mv remaining sorted
            num=$(wc -l sorted | awk '{print $1}')
        done

        mv rPeaks ../tmp.bed
        cd ../
        rm -r scratch

        echo -e "Accessioning rDHSs..."
        sort -k1,1 -k2,2n tmp.bed > sorted.bed
        python $scriptDir/make.cre.accession.py sorted.bed $genome rDHS \
            > $genome-rDHS-Filtered-Summary.txt
        awk '{print $1 "\t" $2 "\t" $3 "\t" $10}' $genome-rDHS-Filtered-Summary.txt \
            > $genome-rDHS-Filtered.bed

        rm sorted.bed
        rm tmp.bed

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action="store_true", default=False)
    parser.add_argument("--DNaseExpAcc", type = str, help = "DNase experiment accession")
    parser.add_argument("--DNaseBigWigAcc", type = str, help = "DNase BigWig file accession")
    return parser.parse_args()

def main():
    args = parseArgs()
    print("args:", args)
    return CreateRDHSs(args).run()


if __name__ == "__main__":
    sys.exit(main())
