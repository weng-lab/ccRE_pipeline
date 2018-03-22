#!/usr/bin/env python3

import os
import sys
import argparse
import tempfile
import shutil

from helpers.utils import Utils, printt, numLines

class ProcessDHSs(object):
    def __init__(self, args):
        self.args = args

    def run(self):
        dhsFnp = os.path.join("/output/Processed-DHSs", self.args.DNaseExpAcc + ".DHSs.bed")
        signalFnp = os.path.join("/data/projects/encode/data", self.args.DNaseExpAcc, self.args.DNaseBigWigAcc + ".bigWig")

        cp $dhs bed
        awk '{print $1 "\t" $2 "\t" $3 "\t" "'$dpeak'-"NR "\t" $4}' bed | sort -k4,4 > new
        awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' new > new.bed
        ~/bin/bigWigAverageOverBed $signal new.bed out.tab
        python $scriptDir/calculate.zscore.sh out.tab | sort -k1,1 > 1
        paste new 1 | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $7 "\t" "." "\t" $8 "\t" 1 "\t" $5}' > output.$j

        mkdir -p $dataDir/Processed-DHSs/
        mv output.$j $dataDir/Processed-DHSs/

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action="store_true", default=False)
    parser.add_argument("--DNaseExpAcc", type = str, help = "DNase experiment accession")
    parser.add_argument("--DNaseBigWigAcc", type = str, help = "DNase BigWig file accession")
    return parser.parse_args()

def main():
    args = parseArgs()
    print("args:", args)
    return ProcessDHSs(args).run()


if __name__ == "__main__":
    sys.exit(main())
