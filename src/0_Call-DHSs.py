#!/usr/bin/env python3

import os
import sys
import argparse
import tempfile

class CallDHSs(object):
    def __init__(self, args):
        self.args = args
        self.genome = args.genome
        self.minP = 4942

    def run(self):
        with tempfile.TemporaryDirectory() as tmpDir:
            try:
                self._run(tmpDir)
                return 0
            except:
                raise
                return 1

    def _run(self, tmpDir):
        enrichFnp = os.path.join("/data/projects/encode/data",
                                 self.args.DNaseExpAcc,
                                 self.args.DNaseBamAcc + ".bam.hotspots",
                                 self.args.DNaseBamAcc + ".allcalls.starch")

~/bin/bedops/unstarch $enrich > $bam.bed

echo "Step 1 ..." >> ~/Job-Logs/jobid_$SLURM_JOBID"_"$jid.error
cp $bam.bed 1
for j in `seq 2 1 $minP`
do
    cutoff=$(awk 'BEGIN{print "1E-'$j'"}')
    echo $cutoff >> ~/Job-Logs/jobid_$SLURM_JOBID"_"$jid.error
    python $scriptDir/filter.long.double.py 1 $cutoff > 2
    $bedtools merge -d 1 -c 5 -o min -i 2 | \
        awk '{if ($3-$2 > 50) print $0}' > $bam.$cutoff.bed
    mv 2 1
    num=$(wc -l $bam.$cutoff.bed | awk '{print $1}')
    echo $cutoff $num
done

echo "Step 2 ..." >> ~/Job-Logs/jobid_$SLURM_JOBID"_"$jid.error
cutoff=1E-2
awk '{if ($3-$2+1 < 350) print $0}' $bam.$cutoff.bed > peaks
for j in `seq 3 1 $minP`
do
    echo -e "\t" $j >> ~/Job-Logs/jobid_$SLURM_JOBID"_"$jid.error
    cutoff=$(awk 'BEGIN{print "1E-'$j'"}')
    $bedtools intersect -v -a $bam.$cutoff.bed -b peaks > tmp
    awk '{if ($3-$2+1 < 350) print $0}' tmp >> peaks
done
mv peaks $bam.DHSs.bed

$bedtools intersect -v -a $bam.1E-$minP.bed -b $bam.DHSs.bed > $bam.Excluded.bed

mkdir -p $dataDir/Processed-DHSs
mv $bam.Excluded.bed $bam.DHSs.bed $dataDir/Processed-DHSs/
#mv * $dataDir/Processed-DHSs/

rm -r /tmp/moorej3/$SLURM_JOBID-$jid

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", type = str, help = "genome")
    parser.add_argument("--DNaseExpAcc", type = str, help = "DNase experiment accession")
    parser.add_argument("--DNaseBamAcc", type = str, help = "DNase BAM file accession")
    parser.add_argument("--DNaseBigWigAcc", type = str, help = "DNase BigWig file accession")
    parser.add_argument("--inputRow", type = int, default=0, help = "row number in input file")
    return parser.parse_args()

def main():
    args = parseArgs()
    print("args:", args)
    return CallDHSs(args).run()


if __name__ == "__main__":
    sys.exit(main())
