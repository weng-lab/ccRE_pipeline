#!/usr/bin/env python3

import os
import sys
import argparse
import tempfile
import shutil
from decimal import Decimal

from helpers.utils import Utils, printt, numLines

def filterByDecimal(inFnp, outFnp, thresholdInt):
    with open(inFnp) as inF:
        with open(outFnp, 'w') as outF:
            for line in inF:
                fdr = str(line.rstrip().split("\t")[4])
                if "e" not in fdr and "E" not in fdr:
                    fdr = str('%.2E' % Decimal(fdr))
                if "+" in fdr:
                    outF.write(line)
                else:
                    order = int(fdr.split("-")[1])
                    if order > thresholdInt:
                        outF.write(line)

class CallDHSs(object):
    def __init__(self, args):
        self.args = args
        self.genome = args.genome

    def run(self):
        tmpDir = tempfile.mkdtemp()
        printt("tmpDir is", tmpDir)
        try:
            self._run(tmpDir)
            if not self.args.debug:
                shutil.rmtree(tmpDir)
            return 0
        except:
            raise

    def _run(self, tmpDir):
        enrichFnp = os.path.join("/data/projects/encode/data",
                                 self.args.DNaseExpAcc,
                                 self.args.DNaseBamAcc + ".bam.hotspots",
                                 self.args.DNaseBamAcc + ".allcalls.starch")

        bedFnp = os.path.join(tmpDir, self.args.DNaseBamAcc + ".bed")
        cmds = ["unstarch",
                enrichFnp,
                '>', bedFnp]
        printt("unstarch-ing", enrichFnp)
        Utils.runCmds(cmds)

        printt("Step 1 ...")
        inputFnp = bedFnp

        for i in range(2, self.args.minP + 1):
            cutoff = "1E-" + str(i)
            printt("starting", cutoff)
            outputFnp = os.path.join(tmpDir, self.args.DNaseBamAcc + '.' + cutoff + ".all.bed")
            filterByDecimal(inputFnp, outputFnp, i)

            mergeAndSizeFilteredFnp = os.path.join(tmpDir, self.args.DNaseBamAcc + '.' + cutoff + ".bed")
            cmds = ["bedtools",
                    "merge -d 1 -c 5 -o min",
                    "-i", outputFnp,
                    '|', """awk '{if ($3-$2 > 50) print $0}'""",
                    '>', mergeAndSizeFilteredFnp]
            Utils.runCmds(cmds)
            num = numLines(mergeAndSizeFilteredFnp)
            printt(cutoff, num)
            inputFnp = outputFnp

        printt("Step 2 ...")
        peaksFnp = os.path.join(tmpDir, "peaks")
        tmpFnp = os.path.join(tmpDir, "tmp")
        cmds = ["""awk '{if ($3-$2+1 < 350) print $0}'""",
                outputFnp,
                '>', peaksFnp]
        Utils.runCmds(cmds)

        for i in range(3, self.args.minP + 1):
            cutoff = "1E-" + str(i)
            printt("starting...", cutoff)
            cmds = ["bedtools",
                    "intersect",
                    "-v",
                    "-a", outputFnp,
                    "-b", peaksFnp,
                    '>',  tmpFnp]
            Utils.runCmds(cmds)

            cmds = [""" awk '{if ($3-$2+1 < 350) print $0}' """,
                    tmpFnp,
                    '>>', peaksFnp]
            Utils.runCmds(cmds)

        dhssFnp = os.path.join(tmpDir, self.args.DNaseBamAcc + ".DHSs.bed")
        shutil.move(peaksFnp, dhssFnp)

        excludedFnp = os.path.join(tmpDir, self.args.DNaseBamAcc + ".excluded.bed")
        cmds = ["bedtools",
                "intersect",
                "-v",
                "-a", os.path.join(tmpDir, self.args.DNaseBamAcc + '.' + "1E-" + str(self.args.minP) + ".bed"),
                "-b", dhssFnp,
                '>', excludedFnp]
        Utils.runCmds(cmds)

        outputDir = "/home/mjp/output/Processed-DHSs"
        Utils.mkdir_p(outputDir)
        shutil.move(excludedFnp, outputDir)
        shutil.move(dhssFnp, outputDir)

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action="store_true", default=False)
    parser.add_argument("--genome", type = str, help = "genome")
    parser.add_argument("--DNaseExpAcc", type = str, help = "DNase experiment accession")
    parser.add_argument("--DNaseBamAcc", type = str, help = "DNase BAM file accession")
    parser.add_argument("--minP", type = int, default=4942,
                        help = "min Hotspot2 p-value")
    return parser.parse_args()

def main():
    args = parseArgs()
    print("args:", args)
    return CallDHSs(args).run()


if __name__ == "__main__":
    sys.exit(main())
