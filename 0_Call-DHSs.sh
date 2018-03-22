#!/bin/bash
#SBATCH --nodes 1
#SBATCH --time=12:00:00
#SBATCH --mem=1G
#SBATCH --array=1-80
#SBATCH --output=/home/moorej3/Job-Logs/jobid_%A_%a.output
#SBATCH --error=/home/moorej3/Job-Logs/jobid_%A_%a.error
#SBATCH --partition=12hours


genome=hg38
jid=$SLURM_ARRAY_TASK_ID

dataDir=~/Lab/ENCODE/Encyclopedia/V5/$genome-DNase/Version4
scriptDir=~/Projects/ENCODE/Encyclopedia/Version5/cRE-Pipeline
hotspots=$dataDir/list
bedtools=~/bin/bedtools2/bin/bedtools

minP=4942


mkdir -p /tmp/moorej3/$SLURM_JOBID-$jid
cd /tmp/moorej3/$SLURM_JOBID-$jid

dset=$(cat $hotspots  | awk '{if (NR == '$jid') print $1}')
bam=$(cat $hotspots  | awk '{if (NR == '$jid') print $2}')

#wget https://www.encodeproject.org/files/$bam/@@download/$bam.bed.gz
#gunzip $bam.bed.gz
#enrich=/data/projects/encode/data/$dset/$bam.bam.hotspots.double/$bam.allcalls.starch

#if [ ! -f $enrich ]; then
enrich=/data/projects/encode/data/$dset/$bam.bam.hotspots/$bam.allcalls.starch
#fi

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
