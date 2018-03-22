

/app/src/0_Call-DHSs.py --genome=hg38 --DNaseExpAcc=ENCSR910OQF --DNaseBamAcc=ENCFF323MYG --debug --minP=3

docker run --rm -it -v /home/mjp/ccRE_pipeline/src:/app/src -v /home/mjp/ccRE_pipeline/info:/app/info -v /scratch/ccRE_pipeline/tmp:/tmp -v /data:/data -v /scratch/output/:/output purcaro/ccre_pipeline:001  /bin/bash
