

/app/src/0_Call-DHSs.py --genome=hg38 --DNaseExpAcc=ENCSR910OQF --DNaseBamAcc=ENCFF323MYG

docker run --rm -it -v /home/mjp/ccRE_pipeline/src:/app/src -v /home/mjp/ccRE_pipeline/info:/app/info -v /data:/data purcaro/ccre_pipeline:001  /bin/bash

