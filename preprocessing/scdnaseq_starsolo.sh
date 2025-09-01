#!/bin/bash

# barcode lists  
bc1=( $(find $PWD -type f -name "bc1_list.txt") )
bc2=( $(find $PWD -type f -name "bc2_list.txt") )
bc3=( $(find $PWD -type f -name "bc3_list.txt") )

# get barcode files
barcode=( $(find $PWD/ -type f -name "*_R2_001.fastq.gz") )

# get genomic files
genomic=( $(find $PWD/ -type f -name "*_R1_001.fastq.gz") )

# run starsolo
for i in "${!genomic[@]}"; do
    
    # get file name
    prefix="${genomic[$i]%%.*}"
    
    echo 'processing: ' $(basename "${prefix}")
    echo 'genomic: ' ${genomic[$i]}
    echo 'barcode: ' ${barcode[$i]}

    STAR \
        --genomeDir 'references/genome/GRCh38_GRCm39_Ensembl_107_rRNA/index' \
        --readFilesIn ${genomic[$i]} ${barcode[$i]} \
        --soloCBwhitelist $bc1 $bc2 $bc3 \
        --runThreadN 16 \
        --outFileNamePrefix ${prefix} \
        --runDirPerm All_RWX \
        --readFilesCommand zcat \
        --outSAMunmapped Within \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CR CB \
        --alignIntronMax 1 \
        --soloType CB_UMI_Complex \
        --soloCBmatchWLtype 1MM \
        --soloCBposition 0_0_0_4 0_7_0_12 0_15_0_19 \
        --soloUMIposition 0_0_0_0
    
    echo 'done'

done
