#!/bin/bash

# get genomic dirs to array
genomic=( $(find $PWD/raw/fastq/ -type f -name "*_2.fastq.gz") )

# get barcode dirs to array
barcode=( $(find $PWD/raw/fastq/ -type f -name "*_1.fastq.gz") )

# run starsolo
for i in "${!genomic[@]}"; do
    
    # get parent directory
    prefix=$(basename -s "_2.fastq.gz" "${genomic[$i]}")
    
    echo 'processing: ' $(basename "${prefix}")
    echo 'genomic: ' ${genomic[$i]}
    echo 'barcode: ' ${barcode[$i]}

    STAR \
        --genomeDir /references/genome/GRCh38_Ensembl_107/index/ \
        --readFilesIn ${genomic[$i]} ${barcode[$i]} \
        --soloCBwhitelist cls1.txt cls2.txt cls3.txt \
        --runThreadN 16 \
        --limitBAMsortRAM 88357076687 \
        --outFileNamePrefix starsolo/${prefix}/ \
        --sjdbOverhang 100 \
        --outSAMunmapped Within \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingBinsN 20 \
        --outSAMattributes NH HI nM AS CR UR CB UB sS sQ sM GX GN \
        --runDirPerm All_RWX \
        --readFilesCommand zcat \
        --soloFeatures GeneFull \
        --soloMultiMappers EM \
        --soloType CB_UMI_Complex \
        --soloCBposition 0_0_0_8 0_21_0_29 0_43_0_51 \
        --soloUMIposition 0_52_0_59 \
        --soloCBmatchWLtype EditDist_2 \
        --soloUMIdedup Exact
    
    echo 'done'

done