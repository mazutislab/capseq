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
        --genomeDir /references/genome/GRCh38_GRCm39_Ensembl_107/index/ \
        --readFilesIn ${genomic[$i]} ${barcode[$i]} \
        --soloCBwhitelist $PWD/raw/barcodes/barcode_1836.txt \
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
        --soloBarcodeReadLength 0 \
        --soloCBposition 0_8_0_17 \
        --soloUMIposition 0_0_0_7 \
        --soloCBmatchWLtype EditDist_2 \
        --soloUMIdedup Exact

    echo 'done'

done