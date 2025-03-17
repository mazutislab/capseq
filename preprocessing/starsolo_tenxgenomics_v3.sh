#!/bin/bash

# get genomic dirs to array
genomic=( $(find $PWD/raw/fastq/ -type f -name "*_S1_R2_001.fastq.gz") )

# get barcode dirs to array
barcode=( $(find $PWD/raw/fastq/ -type f -name "*_S1_R1_001.fastq.gz") )

# run starsolo
for i in "${!genomic[@]}"; do
    
    # get parent directory
    prefix=$(basename -s "_R2_001.fastq.gz" "${genomic[$i]}")
    
    echo 'processing: ' $(basename "${prefix}")
    echo 'genomic: ' ${genomic[$i]}
    echo 'barcode: ' ${barcode[$i]}

    STAR \
        --genomeDir /references/genome/GRCh38_Ensembl_107/index/ \
        --readFilesIn ${genomic[$i]} ${barcode[$i]} \
        --soloCBwhitelist 3M-february-2018.txt \
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
        --soloType CB_UMI_Simple \
        --soloBarcodeReadLength 0 \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIdedup Exact
    
    echo 'done'

done