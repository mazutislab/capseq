#!/bin/bash

# get bams
bams=( $(find "$PWD" -type f -name "*sortedByCoord.out.bam") )

# run picard
for i in "${!bams[@]}"; do
    
    # get file name
    prefix="${bams[$i]%.*}"
    
    echo 'processing: ' "$(basename "${prefix}")"
    echo 'bam: ' "${bams[$i]}"

    tmp="${prefix}.xc_tagged.bam"
    
    # reformat CB tag to XC
    samtools view -h "${bams[$i]}" \
    | awk 'BEGIN{FS=OFS="\t"} 
           /^@/{print; next}
           {
             added=0
             for(i=12;i<=NF;i++){
               if($i ~ /^CB:Z:/){
                 cb=substr($i,6)
                 if(cb != "-"){ gsub(/_/,"-",cb); $0=$0 OFS "XC:Z:" cb }
                 else { $0=$0 OFS "XC:Z:" }
                 added=1; break
               }
             }
             if(!added){ $0=$0 OFS "XC:Z:" }
             print
           }' \
    | samtools view -b -o "$tmp" -

    # run picard on the temp BAM
    picard MarkDuplicates \
        -I "$tmp" \
        -O "${prefix}.dedup.bam" \
        -M "${prefix}.mark_duplicates_metrics.txt" \
        -BARCODE_TAG XC \
        -REMOVE_DUPLICATES true \
        --TMP_DIR /tmp/

    # cleanup
    rm -f "$tmp" "${tmp}.bai"

    echo 'done'
done
