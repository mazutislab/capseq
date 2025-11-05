#!/bin/bash

# species to process
species=(GRCh38 GRCm39)

for sp in "${species[@]}"; do
    
    # get BAMs
    bams=( $(find "$PWD" -type f -name "*.sortedByCoord.out.dedup_${sp}.bam") )

    barcode_file="$PWD/qced_cell_barcodes_${sp}.txt"
    N_CB=$(wc -l < "$barcode_file") 
    
    # split bam
    for i in "${!bams[@]}"; do
    
        prefix="${bams[$i]%.*}"
        out_dir="$(dirname "$prefix")/perCB_${sp}"
        mkdir -p "$out_dir"

        echo "processing: $(basename "$prefix")"
        echo "bam: ${bams[$i]}"
    
        ulimit -n $(( N_CB > 1014 ? N_CB + 10 : 1024 ))
      
        # get qced CBs
        samtools view -bh -D CB:<(cat "$barcode_file") "${bams[$i]}" \
        | samtools split -d CB -M "$N_CB" -u /dev/null -f "${out_dir}/perCB_%!.bam" -

      echo 'done'

  done
done
