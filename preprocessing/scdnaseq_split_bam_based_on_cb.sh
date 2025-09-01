#!/bin/bash

# species to process
species=(GRCh38 GRCm39)

for sp in "${species[@]}"; do

  # get BAMs
  bams=( $(find "$PWD" -type f -name "*.sortedByCoord.out.dedup_${sp}.bam") )

  # split bam
  for i in "${!bams[@]}"; do

      prefix="${bams[$i]%.*}"
      out_dir="$(dirname "$prefix")/perCB_${sp}"
      mkdir -p "$out_dir"

      echo "processing: $(basename "$prefix")"
      echo "bam: ${bams[$i]}"

      # get top CBs
      samtools view -bh -D CB:<(head -n 100 "$PWD/qced_cell_barcodes_${sp}.txt") "${bams[$i]}" \
        | samtools split -d CB -u /dev/null -f "${out_dir}/perCB_%!.bam" -

      echo 'done'

  done
done
