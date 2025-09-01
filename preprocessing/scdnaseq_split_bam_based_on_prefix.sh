#!/bin/bash

bams=( $(find "$PWD" -type f -name "*.bam") )

# split bams by species
for bam in "${bams[@]}"; do
  prefix="${bam%.*}"
  base="$(basename "$prefix")"

  echo "processing: $base"
  echo "bam: $bam"

  # add sp:Z:<species>
  sp_tagged="${prefix}.sp_tagged.bam"
  samtools view -h "$bam" \
  | awk 'BEGIN{FS=OFS="\t"}
         /^@/{print; next}
         {
           sp="";
           if($3!="*" && $3!="="){ n=split($3,a,/_/); if(n>1) sp=a[1] }
           if(sp!="") $0=$0 OFS "sp:Z:" sp
           print
         }' \
  | samtools view -b -o "$sp_tagged" -

  # split by tag
  samtools split -d sp -f "${prefix}_%!.bam" -u "${prefix}_unassigned.bam" "$sp_tagged"

  # for each species, copy and index
  for sp in GRCh38 GRCm39; do
    out="${prefix}_${sp}.bam"
    [[ -s "$out" ]] || continue

    tmp_in="${SLURM_TMPDIR:-/tmp}/$(basename "$out")"
    tmp_fixed="${tmp_in}.fixed"

    cp -f "$out" "$tmp_in"

    # keep only @SQ with SN starting with "<sp>_", strip that prefix in SN,
    # and strip the same prefix in RNAME (col 3) and RNEXT (col 7).
    samtools view -h "$tmp_in" \
    | awk -v SP="${sp}_" 'BEGIN{FS=OFS="\t"}
        # header: keep only species contigs and drop the prefix in SN:
        /^@SQ/{
          sn="";
          for(i=1;i<=NF;i++) if($i ~ /^SN:/) sn=substr($i,4);
          if(sn !~ "^"SP) next;
          for(i=1;i<=NF;i++){
            if($i ~ /^SN:/){
              ref=substr($i,4);
              sub("^"SP,"",ref);
              $i="SN:" ref
            }
          }
          print; next
        }
        /^@/ { print; next }

        # strip prefix in RNAME and RNEXT if present
        {
          if($3!="*" && $3!="=") sub("^"SP,"",$3);
          if($7!="*" && $7!="=") sub("^"SP,"",$7);
          print
        }' \
    | samtools view -b -o "$tmp_fixed" -

    samtools index "$tmp_fixed"
    mv -f "$tmp_fixed" "$out"
    mv -f "$tmp_fixed.bai" "${out}.bai"
    rm -f "$tmp_in"
  done

  # cleanup
  rm -f "$sp_tagged" "${sp_tagged}.bai"

  echo 'done'
done
