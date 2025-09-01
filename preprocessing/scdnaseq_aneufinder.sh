#!/bin/bash

# adjust the bed files
# hg38
awk 'BEGIN{OFS="\t"}
  $1 ~ /^chr([0-9]+|X|Y)$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {
    name = ($4=="" ? "." : $4); gsub(/[[:space:]]+/, "_", name);
    print $1,$2,$3,name,0,"."
  }' \
  "hg38-blacklist.v2.bed" \
  > "hg38-blacklist.v2.bed6"

# mm39
awk 'BEGIN{OFS="\t"}
  $1 ~ /^chr([0-9]+|X|Y)$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {
    name = $6; if (name=="") name="."; gsub(/[[:space:]]+/, "_", name);
    print $1,$2,$3,name,0,"."
  }' \
  "mm39.excluderanges.bed" \
  > "mm39.excluderanges.bed6"

# human
Rscript - <<RSCRIPT
suppressPackageStartupMessages({
  library(AneuFinder)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

Aneufinder(
  inputfolder = "perCB_GRCh38",
  outputfolder = "perCB_GRCh38_out",
  assembly = "hg38",
  correction.method = "GC",
  GC.BSgenome = BSgenome.Hsapiens.UCSC.hg38,
  binsizes = 1e5,
  blacklist = "hg38-blacklist.v2.bed6",
  method = c("HMM"),
  numCPU = 18
)
RSCRIPT

# mouse
Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(AneuFinder)
  library(BSgenome.Mmusculus.UCSC.mm39)
})

Aneufinder(
  inputfolder = "perCB_GRCm39",
  outputfolder = "perCB_GRCm39_out",
  assembly = "mm39",
  correction.method = "GC",
  GC.BSgenome = BSgenome.Mmusculus.UCSC.mm39,
  binsizes = 1e5,
  blacklist = "mm39.excluderanges.bed6",
  method = c("HMM"),
  numCPU = 18
)
RSCRIPT


# deactivate env
conda deactivate

# unload modules
module purge