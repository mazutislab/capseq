#!/bin/bash

# load modules
module load nextflow
module load singularity

#run pipe
nextflow run jsimonas/solo-in-drops \
	--run_dir '/path/to/sequencing/run/directory/' \
	--outdir '/path/to/output/directory/' \
    --sample_sheet '/path/to/extended_sample_sheet.xlsx' \
	--star_index 'references/genome/GRCh38_GRCm39_Ensembl_107/index/' \
    --scrna_protocol 'splitpool' \
    --solo_features 'GeneFull' \
    --sequencer 'nextseq' \
	-profile singularity \
	-ansi-log false \
    -r dev

# unload modules
module purge