#!/bin/bash

# adjust the bed files - hg38
awk 'BEGIN{OFS="\t"}
  $1 ~ /^chr([0-9]+|X|Y)$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {
    name = ($4=="" ? "." : $4); gsub(/[[:space:]]+/, "_", name);
    print $1,$2,$3,name,0,"."
  }' \
  "hg38-blacklist.v2.bed" \
  > "hg38-blacklist.v2.bed6"

# adjust the bed files - mm39
awk 'BEGIN{OFS="\t"}
  $1 ~ /^chr([0-9]+|X|Y)$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {
    name = $6; if (name=="") name="."; gsub(/[[:space:]]+/, "_", name);
    print $1,$2,$3,name,0,"."
  }' \
  "mm39.excluderanges.bed" \
  > "mm39.excluderanges.bed6"

# human (K-562 cells - hypertriploid)
Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(AneuFinder)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(foreach)
  library(doParallel)
  library(rtracklayer)
  library(GenomicRanges)
})

calculate_weighted_gini <- function(x, w = NULL) {
  x <- as.numeric(x)
  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as.numeric(w)
  }
  keep <- !is.na(x) & !is.na(w) & (w > 0)
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) return(NA_real_)
  total_w <- sum(w)
  total_weighted_x <- sum(w * x)
  if (total_weighted_x == 0) return(0)
  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- w[ord]
  cum_w_frac <- c(0, cumsum(w_sorted) / total_w)
  cum_x_frac <- c(0, cumsum(w_sorted * x_sorted) / total_weighted_x)
  dx <- diff(cum_w_frac)
  y_avg <- (head(cum_x_frac, -1) + tail(cum_x_frac, -1)) / 2
  auc_lorenz <- sum(dx * y_avg)
  gini <- 1 - 2 * auc_lorenz
  gini <- max(0, min(1, as.numeric(gini)))
  return(gini)
}

calculate_coverage_metrics <- function(binned_data, blacklist_gr = NULL) {
  if (is(binned_data, "GRangesList")) {
    binned_data <- binned_data[[1]]
  }
  
  if (is(binned_data, "GRanges")) {
    if ("counts" %in% names(mcols(binned_data))) {
      counts <- mcols(binned_data)$counts
    } else {
      counts <- mcols(binned_data)$reads
    }
  }
  
  counts_for_gini <- counts
  n_excluded <- 0
  if (!is.null(blacklist_gr) && length(blacklist_gr) > 0) {
    overlaps <- GenomicRanges::findOverlaps(binned_data, blacklist_gr)
    blacklisted_bins <- unique(S4Vectors::queryHits(overlaps))
    
    if (length(blacklisted_bins) > 0) {
      counts_for_gini <- counts[-blacklisted_bins]
      n_excluded <- length(blacklisted_bins)
    }
  }
  
  gini <- calculate_weighted_gini(counts_for_gini, w = NULL)
  
  return(list(
    gini_coefficient = gini,
    n_bins_excluded_from_gini = n_excluded
  ))
}

calculate_cn_gini_from_model <- function(model) {
  if (is.null(model$segments)) return(NA_real_)
  cn <- model$segments$copy.number
  w <- as.numeric(width(model$segments))
  if (length(cn) == 0 || length(w) == 0) return(NA_real_)
  gini <- calculate_weighted_gini(cn, w = w)
  return(gini)
}

inputfolder <- "path/to/perCB_GRCh38"
outputfolder <- "path/to/perCB_GRCh38_out"
assembly <- "hg38"
GC.BSgenome <- BSgenome.Hsapiens.UCSC.hg38
binsizes <- 5e5
stepsizes <- binsizes
blacklist <- "hg38-blacklist.v2.bed6"
states <- c("0-somy", "1-somy", "2-somy", "3-somy", "4-somy", "5-somy", "6-somy", "7-somy", "8-somy")
most.frequent.state <- "3-somy"
chromosomes <- paste0('chr', c(1:22, 'X', 'Y'))
eps <- 0.001
max.time <- 60
max.iter <- 5000 
num.trials <- 1000
cluster.plots <- TRUE
numCPU <- 18
reuse.existing.files <- TRUE
remove.duplicate.reads <- TRUE
min.mapq <- 10

blacklist_gr <- rtracklayer::import(blacklist, format = "bed")

binpath.uncorrected <- file.path(outputfolder, 'binned')
binpath.corrected <- paste0(binpath.uncorrected, '-GC')
modelpath <- file.path(outputfolder, 'MODELS', 'method-HMM')
qcpath <- file.path(outputfolder, 'QC')

dir.create(outputfolder, recursive = TRUE, showWarnings = FALSE)
dir.create(binpath.uncorrected, recursive = TRUE, showWarnings = FALSE)
dir.create(binpath.corrected, recursive = TRUE, showWarnings = FALSE)
dir.create(modelpath, recursive = TRUE, showWarnings = FALSE)
dir.create(qcpath, recursive = TRUE, showWarnings = FALSE)

cl <- parallel::makeCluster(numCPU)
doParallel::registerDoParallel(cl)
on.exit(parallel::stopCluster(cl))

df.chroms <- GenomeInfoDb::getChromInfoFromUCSC(assembly)
files <- list.files(inputfolder, full.names = TRUE, pattern = '\\.bam$')

firstline <- Rsamtools::scanBamHeader(files[1])[[1]]$targets
if (any(grepl('^chr', names(firstline)))) {
  df <- df.chroms[,c('chrom','size')]
} else {
  df <- df.chroms[,c('chrom','size')]
  df$chrom <- sub('^chr', '', df$chrom)
}

chrom.lengths <- df[,2]
names(chrom.lengths) <- df[,1]
chrom.lengths <- chrom.lengths[!is.na(chrom.lengths) & !is.na(names(chrom.lengths))]
chrom.lengths.df <- data.frame(chromosome = names(chrom.lengths), length = chrom.lengths)

bins <- fixedWidthBins(chrom.lengths = chrom.lengths, chromosomes = chromosomes, 
                      binsizes = binsizes, stepsizes = stepsizes)

write.table(chrom.lengths.df, file = file.path(outputfolder, 'chrominfo.tsv'), 
           sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

if (length(files) == 1) {
  cluster.plots <- FALSE
}

patterns <- paste0('binsize_', format(binsizes, scientific = TRUE, trim = TRUE), 
                  '_stepsize_', format(stepsizes, scientific = TRUE, trim = TRUE), '_')

# binning
parallel.helper <- function(file) {
  existing.binfiles <- grep(basename(file), list.files(binpath.uncorrected), value = TRUE)
  existing.binsizes.stepsizes <- sub('_reads.per.bin.*', '', sub('.*_binsize', 'binsize', existing.binfiles))
  binsizes.stepsizes.todo <- setdiff(sub('_$', '', patterns), existing.binsizes.stepsizes)
  if (length(binsizes.stepsizes.todo) > 0 | !reuse.existing.files) {
    binReads(file = file, assembly = chrom.lengths.df, pairedEndReads = FALSE, 
             binsizes = NULL, variable.width.reference = NULL, 
             bins = bins[as.character(binsizes.stepsizes.todo)],
             chromosomes = chromosomes, remove.duplicate.reads = remove.duplicate.reads, 
             min.mapq = min.mapq, blacklist = blacklist, 
             outputfolder.binned = binpath.uncorrected, save.as.RData = TRUE)
  }
}

temp <- foreach(file = files, .packages = c("AneuFinder")) %dopar% {
  parallel.helper(file)
}

# calculate Gini on raw data
binfiles_raw <- list.files(binpath.uncorrected, pattern = 'RData$', full.names = TRUE)
gini_raw <- data.frame()

for (binfile in binfiles_raw) {
  sample_name <- gsub("\\.RData$", "", basename(binfile))
  binned_data <- get(load(binfile))
  metrics <- calculate_coverage_metrics(binned_data, blacklist_gr = blacklist_gr)
  
  gini_raw <- rbind(gini_raw, data.frame(
    sample = sample_name,
    stage = "raw",
    gini = metrics$gini_coefficient,
    n_bins_excluded = metrics$n_bins_excluded_from_gini
  ))
}

write.table(gini_raw, file = file.path(qcpath, 'gini_raw.tsv'), 
           sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

# GC correction
parallel.helper <- function(pattern) {
  binfiles <- list.files(binpath.uncorrected, pattern = 'RData$', full.names = TRUE)
  binfiles <- grep(gsub('\\+', '\\\\+', pattern), binfiles, value = TRUE)
  binfiles.corrected <- list.files(binpath.corrected, pattern = 'RData$', full.names = TRUE)
  binfiles.corrected <- grep(gsub('\\+', '\\\\+', pattern), binfiles.corrected, value = TRUE)
  binfiles.todo <- setdiff(basename(binfiles), basename(binfiles.corrected))
  if (length(binfiles.todo) > 0 | !reuse.existing.files) {
    binfiles.todo <- paste0(binpath.uncorrected, .Platform$file.sep, binfiles.todo)
    if (length(binfiles.todo) > 0) {
      binned.data.list <- suppressMessages(correctGC(binfiles.todo, GC.BSgenome, same.binsize = TRUE))
      for (i1 in 1:length(binned.data.list)) {
        binned.data <- binned.data.list[[i1]]
        savename <- file.path(binpath.corrected, basename(names(binned.data.list)[i1]))
        save(binned.data, file = savename)
      }
    }
  }
}

temp <- foreach(pattern = patterns, .packages = c("AneuFinder")) %dopar% {
  parallel.helper(pattern)
}

# calculate Gini on GC-corrected data
binfiles_corrected <- list.files(binpath.corrected, pattern = 'RData$', full.names = TRUE)
gini_corrected <- data.frame()

for (binfile in binfiles_corrected) {
  sample_name <- gsub("\\.RData$", "", basename(binfile))
  binned_data <- get(load(binfile))
  metrics <- calculate_coverage_metrics(binned_data, blacklist_gr = blacklist_gr)
  
  gini_corrected <- rbind(gini_corrected, data.frame(
    sample = sample_name,
    stage = "GC_corrected",
    gini = metrics$gini_coefficient,
    n_bins_excluded = metrics$n_bins_excluded_from_gini
  ))
}

write.table(gini_corrected, file = file.path(qcpath, 'gini_GC_corrected.tsv'), 
           sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

gini_all <- rbind(gini_raw, gini_corrected)
write.table(gini_all, file = file.path(qcpath, 'gini_all.tsv'), 
           sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

# HMM analysis
binfiles <- list.files(binpath.corrected, full.names = TRUE, pattern = '.RData$')
binfiles <- grep(paste(gsub('\\+', '\\\\+', patterns), collapse = '|'), binfiles, value = TRUE)

cn_metrics <- data.frame()

parallel.helper <- function(file) {
  savename <- file.path(modelpath, basename(file))
  if (!file.exists(savename) | !reuse.existing.files) {
    model <- findCNVs(file, method = 'HMM', eps = eps, max.time = max.time, 
                     max.iter = max.iter, num.trials = num.trials, 
                     states = states, most.frequent.state = most.frequent.state)
    model$breakpoints <- getBreakpoints(model, fragments = NULL, confint = NULL)
    model$cn_gini <- calculate_cn_gini_from_model(model)
    save(model, file = savename)
    return(model)
  } else {
    return(get(load(savename)))
  }
}

models <- foreach(file = binfiles, .packages = c("AneuFinder")) %dopar% {
  parallel.helper(file)
}

for (model in models) {
  if (!is.null(model) && !is.null(model$ID)) {
    cn_metrics <- rbind(cn_metrics, data.frame(
      sample = model$ID,
      cn_gini = model$cn_gini,
      num_segments = length(model$segments),
      num_breakpoints = if(!is.null(model$breakpoints)) length(model$breakpoints) else 0
    ))
  }
}

write.table(cn_metrics, file = file.path(qcpath, 'cn_gini_metrics.tsv'), 
           sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
RSCRIPT

# mouse (NIH3T3 cells - diploid)
Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(AneuFinder)
  library(BSgenome.Mmusculus.UCSC.mm39)
  library(foreach)
  library(doParallel)
  library(rtracklayer)
  library(GenomicRanges)
})

calculate_weighted_gini <- function(x, w = NULL) {
  x <- as.numeric(x)
  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as.numeric(w)
  }
  keep <- !is.na(x) & !is.na(w) & (w > 0)
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) return(NA_real_)
  total_w <- sum(w)
  total_weighted_x <- sum(w * x)
  if (total_weighted_x == 0) return(0)
  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- w[ord]
  cum_w_frac <- c(0, cumsum(w_sorted) / total_w)
  cum_x_frac <- c(0, cumsum(w_sorted * x_sorted) / total_weighted_x)
  dx <- diff(cum_w_frac)
  y_avg <- (head(cum_x_frac, -1) + tail(cum_x_frac, -1)) / 2
  auc_lorenz <- sum(dx * y_avg)
  gini <- 1 - 2 * auc_lorenz
  gini <- max(0, min(1, as.numeric(gini)))
  return(gini)
}

calculate_coverage_metrics <- function(binned_data, blacklist_gr = NULL) {
  if (is(binned_data, "GRangesList")) {
    binned_data <- binned_data[[1]]
  }
  
  if (is(binned_data, "GRanges")) {
    if ("counts" %in% names(mcols(binned_data))) {
      counts <- mcols(binned_data)$counts
    } else {
      counts <- mcols(binned_data)$reads
    }
  }
  
  counts_for_gini <- counts
  n_excluded <- 0
  if (!is.null(blacklist_gr) && length(blacklist_gr) > 0) {
    overlaps <- GenomicRanges::findOverlaps(binned_data, blacklist_gr)
    blacklisted_bins <- unique(S4Vectors::queryHits(overlaps))
    
    if (length(blacklisted_bins) > 0) {
      counts_for_gini <- counts[-blacklisted_bins]
      n_excluded <- length(blacklisted_bins)
    }
  }
  
  gini <- calculate_weighted_gini(counts_for_gini, w = NULL)
  
  return(list(
    gini_coefficient = gini,
    n_bins_excluded_from_gini = n_excluded
  ))
}

calculate_cn_gini_from_model <- function(model) {
  if (is.null(model$segments)) return(NA_real_)
  cn <- model$segments$copy.number
  w <- as.numeric(width(model$segments))
  if (length(cn) == 0 || length(w) == 0) return(NA_real_)
  gini <- calculate_weighted_gini(cn, w = w)
  return(gini)
}

inputfolder <- "path/to/perCB_GRCm39"
outputfolder <- "path/to/perCB_GRCm39_out"
assembly <- "mm39"
GC.BSgenome <- BSgenome.Mmusculus.UCSC.mm39
binsizes <- 5e5
stepsizes <- binsizes
blacklist <- "mm39.excluderanges.bed6"
states <- c("0-somy", "1-somy", "2-somy", "3-somy", "4-somy", "5-somy")
most.frequent.state <- "2-somy"
chromosomes <- paste0('chr', c(1:19, 'X', 'Y'))
eps <- 0.001
max.time <- 60
max.iter <- 5000
num.trials <- 1000
cluster.plots <- TRUE
numCPU <- 18
reuse.existing.files <- TRUE
remove.duplicate.reads <- TRUE
min.mapq <- 10

blacklist_gr <- rtracklayer::import(blacklist, format = "bed")

binpath.uncorrected <- file.path(outputfolder, 'binned')
binpath.corrected <- paste0(binpath.uncorrected, '-GC')
modelpath <- file.path(outputfolder, 'MODELS', 'method-HMM')
qcpath <- file.path(outputfolder, 'QC')

dir.create(outputfolder, recursive = TRUE, showWarnings = FALSE)
dir.create(binpath.uncorrected, recursive = TRUE, showWarnings = FALSE)
dir.create(binpath.corrected, recursive = TRUE, showWarnings = FALSE)
dir.create(modelpath, recursive = TRUE, showWarnings = FALSE)
dir.create(qcpath, recursive = TRUE, showWarnings = FALSE)

cl <- parallel::makeCluster(numCPU)
doParallel::registerDoParallel(cl)
on.exit(parallel::stopCluster(cl))

df.chroms <- GenomeInfoDb::getChromInfoFromUCSC(assembly)
files <- list.files(inputfolder, full.names = TRUE, pattern = '\\.bam$')

firstline <- Rsamtools::scanBamHeader(files[1])[[1]]$targets
if (any(grepl('^chr', names(firstline)))) {
  df <- df.chroms[,c('chrom','size')]
} else {
  df <- df.chroms[,c('chrom','size')]
  df$chrom <- sub('^chr', '', df$chrom)
}

chrom.lengths <- df[,2]
names(chrom.lengths) <- df[,1]
chrom.lengths <- chrom.lengths[!is.na(chrom.lengths) & !is.na(names(chrom.lengths))]
chrom.lengths.df <- data.frame(chromosome = names(chrom.lengths), length = chrom.lengths)

bins <- fixedWidthBins(chrom.lengths = chrom.lengths, chromosomes = chromosomes, 
                      binsizes = binsizes, stepsizes = stepsizes)

write.table(chrom.lengths.df, file = file.path(outputfolder, 'chrominfo.tsv'), 
           sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

if (length(files) == 1) {
  cluster.plots <- FALSE
}

patterns <- paste0('binsize_', format(binsizes, scientific = TRUE, trim = TRUE), 
                  '_stepsize_', format(stepsizes, scientific = TRUE, trim = TRUE), '_')

# binning
parallel.helper <- function(file) {
  existing.binfiles <- grep(basename(file), list.files(binpath.uncorrected), value = TRUE)
  existing.binsizes.stepsizes <- sub('_reads.per.bin.*', '', sub('.*_binsize', 'binsize', existing.binfiles))
  binsizes.stepsizes.todo <- setdiff(sub('_$', '', patterns), existing.binsizes.stepsizes)
  if (length(binsizes.stepsizes.todo) > 0 | !reuse.existing.files) {
    binReads(file = file, assembly = chrom.lengths.df, pairedEndReads = FALSE, 
             binsizes = NULL, variable.width.reference = NULL, 
             bins = bins[as.character(binsizes.stepsizes.todo)],
             chromosomes = chromosomes, remove.duplicate.reads = remove.duplicate.reads, 
             min.mapq = min.mapq, blacklist = blacklist, 
             outputfolder.binned = binpath.uncorrected, save.as.RData = TRUE)
  }
}

temp <- foreach(file = files, .packages = c("AneuFinder")) %dopar% {
  parallel.helper(file)
}

# calculate Gini on raw data
binfiles_raw <- list.files(binpath.uncorrected, pattern = 'RData$', full.names = TRUE)
gini_raw <- data.frame()

for (binfile in binfiles_raw) {
  sample_name <- gsub("\\.RData$", "", basename(binfile))
  binned_data <- get(load(binfile))
  metrics <- calculate_coverage_metrics(binned_data, blacklist_gr = blacklist_gr)
  
  gini_raw <- rbind(gini_raw, data.frame(
    sample = sample_name,
    stage = "raw",
    gini = metrics$gini_coefficient,
    n_bins_excluded = metrics$n_bins_excluded_from_gini
  ))
}

write.table(gini_raw, file = file.path(qcpath, 'gini_raw.tsv'), 
           sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

# GC correction
parallel.helper <- function(pattern) {
  binfiles <- list.files(binpath.uncorrected, pattern = 'RData$', full.names = TRUE)
  binfiles <- grep(gsub('\\+', '\\\\+', pattern), binfiles, value = TRUE)
  binfiles.corrected <- list.files(binpath.corrected, pattern = 'RData$', full.names = TRUE)
  binfiles.corrected <- grep(gsub('\\+', '\\\\+', pattern), binfiles.corrected, value = TRUE)
  binfiles.todo <- setdiff(basename(binfiles), basename(binfiles.corrected))
  if (length(binfiles.todo) > 0 | !reuse.existing.files) {
    binfiles.todo <- paste0(binpath.uncorrected, .Platform$file.sep, binfiles.todo)
    if (length(binfiles.todo) > 0) {
      binned.data.list <- suppressMessages(correctGC(binfiles.todo, GC.BSgenome, same.binsize = TRUE))
      for (i1 in 1:length(binned.data.list)) {
        binned.data <- binned.data.list[[i1]]
        savename <- file.path(binpath.corrected, basename(names(binned.data.list)[i1]))
        save(binned.data, file = savename)
      }
    }
  }
}

temp <- foreach(pattern = patterns, .packages = c("AneuFinder")) %dopar% {
  parallel.helper(pattern)
}

# calculate Gini on GC-corrected data
binfiles_corrected <- list.files(binpath.corrected, pattern = 'RData$', full.names = TRUE)
gini_corrected <- data.frame()

for (binfile in binfiles_corrected) {
  sample_name <- gsub("\\.RData$", "", basename(binfile))
  binned_data <- get(load(binfile))
  metrics <- calculate_coverage_metrics(binned_data, blacklist_gr = blacklist_gr)
  
  gini_corrected <- rbind(gini_corrected, data.frame(
    sample = sample_name,
    stage = "GC_corrected",
    gini = metrics$gini_coefficient,
    n_bins_excluded = metrics$n_bins_excluded_from_gini
  ))
}

write.table(gini_corrected, file = file.path(qcpath, 'gini_GC_corrected.tsv'), 
           sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

gini_all <- rbind(gini_raw, gini_corrected)
write.table(gini_all, file = file.path(qcpath, 'gini_all.tsv'), 
           sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

# HMM analysis
binfiles <- list.files(binpath.corrected, full.names = TRUE, pattern = '.RData$')
binfiles <- grep(paste(gsub('\\+', '\\\\+', patterns), collapse = '|'), binfiles, value = TRUE)

cn_metrics <- data.frame()

parallel.helper <- function(file) {
  savename <- file.path(modelpath, basename(file))
  if (!file.exists(savename) | !reuse.existing.files) {
    model <- findCNVs(file, method = 'HMM', eps = eps, max.time = max.time, 
                     max.iter = max.iter, num.trials = num.trials, 
                     states = states, most.frequent.state = most.frequent.state)
    model$breakpoints <- getBreakpoints(model, fragments = NULL, confint = NULL)
    model$cn_gini <- calculate_cn_gini_from_model(model)
    save(model, file = savename)
    return(model)
  } else {
    return(get(load(savename)))
  }
}

models <- foreach(file = binfiles, .packages = c("AneuFinder")) %dopar% {
  parallel.helper(file)
}

for (model in models) {
  if (!is.null(model) && !is.null(model$ID)) {
    cn_metrics <- rbind(cn_metrics, data.frame(
      sample = model$ID,
      cn_gini = model$cn_gini,
      num_segments = length(model$segments),
      num_breakpoints = if(!is.null(model$breakpoints)) length(model$breakpoints) else 0
    ))
  }
}

write.table(cn_metrics, file = file.path(qcpath, 'cn_gini_metrics.tsv'), 
           sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
RSCRIPT
