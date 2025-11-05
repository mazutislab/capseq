#libs
library(tidyverse)
library(AneuFinder)

## data
h_cna <- "path/to/perCB_GRCh38_out/MODELS/method-HMM"
m_cna <- "path/to/perCB_GRCm39_out/MODELS/method-HMM"

###
### human
###

h_hmm_ls <- list.files(
  h_cna, full.names = TRUE
)

h_models <- loadFromFiles(h_hmm_ls)

h_models <- lapply(h_models, function(x){
  x$ID <- str_extract(x$ID, "[A-Z]{5}_[A-Z]{6}_[A-Z]{5}")
  return(x)
})

# calculate coverage
h_cov_tb <- lapply(h_models, function(x){
  x$bincounts %>% as_tibble() %>% 
    mutate(
      origin = "GRCh38",
      chrom = seqnames,
      bin = row_number(),
      count = counts,
      coverage = counts/width,
    ) %>% 
    select(origin, chrom, bin, start, end, count, coverage)
  }) %>% 
  bind_rows(
    .id = "CB"
  ) %>% 
  mutate(
    CB = str_extract(CB, "[A-Z]{5}_[A-Z]{6}_[A-Z]{5}")
  )

###
### mouse
###

m_hmm_ls <- list.files(
  m_cna, full.names = TRUE
)

m_models <- loadFromFiles(m_hmm_ls)

m_models <- lapply(m_models, function(x){
  x$ID <- str_extract(x$ID, "[A-Z]{5}_[A-Z]{6}_[A-Z]{5}")
  return(x)
})

# calculate coverage
m_cov_tb <- lapply(m_models, function(x){
  x$bincounts %>% as_tibble() %>% 
    mutate(
      origin = "GRCm39",
      chrom = seqnames,
      bin = row_number(),
      count = counts,
      coverage = counts/width,
    ) %>% 
    select(origin, chrom, bin, start, end, count, coverage)
  }) %>% 
  bind_rows(
    .id = "CB"
  ) %>% 
  mutate(
    CB = str_extract(CB, "[A-Z]{5}_[A-Z]{6}_[A-Z]{5}")
  )

sc_cov <- bind_rows(
  h_cov_tb,
  m_cov_tb
)



