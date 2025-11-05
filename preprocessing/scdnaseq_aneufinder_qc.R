#libs
library(tidyverse)
library(AneuFinder)

# path to aneufinder results
h_cna <- "path/to/perCB_GRCh38_out/MODELS/method-HMM"
m_cna <- "path/to/perCB_GRCm39_out/MODELS/method-HMM"

##
## human
##

h_hmm_ls <- list.files(
  h_cna, full.names = TRUE
)

h_models <- loadFromFiles(h_hmm_ls)

h_models <- lapply(h_models, function(x){
  x$ID <- str_extract(x$ID, "[A-Z]{5}_[A-Z]{6}_[A-Z]{5}")
  return(x)
})

h_qc <- getQC(h_models) %>% 
  as_tibble(rownames = "sample") %>% 
  filter(
    total.read.count > 5e4
  ) %>% 
  left_join(
    ., lapply(h_models, function(x){
      x$segments %>% as_tibble() %>% 
        summarise(
          avg.copy.number = mean(copy.number)
        )
    }) %>% 
      bind_rows(.id = "sample")
  ) %>% 
  filter(
    !avg.copy.number %in% boxplot(avg.copy.number, plot = FALSE)$out
  )

h_models_qced <- h_models[h_qc$sample]

h_cl <- clusterByQuality(
  h_models_qced,
  measures=c('spikiness','num.segments','entropy','bhattacharyya','sos'),
  orderBy = "sos"
)

selected_files <- unlist(h_cl$classification[c(1:2)])

set.seed(50)
selected_files_50 <- sample(selected_files, 50)

heatmapGenomewide(
  h_models_qced[selected_files_50]
)

##
## mouse
##

m_hmm_ls <- list.files(
  m_cna, full.names = TRUE
)

m_models <- loadFromFiles(m_hmm_ls)

m_models <- lapply(m_models, function(x){
  x$ID <- str_extract(x$ID, "[A-Z]{5}_[A-Z]{6}_[A-Z]{5}")
  return(x)
})

m_qc <- getQC(m_models) %>% 
  as_tibble(rownames = "sample") %>% 
  filter(
    total.read.count > 5e4
  ) %>% 
  left_join(
    ., lapply(m_models, function(x){
      x$segments %>% as_tibble() %>% 
        summarise(
          avg.copy.number = mean(copy.number)
        )
    }) %>% 
      bind_rows(.id = "sample")
  ) %>% 
  filter(
    !avg.copy.number %in% boxplot(avg.copy.number, plot = FALSE)$out
  )

m_models_qced <- m_models[m_qc$sample]

m_cl <- clusterByQuality(
  m_models_qced,
  measures=c('spikiness','num.segments','entropy','bhattacharyya','sos'),
  orderBy = "sos"
)

m_selected_files <- unlist(m_cl$classification[c(1:5)])

set.seed(50)
m_selected_files_50 <- sample(m_selected_files, 50)


heatmapGenomewide(
  m_models_qced[m_selected_files_50]
)
