if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}

load("./data/200811_HSselectedSNPs_usedInCorEtc.RData")
freqs = freqs.merge_union_190426.hsSel
freq_G1 = freqs %>%
  select(matches("G1_N[456]")) 
g1_AFS  = abs(rowMeans(freq_G1) - 0.5)
hist(g1_AFS, breaks = 100)

