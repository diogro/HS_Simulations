library(tidyverse)
library(data.table)

snp2Num = function(x){
  num = substr(x, start = 1, stop = 3)
  sapply(num, snpConvert)
}
snpConvert = function(x){
  if(x == '1/1')
    return(2)
  else if(x == '0/0')
    return(0)
  else if(x == '1/0' | x=='0/1')
    return(1)
  else 
    return(NA_integer_)
}
na.prop = function(x) sum(is.na(x))/length(x)

load("data/200615_binGlmScan_contTimeParam1234_slopePerTreat.RData")
load("data/200811_HSselectedSNPs_usedInCorEtc.RData")
glm_scan = binGlmScan_contTime_param1234_200615 %>%
  rename(interaction_p = "generation:treatment_p", 
         interaction = "generation:treatment") %>%
  filter(CHROM != "211000022280564")
freq_tab = freqs.merge_union_190426.hsSel
rm(binGlmScan_contTime_param1234_200615, freqs.merge_union_190426.hsSel)

top_snps = glm_scan %>% as_tibble() %>%
  filter(interaction_p < 1e-9 & emTrend.HS_p < 1e-4 ) %>%
  select(CHROM, POS, interaction_p, emTrend.HS, emTrend.HS_p) %>%
  filter(CHROM != "211000022280564")
#write_delim(top_snps[c("CHROM", "POS")], file = "top_snps_positions.txt", delim = "\t")
#Run filtering in server...
topVCF = fread("data/topSnpsSmall.recode.vcf", skip = "#CHROM") %>% 
  mutate_at(vars(matches("^N")), snp2Num) %>%
  rename(CHROM = '#CHROM') %>%
  filter(CHROM != "211000022280564")  

topVCF_list = lapply(1:6, function(pop) topVCF %>%
  select(POS, CHROM, matches(paste0("^N[", pop, "]"))))

x = lapply(topVCF_list, function (x) apply(x[,-c(1, 2)], 1, na.prop)) %>% as.data.frame()
names(x) = paste0("N", 1:6)
if(!require(GGally)){install.packages("GGally"); library(GGally)}
ggpairs(x)

freq_tab = freq_tab %>%
  filter(POS %in% topVCF$POS) %>%
  filter(CHROM != "211000022280564")
setkey(freq_tab, CHROM, POS)
glm_scan = glm_scan  %>%
  arrange(desc(interaction_p))%>%
  filter(POS %in% topVCF$POS & interaction_p < 1e-20 & emTrend.HS_p < 1e-4)

topSnptbl = tibble(CHROM = character(), POS = integer(), interaction_p = numeric())
for(chr in unique(glm_scan$CHROM)){
  dt = glm_scan[CHROM==chr]
  for(j in 1:20){
    index = 1
    dt = dt[order(interaction_p)]
    while(TRUE){
      current_top = dt[index,.(CHROM, POS, interaction_p), key = POS]
      na_proportions = c(na.prop(topVCF[current_top, .SD, .SDcols = names(topVCF) %like% "N[4]" , on = "POS"]),
                         na.prop(topVCF[current_top, .SD, .SDcols = names(topVCF) %like% "N[5]" , on = "POS"]),
                         na.prop(topVCF[current_top, .SD, .SDcols = names(topVCF) %like% "N[6]" , on = "POS"]))
      print(c(index, mean(na_proportions)))
      index = index + 1
      if(mean(na_proportions) < 0.3) {print("found it"); break}
      if(mean(na_proportions) == 1) stop("No good deed.")
    }
    topSnptbl = rbind(topSnptbl, current_top)
    dt = dt[abs(current_top[,POS] - POS) > 1e6]
    if(nrow(dt) < 2) break
  }
}
setkey(topSnptbl, CHROM, POS)

ggplot(glm_scan, aes(POS, -log10(interaction_p))) + geom_point() + 
  facet_wrap(~CHROM, nrow = 1) + 
  geom_text_repel(data= topSnptbl, aes(label = POS))


topVCF = lapply(topVCF, function(x) x[topSnptbl])
topfreq_tab = freq_tab[topSnptbl]

target = t(topVCF[POS==top_snp, .SD, .SDcols = names(topVCF) %like% "^N[4]"])
partner = t(topVCF[6000, .SD, .SDcols = names(topVCF) %like% "^N[4]"])

for(i in 1:nrow(topSnptbl)){
  for(j in 1:nrow(topSnptbl)){
    if(i!=j){
      t_pos = topSnptbl[i, POS]
      target = t(topVCF[POS==t_pos, .SD, .SDcols = names(topVCF) %like% "^N[4]"])
      partner = t(topVCF[j, .SD, .SDcols = names(topVCF) %like% "^N[4]"])
      table(target, partner)
    }
  }
}
