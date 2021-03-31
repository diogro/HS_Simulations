if(!require(plyr)){install.packages("plyr"); library(plyr)}
library(tidyverse)
library(data.table)
if(!require(patchwork)){install.packages("patchwork"); library(patchwork)}
if(!require(ggthemes)){install.packages("ggthemes"); library(ggthemes)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}


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
getP = function(x){
  A1 = 2*x[1,1] + x[2,1]
  B1 = 2*x[3,1] + x[2,1]
  
  A2 = 2*x[1,3] + x[2,3]
  B2 = 2*x[3,3] + x[2,3]
  
  dat = data.frame(A = c(A1, A2), B = c(B1, B2), partner = c(0, 1))
  m = glm(cbind(A, B) ~ partner, family = "binomial", data = dat)
  return(m)
}
p_change = function(s, p1){
  p_t = numeric(100)
  p_t[1] = p1
  for(t in 2:100){
    p = p_t[t-1]
    p_t[t] = p + (-((1/2) * (s*p*(1-p))) / (1-s*p))
  }
  p_t
}
LL = function(s, p1, p100){
  p_t = p_change(s, p1)
  abs(p_t[100] - p100)^2
}
estimate_s <- function(p100, p1) {
  possible_s = seq(-1, 1, by=0.0001)
  s_vec = sapply(possible_s, LL, p1, p100)
  s_hat = possible_s[which(min(s_vec) == s_vec)]
  #LL(s_hat, p1, p100)
  #plot(p_change(s_hat, p1), ylim = c(0, 1), pch = 19)
  #points(data.frame(c(1, 100), c(p1, p100)), col = "red", pch = 19)
  s_fit = optim(c(s = s_hat), function(s) LL(s, p1, p100), lower = -1, upper = 1, method = "Brent")$par
  s_fit
}
logit = function(p) log(p/(1-p))
invlogit = function(x) exp(x)/(1+exp(x))

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
topVCF = fread("data/topSnps.recode.vcf", skip = "#CHROM") %>% 
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
  for(j in 1:100){
    index = 1
    dt = dt[order(interaction_p)]
    while(TRUE){
      current_top = dt[index,.(CHROM, POS, interaction_p)]
      na_proportions = c(na.prop(topVCF[current_top, .SD, .SDcols = names(topVCF) %like% "N[4]" , on = "POS"]),
                         na.prop(topVCF[current_top, .SD, .SDcols = names(topVCF) %like% "N[5]" , on = "POS"]),
                         na.prop(topVCF[current_top, .SD, .SDcols = names(topVCF) %like% "N[6]" , on = "POS"]))
      index = index + 1
      if(mean(na_proportions) < 0.25) {add = TRUE; break}
      if(mean(na_proportions) == 1){ add = FALSE; break; } 
    }
    if(add){
      topSnptbl = rbind(topSnptbl, current_top)
    }
    dt = dt[abs(current_top[,POS] - POS) > 1e5]
    if(nrow(dt) < 2) break
  }
}
setkey(topSnptbl, CHROM, POS)

if(!require(ggrepel)){install.packages("ggrepel"); library(ggrepel)}

ggplot(glm_scan, aes(POS, -log10(interaction_p))) + geom_point() + 
  facet_wrap(~CHROM, nrow = 1) + 
  geom_text_repel(data= topSnptbl, aes(label = POS))

topVCF_list = lapply(topVCF_list, function(x) x[topSnptbl, on = "POS"])
topfreq_tab = freq_tab[topSnptbl]

getFreqTable = function(pos1, pos2, pop = '456'){
  target  = t(topVCF[POS==pos1, .SD, .SDcols = names(topVCF) %like% paste0("^N[", pop, "]")])
  partner = t(topVCF[POS==pos2, .SD, .SDcols = names(topVCF) %like% paste0("^N[", pop, "]")])
  table(target, partner)
}
getSTable = function(pop = '456'){
  sig_pairs = tibble(target_CHROM = character(), target_POS = integer(), 
                     partner_CHROM = character(), partner_POS = integer(), 
                     p_val = numeric(), chisq.p_val = numeric(), chisq = numeric(),
                     p0 = numeric(), p2 = numeric())
  sig_names = names(sig_pairs)
  for(i in 1:nrow(topSnptbl)){
    for(j in 1:nrow(topSnptbl)){
      if(i!=j){
        t_pos = topSnptbl[i, POS]
        p_pos = topSnptbl[j, POS] 
        tb = getFreqTable(t_pos, p_pos, pop = pop)
        if(!all(dim(tb) == c(3, 3))) break
        chi_sq = chisq.test(tb)
        p_model = getP(tb) # Allele frequency of the 0 allele in both partner homozygotes
        p_val = summary(p_model)$coefficients[,"Pr(>|z|)"]["partner"]
        if(!is.na(p_val))
          if(p_val<0.1){ 
            sig_pairs = rbind(sig_pairs, data.frame(topSnptbl[POS==t_pos, .(CHROM, POS)],
                                                    topSnptbl[POS==p_pos, .(CHROM, POS)], p_val,
                                                    chisq.p_val = chi_sq$p.value, 
                                                    chisq = chi_sq$statistic,
                                                    p0 = invlogit(p_model$coefficients[1]), 
                                                    p2 = invlogit(sum(p_model$coefficients))))
          }
      }
    }
  }
  names(sig_pairs) = sig_names
  sig_pairs
}
sig_pairs_HS = getSTable('456')

ggplot(sig_pairs_HS, aes(-log10(chisq.p_val), -log10(p_val))) + 
  geom_point() + geom_abline(intercept = 0, slope = 1) + geom_vline(xintercept = 4 )

x = sig_pairs_HS %>%
  filter(target_CHROM != partner_CHROM, chisq.p_val < 1e-4) %>%
  arrange(chisq.p_val)

pop = '456'
target = 2955361
partner = 1308031

getS <- function(target, partner, pop = '456') {
  (tb = getFreqTable(target, partner))
  m = getP(tb)
  p100 = invlogit(c(m$coefficients[1], m$coefficients[1] + m$coefficients[2]))
  p0 = mean(as.numeric(topfreq_tab[POS==target, .SD, .SDcols = names(topfreq_tab) %like% paste0("^G1_N[", pop, "]")]))
  c(s_0 = estimate_s(p100[1], p0),
    s_2 = estimate_s(p100[2], p0),
    p_0 = p0)
}
getS(target, partner)

x$ID = 1:nrow(x)
s_df_estimates = ddply(x, .(ID), 
                       function(x) getS(x$target_POS, x$partner_POS, '456'), 
                       .progress = "text")
y = inner_join(x, s_df_estimates)

stat = ggplot(y, aes(chisq, abs(s_0 - s_2))) + geom_point(size = 2) + ggtitle("Difference in S vs. ChiSq-Stat") + theme_classic() 
p_value_plot = ggplot(y, aes(-log10(chisq.p_val), -log10(p_val))) + geom_point(size = 2) + ggtitle("Allele freq diff p-value vs. ChiSq p-value") + theme_classic() + geom_abline(slope = 1, intercept = 0)

stat + p_value_plot

getFreqEvolDf = function(line){
  data.frame(Generation = 1:100,
             freq_0 = p_change(line$s_0, line$p_0),
             freq_2 = p_change(line$s_2, line$p_0)) %>% 
    gather(key = "Stratification", value = "Frequency", freq_0:freq_2)
}
frequency_df = ddply(y[1:4,], .(ID), getFreqEvolDf) %>%
  mutate(ID = factor(ID)) %>%
  dplyr::rename(Locus = ID)
allele_plot_v1 = ggplot(frequency_df, aes(Generation, Frequency, 
                         group = interaction(Stratification, Locus),
                         color = Locus)) + 
  geom_line() + theme_cowplot() + theme(legend.position = "none")

allele_plot_v2 = ggplot(frequency_df, aes(Generation, Frequency, 
                         group = interaction(Stratification, Locus),
                         color = Locus)) + 
  geom_line() + facet_wrap(~Locus) + theme_cowplot() + theme(legend.position = "none")
save_plot("./HS_simulation_data/plots/allele_freq_evolution_locus_v1.png", 
          allele_plot_v1, base_height = 5, base_asp = 1)  
save_plot("./HS_simulation_data/plots/allele_freq_evolution_locus_v2.png", 
          allele_plot_v2, base_height = 5, base_asp = 1, nrow = 2, ncol = 2)  

sign_change = which(sign(y$s_0) != sign(y$s_2))
frequency_df = ddply(y[sign_change[1:4],], .(ID), getFreqEvolDf) %>%
  mutate(ID = factor(ID)) %>%
  dplyr::rename(Locus = ID)
allele_plot_v3 = ggplot(frequency_df, aes(Generation, Frequency, 
                                          group = interaction(Stratification, Locus),
                                          color = Locus)) + 
  geom_line() + theme_cowplot() + theme(legend.position = "none")
save_plot("./HS_simulation_data/plots/allele_freq_evolution_locus_withSignChange_v3.png", 
          allele_plot_v3, base_height = 5, base_asp = 1)  

y %>%
  dplyr::mutate(target = paste(target_CHROM, target_POS, sep = ":"),
                partner = paste(partner_CHROM, partner_POS, sep = ":"),
                sign_change = sign(y$s_0) != sign(y$s_2)) %>%
  dplyr::rename(Gen1_freq = p_0, 
                freq_p0 = p0,
                freq_p2 = p2) %>%
  dplyr::select(ID, target, partner, chisq.p_val, sign_change, Gen1_freq, freq_p0, freq_p2) %>%
  write_csv("data/putative_epistatic_partners.csv")
