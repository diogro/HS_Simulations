#if(!require(pegas)){install.packages("pegas"); library(pegas)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(gaston)){install.packages("gaston"); library(gaston)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
if(!require(animation)){install.packages("animation"); library(animation)}
if(!require(ggrepel)){install.packages("ggrepel"); library(ggrepel)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(yamdar)){remotes::install_github("diogro/yamda-r", subdir = "package")
; library(yamdar)}
if(!require(vcfR)){install.packages("vcfR"); library(vcfR)}
if(!require(genetics)){install.packages("genetics"); library(genetics)}

# bed.mat = gen1
# lim = c(1, 100)
# lim2 = NULL
LD.chisq = function(bed.mat, lim, lim2 = NULL, stat = FALSE,...){
  if(is.null(lim2)) lim2 = lim
  x = as.matrix(bed.mat)
  n = nrow(x)
  p = ncol(x)
  LD = array(NA, dim = c(lim[2]-lim[1]+1, lim2[2]-lim2[1]+1))
  for(snp1 in lim[1]:lim[2])
    for(snp2 in lim2[1]:lim2[2]){
      tab <- table(x[, c(snp1, snp2)])
      if(length(tab) > 1){
        csq_t = chisq.test(tab)
        test = ifelse(stat, csq_t$statistic, -log10(csq_t$p.value))
      }
      else test = NA
      LD[snp1 - lim[1] + 1, snp2 - lim2[1] + 1] = test
    }
  LD
}

# vcf = vcfR::read.vcfR(file.path(sim_folder_epi, "Gen_001_Sample.vcf"))
vcfToGenotype = function(vcf){
  x = gsub("\\|", "/", vcf@gt[,-1])
  p = nrow(x)
  df_list = lapply(1:p, function(i) genotype(x[i,], alleles = c("0", "1")))
  df_snps = as.data.frame(df_list)
  colnames(df_snps) = vcf@fix[,2]
  df_snps = df_snps[,-1]
}
vcfToNumeric = function(vcf){
  x = gsub("1\\|1", "2", vcf@gt[,-1]) 
  x = gsub("0\\|0", "0", x) 
  x = gsub("1\\|0", "1", x) 
  x = gsub("0\\|1", "1", x) 
  p = nrow(x)
  df_list = lapply(1:p, function(i) x[i,])
  df_snps = as.data.frame(df_list)
  colnames(df_snps) = vcf@fix[,2]
  df_snps = df_snps[,-1]
}

# pairLD = function(geno, snp1, snp2){
#   LD_object = LD(geno[,c(snp1, snp2)])
#   LD_object$`X^2`[1,2]
# }
# geno = gen1
# snp1 = as.character(epistatic_pairs$pos1[1])
# snp2 = as.character(epistatic_pairs$pos2[1])

pairLD = function(geno, snp1, snp2){
  tab <- table(geno[, c(snp1, snp2)])
  if(length(tab) < 2 )
    return(NA)
  LD_object = chisq.test(tab)
  -log10(LD_object$p.value)
}
mean_pairLD_list = function(geno, snp_list1, snp_list2){
  mean(unlist(Map(function(x, y) pairLD(geno, x, y),
                  snp_list1, 
                  snp_list2)), na.rm = TRUE)
}
lt = function(x) x[lower.tri(x)]
#x = vcfs_C[[1]]
# sim = 1
# experimental = 1
measure_LD = function(sim, experimental){

  sim_folder_epi = file.path(paste0("data/epistatic_tests/NQTL1000f_E200/SimRep", sim, "/ExpRepPlus", experimental))
  sim_folder_add = file.path(paste0("data/epistatic_tests/NQTL1000f/SimRep", sim, "/ExpRepPlus", experimental))
  
  trait_files = dir(sim_folder_add, full.names = T, pattern = "Trait.txt")
  traits = lapply(trait_files, read.table, sep = " ", header = FALSE)
  traits_df = do.call(cbind, traits)
  colnames(traits_df) = paste0("Gen", c(1, seq(10, 90, 10),100))

  epistatic_pairs = read_delim(file.path(sim_folder_epi, "_EpiQTL_list.txt"), delim = " ")
  non_epistatic_pairs = read_delim(file.path(sim_folder_epi, "_EpiQTL_list.txt"), delim = " ")
  
  gen1 = vcfToNumeric(vcfR::read.vcfR(file.path(sim_folder_epi, "Gen_001_Sample.vcf")))
  gen100_epi = vcfToNumeric(vcfR::read.vcfR(file.path(sim_folder_epi, "Gen_100_Sample.vcf")))
  gen100_add = vcfToNumeric(vcfR::read.vcfR(file.path(sim_folder_add, "Gen_100_Sample.vcf")))
  cat(sim)
  mask = intersect(intersect(colnames(gen1),colnames(gen100_epi)),colnames(gen100_add))
  
  gen1 = gen1[mask]
  gen100_epi = gen100_epi[mask]
  gen100_add = gen100_add[mask]
  
  epistatic_pairs = epistatic_pairs[epistatic_pairs$pos1 %in% colnames(gen1),]
  epistatic_pairs = epistatic_pairs[epistatic_pairs$pos2 %in% colnames(gen1),]
  mean_LD1 = mean_pairLD_list(gen1, 
                              as.character(epistatic_pairs$pos1), 
                              as.character(epistatic_pairs$pos2))
  
  epistatic_pairs = epistatic_pairs[epistatic_pairs$pos1 %in% colnames(gen100_epi),]
  epistatic_pairs = epistatic_pairs[epistatic_pairs$pos2 %in% colnames(gen100_epi),]
  mean_LD100_epi = mean_pairLD_list(gen100_epi, 
                                    as.character(epistatic_pairs$pos1), 
                                    as.character(epistatic_pairs$pos2))
  
  non_epistatic_pairs = non_epistatic_pairs[non_epistatic_pairs$pos1 %in% colnames(gen100_add),]
  non_epistatic_pairs = non_epistatic_pairs[non_epistatic_pairs$pos2 %in% colnames(gen100_add),]
  mean_LD100_add = mean_pairLD_list(gen100_add, 
                                    as.character(epistatic_pairs$pos1), 
                                    as.character(epistatic_pairs$pos2))
  
  n_rep = 10000
  replicate_1 = numeric(n_rep)
  replicate_100_epi = numeric(n_rep)
  replicate_100_add = numeric(n_rep)
  cutoff = table(as.numeric(colnames(gen1)) > 30000000/2)[1]
  n_snps = ncol(gen1)
  for(i in 1:n_rep){
    # cat(i)
    snp1 = sample(1:cutoff, nrow(epistatic_pairs))
    snp2 = sample(cutoff:n_snps, nrow(epistatic_pairs))
    replicate_1[i]       = mean_pairLD_list(gen1,       snp1, snp2)
    replicate_100_epi[i] = mean_pairLD_list(gen100_epi, snp1, snp2)
    replicate_100_add[i] = mean_pairLD_list(gen100_add, snp1, snp2)
    if(is.na(replicate_1[i]) | is.na(replicate_100_epi[i]) | is.na(replicate_100_add[i])){
      stop("Deu merda.")
      break
    } 
  }
  list(mean_LD1=mean_LD1, 
       mean_LD100_epi=mean_LD100_epi, 
       mean_LD100_add=mean_LD100_add, 
       replicate_1=replicate_1, 
       replicate_100_epi=replicate_100_epi, 
       replicate_100_add=replicate_100_add,
       traits_df=traits_df)
}
plotHist = function (x){
  par(mfrow = c(3, 1))
  xlim = c(min(c(x$replicate_1, x$replicate_100_epi, x$replicate_100_add, x$mean_LD1, x$mean_LD100_add, x$mean_LD100_epi)) - 0.001, 
           max(c(x$replicate_1, x$replicate_100_epi, x$replicate_100_add, x$mean_LD1, x$mean_LD100_add, x$mean_LD100_epi)) + 0.001)
  hist(x$replicate_1, 
       xlim = xlim, 
       breaks = 100, 
       main = "Distribution of LD between random pairs of loci in different chr\nGeneration 1", xlab = "LD")
  abline(v = x$mean_LD1, pch = 19, col = "red", lwd = 2)
  hist(x$replicate_100_epi, 
       xlim = xlim, breaks = 100, 
       main = "Generation 100 - Epistatic", xlab = "LD")
  mtext("LD between epistatic pairs", side = 3, col =  "red", lwd =2, adj = 0.64)
  abline(v = x$mean_LD100_epi, pch = 19, col = "red", lwd = 2)
  hist(x$replicate_100_add, xlim = xlim, breaks = 100, main = "Generation 100 - Additive", xlab = "LD")
  mtext("LD between non-epistatic pairs", side = 3, col =  "red", lwd =2, adj = 0.48)
  abline(v = x$mean_LD100_add, pch = 19, col = "red", lwd = 2)
}

#x = Map(measure_LD, 1:50, sample(1:3, 50, T))
#saveRDS(x, file = "data/LD_simulation_ChiSqPvalue.rds")
x = readRDS(file = "data/LD_simulation_ChiSqPvalue.rds")
# plotHist(x[[2]])

LD_boxplot = ldply(x, function(x) data.frame(value = c(x$mean_LD1, x$mean_LD100_epi, x$mean_LD100_add,
                                          x$replicate_1,  
                                          x$replicate_100_epi,
                                          x$replicate_100_add),
                                key = c("start", "epi", "add", rep(c("start", "epi", "add"), each = length(x$replicate_1))),
                                group = c(rep("Pairs", 3), 
                                          rep("Average", 3*length(x$replicate_1)))
                                )) %>%
  mutate(key = factor(key, levels = c("start", "add", "epi"))) %>%
  ggplot(aes(key, value, group = interaction(group, key), fill = group)) + 
  geom_boxplot() + labs(y = "LD (X^2 -log10(p-value))", x = "Generation and Scenario") + 
  scale_fill_discrete(labels = c("Average between random SNP pairs", "Average between QTL pairs"), name = "") +
  scale_x_discrete(labels = c("Starting Generation", "Generation 100 \n Additive Scenario", "Generation 100 \n Epistatic Scenario")) +
  theme_cowplot() + theme(legend.position = "bottom")
save_plot(file= "HS_simulation_data/plots/LD_boxplot_epistatic_Additive_ChiSqP-value.svg", LD_boxplot, base_height = 8, base_asp = 1.2)

# png("HS_simulation_data/plots/LD_truncation_selection_epistasis.png", width=15, height=15, units="in", res=300, pointsize=20)
# 
# dev.off()

