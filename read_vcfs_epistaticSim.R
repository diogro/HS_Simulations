#if(!require(pegas)){install.packages("pegas"); library(pegas)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(gaston)){install.packages("gaston"); library(gaston)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
# if(!require(animation)){install.packages("animation"); library(animation)}
if(!require(ggrepel)){install.packages("ggrepel"); library(ggrepel)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(yamdar)){remotes::install_github("diogro/yamda-r", subdir = "package")
; library(yamdar)}
if(!require(vcfR)){install.packages("vcfR"); library(vcfR)}
if(!require(latex2exp)){install.packages("latex2exp"); library(latex2exp)}

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

vcfToGenotype = function(vcf){
  x = gsub("\\|", "/", vcf@gt[,-1])
  p = nrow(x)
  df_list = lapply(1:p, function(i) genotype(x[i,], alleles = c("0", "1")))
  df_snps = as.data.frame(df_list)
  colnames(df_snps) = vcf@fix[,2]
  df_snps = df_snps[,-1]
}
LD_object = LD(df_snps)


lt = function(x) x[lower.tri(x)]

measure_LD = function(sim, experimental){
  sim_folder_epi = file.path(paste0("data/epistatic_tests/NQTL1000f_E200/SimRep", sim, "/ExpRepPlus", experimental))
  sim_folder_add = file.path(paste0("data/epistatic_tests/NQTL1000f/SimRep", sim, "/ExpRepPlus", experimental))
  
  trait_files = dir(sim_folder_add, full.names = T, pattern = "Trait.txt")
  traits = lapply(trait_files, read.table, sep = " ", header = FALSE)
  traits_df = do.call(cbind, traits)
  colnames(traits_df) = paste0("Gen", c(1, seq(10, 90, 10),100))

  epistatic_pairs = read_delim(file.path(sim_folder_epi, "_EpiQTL_list.txt"), delim = " ")
  non_epistatic_pairs = read_delim(file.path(sim_folder_epi, "_EpiQTL_list.txt"), delim = " ")
  
  gen1 = gaston::read.vcf(file.path(sim_folder_epi, "Gen_001_Sample.vcf"))
  gen1 = gaston::read.vcf(file.path(sim_folder_epi, "Gen_001_Sample.vcf"))
  
  gen100_epi = gaston::read.vcf(file.path(sim_folder_epi, "Gen_100_Sample.vcf"))
  gen100_add = gaston::read.vcf(file.path(sim_folder_add, "Gen_100_Sample.vcf"))
  
  
  LD1 = gaston::LD(gen1, c(1, nrow(gen1@snps)))
  LD100_epi = gaston::LD(gen100_epi, c(1, nrow(gen100_epi@snps)))
  LD100_add = gaston::LD(gen100_add, c(1, nrow(gen100_add@snps)))
  
  dimnames(LD1) = list(gen1@snps$pos, gen1@snps$pos)
  dimnames(LD100_epi) = list(gen100_epi@snps$pos, gen100_epi@snps$pos)
  #superheat(LD100_epi)
  
  dimnames(LD100_add) = list(gen100_add@snps$pos, gen100_add@snps$pos)
  #superheat(LD100_add)
  
  LD1 = LD1[which(gen1@snps$maf > 0), which(gen1@snps$maf > 0)]
  LD100_epi = LD100_epi[which(gen100_epi@snps$maf > 0), which(gen100_epi@snps$maf > 0)]
  LD100_add = LD100_add[which(gen100_add@snps$maf > 0), which(gen100_add@snps$maf > 0)]
  
  mask = intersect(intersect(rownames(LD1),rownames(LD100_add)),rownames(LD100_epi))
  
  LD1 = LD1[mask, mask]
  LD100_epi = LD100_epi[mask, mask]
  LD100_add = LD100_add[mask, mask]
  
  epistatic_pairs = epistatic_pairs[epistatic_pairs$pos1 %in% rownames(LD100_epi),]
  epistatic_pairs = epistatic_pairs[epistatic_pairs$pos2 %in% rownames(LD100_epi),]
  
  mean_LD1   = mean(diag(  LD1[as.character(epistatic_pairs$pos1), as.character(epistatic_pairs$pos2)]), na.rm = TRUE)
  
  epistatic_pairs = epistatic_pairs[epistatic_pairs$pos1 %in% rownames(LD100_epi),]
  epistatic_pairs = epistatic_pairs[epistatic_pairs$pos2 %in% rownames(LD100_epi),]
  mean_LD100_epi = mean(diag(LD100_epi[as.character(epistatic_pairs$pos1), as.character(epistatic_pairs$pos2)]), na.rm = TRUE)
  
  non_epistatic_pairs = non_epistatic_pairs[non_epistatic_pairs$pos1 %in% rownames(LD100_add),]
  non_epistatic_pairs = non_epistatic_pairs[non_epistatic_pairs$pos2 %in% rownames(LD100_add),]
  mean_LD100_add = mean(diag(LD100_add[as.character(non_epistatic_pairs$pos1), 
                                       as.character(non_epistatic_pairs$pos2)]), na.rm = TRUE)
  
  n_rep = 10000
  replicate_1 = numeric(n_rep)
  replicate_100_epi = numeric(n_rep)
  replicate_100_add = numeric(n_rep)
  cutoff = table(as.numeric(rownames(LD1)) > 30000000/2)[1]
  for(i in 1:n_rep){
    snp1 = sample(1:cutoff, nrow(epistatic_pairs))
    snp2 = sample(cutoff:nrow(LD1), nrow(epistatic_pairs))
    replicate_1[i]   = mean(diag(LD1[snp1, snp2]))
    replicate_100_epi[i] = mean(diag(LD100_epi[snp1, snp2]))
    replicate_100_add[i] = mean(diag(LD100_add[snp1, snp2]))
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

# x = Map(measure_LD, 1:50, sample(1:3, 50, T))
# rio::export(x, "data/LD_simulation_r2.rds")
x = rio::import("data/LD_simulation_r2.rds")

axis.font.size = 9*1.7
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
  geom_boxplot(outlier.size = 0.1, size = 0.3) + labs(x = "Generation and Scenario") + 
  ylab(TeX("Gametic Disequilibrium ($R^2$)")) +
  scale_fill_discrete(labels = c("\nAverage between\nrandom SNP pairs\n", "\nAverage between\nQTL pairs\n"), name = "") +
  scale_x_discrete(labels = c("Generation 1\nBoth", "Generation 100 \n Additive", "Generation 100 \n Epistatic")) +
   theme_classic() + theme(legend.position = "right") +
     theme(axis.title.x = element_text(size = axis.font.size*1.1), 
           axis.title.y = element_text(size = axis.font.size*1.2),
           axis.text = element_text(size = axis.font.size), 
           legend.text = element_text(size = axis.font.size))
save_plot(file= "figures/LD_boxplot_epistatic_Additive_r2.png", LD_boxplot, base_width = 1.5*6  , base_height = 1.5*3.5)

saveRDS(LD_boxplot, "figures/LD_boxplots.rds")

png("HS_simulation_data/plots/LD_truncation_selection_epistasis.png", width=15, height=15, units="in", res=300, pointsize=20)

dev.off()

