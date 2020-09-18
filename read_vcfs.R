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

lt = function(x) x[lower.tri(x)]
#x = vcfs_C[[1]]
last = function(x) x[[length(x)]]
LDplot = function(x, main = NULL){
  print(1)
  n_snp = dim(x)[2]
  LDmap = LD(x, lim = c(1, n_snp))
  trim_snps = floor(seq(1, n_snp, length.out = 600))
  LDmap = LDmap[trim_snps, trim_snps]
  chr = floor(4*(x@snps$pos/100000)) + 1
  names(chr) = x@snps$pos
  trim_chr = chr[trim_snps]
  if(!any(is.na(LDmap))){
    hypot = toHypotMatrix(trim_chr)
    meanLD = TestModularity(LDmap, hypot, permutations = 1)[5,4:5]  
  }
  else {
    exclude = which(apply(LDmap, 1, function(x) sum(is.na(x)) > 2))
    trim_chr = trim_chr[-exclude]
    LDmap = LDmap[-exclude, -exclude]
    hypot = toHypotMatrix(trim_chr)
    meanLD = TestModularity(LDmap, hypot, permutations = 1)[5,4:5]   
  }
  p = superheat(LDmap, heat.lim = c(0, 1), title = main, 
                membership.rows = trim_chr, membership.cols = trim_chr, 
                left.label = "none",
                bottom.label = "none", legend = FALSE)
  return(list(LD = LDmap, chr = chr, meanLD = meanLD))
}
read_SimVCFs = function(pattern){
  vcfs_files = dir("HS_simulation_data/outputs/", pattern = pattern, full.names = T)
  vcfs = lapply(vcfs_files, read.vcf)
  vcfs = lapply(vcfs, select.snps, maf > 0.05)
  names(vcfs) = paste(rep(paste0("HS", str_pad(c(1, seq(10, 100, 10)), 3, "0", 
                                                  side = "left")), each = 3),
                         1:3, sep = "_")
  all_pos = lapply(vcfs, function(x) x@snps$pos)
  possible_pos = unique(unlist(all_pos))
  for(i in 1:length(possible_pos)){
    current_pos = possible_pos[i]
    mask = sapply(all_pos, function(x) current_pos %in% x)
    if(!all(mask)) possible_pos[i] = NA
    if(!any(mask)) print(current_pos)
  }
  possible_pos = na.omit(possible_pos)
  vcfs = lapply(vcfs, select.snps, pos %in% possible_pos) %>%
    lapply(., SNP.rm.duplicates)
  vcfs
}

vcfs_HS = read_SimVCFs("ssHS_\\d{3}_\\d{1}.vcf")
vcfs_C = read_SimVCFs("ssC_\\d{3}_\\d{1}.vcf")


pdf(file = "HS_simulation_data/plots/HS.pdf")
vcfs_HS_LD = lapply(vcfs_HS, LDplot, "HS")
dev.off()

png("plots/HS_gen001.png", height = 1080, width = 1080)
LDplot(vcfs_HS[[1]])
dev.off()
png("plots/HS_gen100.png")
LDplot(last(vcfs_HS))
dev.off()

pdf(file = "HS_simulation_data/plots/control.pdf")
vcfs_C_LD = lapply(vcfs_C, LDplot, "Control")
dev.off()
png("plots/C_gen001.png")
LDplot(vcfs_C[[1]])
dev.off()
png("plots/C_gen100.png")
LDplot(vcfs_C[[11]])
dev.off()


maf_list_C = llply(names(vcfs_C), 
                 function(x, vcfs){ 
                   out = data.frame(pos = vcfs[[x]]@snps$pos, x = vcfs[[x]]@p)
                   names(out) = c("pos", x)
                   out
                   }, 
                 vcfs_C) 
maf_list_HS = llply(names(vcfs_HS), 
                   function(x, vcfs){ 
                     out = data.frame(pos = vcfs[[x]]@snps$pos, x = vcfs[[x]]@p)
                     names(out) = c("pos", x)
                     out
                   }, 
                   vcfs_HS) 
maf_df = c(maf_list_C, maf_list_HS) %>% reduce(left_join, by = "pos")
maf_df = maf_df[complete.cases(maf_df),]
rownames(maf_df) = maf_df$pos
chr = vcfs_C_LD$C001_1$chr[rownames(maf_df)]
maf_df$chr = chr
write_csv("./HS_simulation_data/outputs/allele_freq_table_sharedSelection.csv", x = maf_df)
maf_array = maf_df
maf_array$pos = NULL
maf_array$chr = NULL
maf_array = t(maf_array)
corr_freq = cor(maf_array)
PrintMatrix(corr_freq, "./HS_simulation_data/outputs/allele_freq_corr_matrix_sharedSelection.csv")
png("./HS_simulation_data/plots/allele_freq_corr_sharredSelection.png", width = 1440, height = 1440)
superheat(corr_freq,
          membership.rows = chr, membership.cols = chr, 
          left.label = "none",
          bottom.label = "none", legend = FALSE,  grid.hline.size = 1.5,
          grid.vline.size = 1.5)
dev.off()

freq_corr_hist = data.frame(correlation = lt(corr_freq)) %>%
  ggplot(aes(correlation)) + geom_histogram(bins = 100) + theme_bw() +
  labs(x = 'Allele-frequency correlation')
save_plot("./plots/allele_freq_corr_histogram_sharredSelection.png", 
          freq_corr_hist, base_height = 6, base_asp = 1.5)

x = maf_list_HS %>% reduce(left_join, by = "pos")
x = x[complete.cases(x),]
rownames(x) = x$pos
x$pos = NULL
x = t(x)
selected_corr_freq = cor(x)
HS_freq_corr_hist = data.frame(correlation = lt(selected_corr_freq)) %>%
  ggplot(aes(correlation)) + geom_histogram(bins = 100) + theme_bw() +
  labs(x = 'Allele-frequency correlation')
save_plot("./plots/allele_freq_corr_histogram_HS.png", 
          HS_freq_corr_hist, base_height = 6, base_asp = 1.5)

x = maf_list_C %>% reduce(left_join, by = "pos")
x = x[complete.cases(x),]
rownames(x) = x$pos
x$pos = NULL
x = t(x)
C_corr_freq = cor(x)
C_freq_corr_hist = data.frame(correlation = lt(C_corr_freq)) %>%
  ggplot(aes(correlation)) + geom_histogram(bins = 100) + theme_bw() +
  labs(x = 'Allele-frequency correlation')
save_plot("./plots/allele_freq_corr_histogram_HS.png", 
          HS_freq_corr_hist, base_height = 6, base_asp = 1.5)


a_freq_pca = data.frame(loadings, generation = rownames(loadings), 
                        type = rep(c("control", "selection"), each = 33)) %>%
  mutate(replica = generation) %>% separate(replica, c("ID", "replica"), "_") %>%
  ggplot(aes(X1, X2, color = type, group = interaction(type, replica), shape = replica)) +
  geom_line() + geom_point(size = 3) +
  geom_label_repel(aes(label = ID)) +
  labs(x = "PC1", y = "PC2") + theme_bw() + ggtitle("Allele Frequency PCA")
save_plot("plots/allele_frequency_PCA_sharedSelection.png", a_freq_pca, base_height = 8, base_asp = 1.4)
