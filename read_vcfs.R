#if(!require(pegas)){install.packages("pegas"); library(pegas)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(gaston)){install.packages("gaston"); library(gaston)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
if(!require(animation)){install.packages("animation"); library(animation)}
if(!require(ggrepel)){install.packages("ggrepel"); library(ggrepel)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}

lt = function(x) x[lower.tri(x)]
LDplot = function(x, main = NULL){
  n_snp = dim(x)[2]
  LDmap = LD(x, lim = c(1, n_snp))
  trim_snps = floor(seq(1, n_snp, length.out = 600))
  LDmap = LDmap[trim_snps, trim_snps]
  chr = floor(4*(x@snps$pos/100000)) + 1
  chr = chr[trim_snps]
  p = superheat(LDmap, heat.lim = c(0, 1), title = main, 
                membership.rows = chr, membership.cols = chr, 
                left.label = "none",
                bottom.label = "none", legend = FALSE)
  return(list(p, LDmap))
}

vcfs_files_HS = dir("outputs/", pattern = "HS_\\d{3}.vcf",full.names = T)
vcfs_HS = lapply(vcfs_files_HS, read.vcf)
vcfs_HS = lapply(vcfs_HS, select.snps, maf > 0.05)
names(vcfs_HS) = paste0("HS", str_pad(c(1, seq(10, 100, 10)), 3, "0", side = "left"))
all_pos = lapply(vcfs_HS, function(x) x@snps$pos)
possible_pos = unique(unlist(all_pos))
for(i in 1:length(possible_pos)){
  current_pos = possible_pos[i]
  mask = sapply(all_pos, function(x) current_pos %in% x)
  if(!all(mask)) possible_pos[i] = NA
  if(!any(mask)) print(current_pos)
}
possible_pos = na.omit(possible_pos)
vcfs_HS = lapply(vcfs_HS, select.snps, pos %in% possible_pos) %>%
  lapply(., SNP.rm.duplicates)
pdf(file = "HS.pdf")
lapply(vcfs_HS, LDplot, "HS")
dev.off()

png("plots/HS_gen001.png", height = 1080, width = 1080)
LDplot(vcfs_HS[[1]])
dev.off()
png("plots/HS_gen100.png")
LDplot(vcfs_HS[[11]])
dev.off()

vcfs_files_C = dir("outputs/", pattern = "C_\\d{3}.vcf",full.names = T)
vcfs_C = lapply(vcfs_files_C, read.vcf)
names(vcfs_C) = paste0("C", str_pad(c(1, seq(10, 100, 10)), 3, "0", side = "left"))
vcfs_C = lapply(vcfs_C, select.snps, maf > 0.05)
all_pos = lapply(vcfs_C, function(x) x@snps$pos)
possible_pos = unique(unlist(all_pos))
for(i in 1:length(possible_pos)){
  current_pos = possible_pos[i]
  mask = sapply(all_pos, function(x) current_pos %in% x)
  if(!all(mask)) possible_pos[i] = NA
  if(!any(mask)) print(current_pos)
}
possible_pos = na.omit(possible_pos)
vcfs_C = lapply(vcfs_C, select.snps, pos %in% possible_pos) %>%
  lapply(., SNP.rm.duplicates)
pdf(file = "control.pdf")
lapply(vcfs_C, LDplot, "Control")
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
x = c(maf_list_C, maf_list_HS) %>% reduce(left_join, by = "pos")
x = x[complete.cases(x),]
rownames(x) = x$pos
x$pos = NULL
x = t(x)
eigen_Freq = eigen(cov(x))
loadings = (x %*% eigen_Freq$vectors)[,1:2]
plot(loadings)
a_freq_pca = data.frame(loadings, generation = rownames(loadings), type = rep(c("control", "selection"), each = 11)) %>%
  ggplot(aes(X1, X2, color = type)) + geom_point() + geom_label_repel(aes(label = generation)) +
  labs(x = "PC1", y = "PC2") + theme_bw() + ggtitle("Allele Frequency PCA")
save_plot("plots/allele_frequency_PCA.png", a_freq_pca, base_height = 6, base_asp = 1.4)
