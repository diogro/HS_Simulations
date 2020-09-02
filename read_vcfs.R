#if(!require(pegas)){install.packages("pegas"); library(pegas)}
if(!require(gaston)){install.packages("gaston"); library(gaston)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
if(!require(animation)){install.packages("animation"); library(animation)}

lt = function(x) x[lower.tri(x)]
LDplot = function(x, main){
  n_snp = dim(x)[2]
  LDmap = LD(x, lim = c(1, n_snp))
  trim_snps = floor(seq(1, n_snp, length.out = 400))
  LDmap = LDmap[trim_snps, trim_snps]
  p = superheat(sqrt(LDmap), heat.lim = c(0, 1), title = main)
  return(list(p, LDmap))
}

vcfs_files = dir("outputs/", pattern = "HS_\\d{3}.vcf",full.names = T)
vcfs = lapply(vcfs_files, read.vcf)
vcfs = lapply(vcfs, select.snps, maf > 0.05)
all_pos = lapply(vcfs, function(x) x@snps$pos)
possible_pos = unique(unlist(all_pos))
possible_pos = possible_pos[sapply(possible_pos, function(x) all(sapply(all_pos, function(y) x %in% y)))]
vcfs = lapply(vcfs, select.snps, pos %in% possible_pos)
lapply(vcfs, LDplot, "HS")


vcfs_files = dir("outputs/", pattern = "C_\\d{3}.vcf",full.names = T)
vcfs = lapply(vcfs_files, read.vcf)
vcfs = lapply(vcfs, select.snps, maf > 0.05)
all_pos = lapply(vcfs, function(x) x@snps$pos)
possible_pos = unique(unlist(all_pos))
possible_pos = possible_pos[sapply(possible_pos, function(x) all(sapply(all_pos, function(y) x %in% y)))]
vcfs = lapply(vcfs, select.snps, pos %in% possible_pos)
lapply(vcfs, LDplot, "Control")
