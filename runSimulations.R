if(!require(rhdf5)){
  if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")
  BiocManager::install("rhdf5"); library(rhdf5)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(gaston)){install.packages("gaston"); library(gaston)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
if(!require(animation)){install.packages("animation"); library(animation)}
if(!require(ggrepel)){install.packages("ggrepel"); library(ggrepel)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(fs)){install.packages("fs"); library(fs)}
if(!require(yamdar)){ if (!requireNamespace("remotes", quietly = TRUE)) 
  install.packages("remotes")
  remotes::install_github("diogro/yamda-r", subdir = "package")
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
folder = out_folder
pattern = "ssC_\\d{3}_\\d{1}.vcf"
type = "C"
read_SimVCFs = function(pattern, folder, type){
  vcfs_files = dir(folder, pattern = pattern, full.names = T)
  vcfs = lapply(vcfs_files, read.vcf)
  vcfs = lapply(vcfs, select.snps, maf > 0.05)
  names(vcfs) = paste(rep(paste0(type, str_pad(c(1, seq(10, 100, 10)), 3, "0", 
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

makeFreqTable = function(folder, type, n_pop){
  mutations = read_delim(file.path(folder, paste0("ss", type, ".mut")), delim = " ", 
             col_names = c("out", "gen", "tracked", "pop", "id", "mut", 
                           "pos", "s", "d", "origin_pop", "origin_gen", 
                           "prevalence")) %>% 
    mutate(gen = gen - 10100,
           type = type,
           pop = gsub("p2", "", pop),
           freq = prevalence/n_pop) %>%
    select(id, type, gen, pop, mut, pos, s, freq) 

  mutations
}
makeFreqPCA = function(freq){
  snp_list = freq %>%
    select(id, type, gen, pos, freq, pop) %>%
    dlply(.(id)) 
  mask <- sapply(snp_list, nrow)
  mask = mask == max(mask)
  c_snps = snp_list[mask]
  gen_id = paste(c_snps[[1]]$type, c_snps[[1]]$gen, c_snps[[1]]$pop, sep = "_")
  gen_df = data.frame(gen_id, type = c_snps[[1]]$type, gen = c_snps[[1]]$gen, rep = c_snps[[1]]$pop)
  freq_array = t(laply(c_snps, function(x) x$freq))
  eVec = prcomp(freq_array)
  loadings = freq_array %*% eVec$rotation[,1:2]
  freq_pca = data.frame(gen_df, loadings)
  freq_pca$gen = factor(freq_pca$gen, levels = c(0, seq(10, 100, 10)))
  a_freq_pca = freq_pca %>%
    ggplot(aes(PC1, PC2, color = type, group = interaction(type, rep), shape = rep)) +
    geom_point(size = 3) +
    geom_label_repel(aes(label = gen)) +
    labs(x = "PC1", y = "PC2") + theme_bw() + ggtitle("Allele Frequency PCA")
  print(a_freq_pca)
  list(pca = freq_pca, plot = a_freq_pca)
}

call_slim = function(seed = NULL, out = NULL, 
                     n_founders = 20, n_pop = 1000, n_sample = 100, 
                     s_gauss = 1, sel_var = 0.1, shared = "T",
                     selected_fraction = 0.2){
  sim_name = file.path(paste0("sSD-", sel_var, 
                   "_s_gauss-", s_gauss,
                   "_shared-", shared,
                   "_nsample-", n_sample,
                   "_npop-", n_pop,
                   "_nfounders-", n_founders,
                   "_seed-", seed))
  if(is.null(out))
    out = sim_name
  else 
    out = file.path(out, sim_name) 
  out_folder = file.path("~/projects/HS_simulations/HS_simulation_data/outputs/", out)
  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
  
  args = list(seed, n_founders, n_pop, n_sample, s_gauss, sel_var, shared,
              selected_fraction)
  names(args) = c("seed", "n_founders", "n_pop", "n_sample", "s_gauss",
                  "sel_var", "shared", "selected_fraction")
  mask = sapply(args, is.null)
  args = args[!mask]
  
  args["out"] = digest::digest2int(digest::digest(args, "md5"))
  slim_out = paste0("~/projects/HS_simulations/HS_simulation_data/outputs/", args["out"])
  dir.create(slim_out, showWarnings = FALSE)
  
  call_control = paste0("slim -d control=T ", 
               paste("-d ", names(args), "=", args, " ", sep = "", collapse = ""),
               "sharedSelection.slim")
  system(call_control)
  
  call_HS = paste0("slim -d control=F ", 
                paste("-d ", names(args), "=", args, " ", sep = "", collapse = ""),
                "sharedSelection.slim")
  system(call_HS)
  write.table(as.data.frame(args), file.path(slim_out, "args.tsv"), 
              sep = "\t", row.names = F)
  
  for(i in dir(slim_out, full.names = T, recursive = T))
    fs::file_move(i, out_folder)
  file_delete(slim_out)
  return(out_folder)
}

out_folder = gaussian_out[[i]]
readSimulation = function(out_folder){
  vcfs_HS = read_SimVCFs(pattern = "ssHS_\\d{3}_\\d{1}.vcf", folder = out_folder , type = "HS")
  vcfs_C = read_SimVCFs(pattern = "ssC_\\d{3}_\\d{1}.vcf", folder = out_folder, "C")
  vcfs = list(HS = vcfs_HS, C = vcfs_C)
  
  n_pop = as.numeric(strsplit(grep("npop", strsplit(out_folder, "_")[[1]], value = T), "-")[[1]][2])
  freq_table_C = makeFreqTable(out_folder, "C", n_pop)
  freq_table_HS = makeFreqTable(out_folder, "HS", n_pop)
  freq_table = rbind(freq_table_C, freq_table_HS)
  pca = makeFreqPCA(freq_table)
  list(vcfs = vcfs, freq = freq_table, pca = pca$pca, plot_pca = pca$plot, folder = out_folder)
}

gaussian_out = vector("list", 10)
gaussian_sim = vector("list", 10)

s_var = seq(0.01, 0.1, length.out = 10)
for(i in 1:10){
  #gaussian_out[[i]] = call_slim(seed = i, shared = "T", s_gauss=1, sel_var = s_var[i], out = "GaussianSelection")
  gaussian_sim[[i]] = readSimulation(gaussian_out[[i]])
}

p1 = gaussian_sim[[8]]$plot_pca + ggtitle("Allele Frequency PCA - Strong Selection - Var_s = 0.08")
save_plot("gaussian_effects-Var_s-0.08.png", p1, base_height = 7)
p2 = gaussian_sim[[1]]$plot_pca + ggtitle("Allele Frequency PCA - Weak Selection - Var_s = 0.01")
save_plot("gaussian_effects-Var_s-0.01.png", p2, base_height = 7)

beta_out = vector("list", 10)
beta_sim = vector("list", 10)

s_a = seq(50, 150, length.out = 10)
for(i in 1:10){
  beta_out[[i]] = call_slim(seed = i, shared = "T", s_gauss=0, sel_var = s_a[i], out = "BetaSelection")
  beta_sim[[i]] = readSimulation(beta_out[[i]])
}
p1 = beta_sim[[1]]$plot_pca + ggtitle("Allele Frequency PCA - Strong Selection - Beta(50, 1)")
save_plot("beta_effects-Beta_50.png", p1, base_height = 7)
p2 = beta_sim[[10]]$plot_pca + ggtitle("Allele Frequency PCA - Weak Selection - Beta(150, 1)")
save_plot("beta_effects-Beta_150.png", p2, base_height = 7)


#out_folder = "./HS_simulation_data/outputs/test/sSD-0.01_shared-T_nsample-100_npop-1000_nfounders-20_seed-1/"
sim_out = readSimulation(beta_out[[10]])
sim_out = readSimulation(gaussian_out[[10]])

snp_list = sim_out$freq %>%
  select(id, type, gen, pos, freq, pop) %>%
  dlply(.(id)) 

