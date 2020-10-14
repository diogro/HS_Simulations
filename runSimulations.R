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

makeFreqTable = function(folder, type){
  mutations = read_delim(file.path(folder, paste0("ss", type, ".mut")), delim = " ", 
             col_names = c("out", "gen", "tracked", "pop", "id", "mut", 
                           "pos", "s", "d", "origin_pop", "origin_gen", 
                           "prevalence")) %>% 
    select(id, gen, pop, mut, pos, s, prevalence) %>%
    mutate(gen = gen - 10100,
           type = type,
           pop = gsub("p2", "", pop))
  mutations
}

call_slim = function(seed = NULL, out = NULL, 
                     n_founders = 20, n_pop = 1000, n_sample = 100, 
                     sel_var = 0.1, shared = "T",
                     selected_fraction = 0.2){
  sim_name = file.path(paste0("sSD-", sel_var, 
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
  
  args = list(seed, n_founders, n_pop, n_sample, sel_var, shared,
              selected_fraction)
  names(args) = c("seed", "n_founders", "n_pop", "n_sample", 
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
  
  vcfs_HS = read_SimVCFs(out_folder, pattern = "ssHS_\\d{3}_\\d{1}.vcf", "HS")
  vcfs_C = read_SimVCFs(out_folder, pattern = "ssC_\\d{3}_\\d{1}.vcf", "C")
  vcfs = list(HS = vcfs_HS, C = vcfs_C)
  
  freq_table_C = makeFreqTable(out_folder, "C")
  freq_table_HS = makeFreqTable(out_folder, "HS")
  freq_table = rbind(freq_table_C, freq_table_HS)
  
  list(vcfs = vcfs, freq = freq_table)
}
slim_out = "./HS_simulation_data/outputs/test/sSD-0.01_shared-T_nsample-1000_npop-1000_nfounders-20_seed-1/"

out = call_slim(seed = 1, shared = "T", sel_var = 0.01, out = "test")
out$freq %>%
  dlply(.(id))
