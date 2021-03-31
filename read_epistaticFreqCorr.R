#if(!require(pegas)){install.packages("pegas"); library(pegas)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
if(!require(animation)){install.packages("animation"); library(animation)}
if(!require(ggrepel)){install.packages("ggrepel"); library(ggrepel)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(yamdar)){remotes::install_github("diogro/yamda-r", subdir = "package")
; library(yamdar)}

file = QTL_files_add[[1]]
read_mutation = function(file){
  read_delim(file, 
           col_names = c("out", "gen", "tracked", 
                         "pop", "id", "type", 
                         "pos", "ad", "dm", "pop_org",
                         "gen_org", "prev"),
           delim = " ") %>% 
    select(gen, pos, prev) %>% unique %>%
    mutate(gen = gen - 49999,
           freq = prev / 10000,
           pos = pos + 1,
           epi = ifelse(pos %in% c(epistatic_pairs$pos1, 
                                   epistatic_pairs$pos2), 1, 0)) %>%
    arrange(pos)
}

sim = 41
experimental = 1
measure_CorrPairs = function(sim, experimental){
  print(sim)
  print(experimental)
  sim_folder_epi = file.path(paste0("data/epistatic_tests/NQTL1000f_E200/SimRep", sim, "/ExpRepPlus", experimental))
  sim_folder_add = file.path(paste0("data/epistatic_tests/NQTL1000f/SimRep", sim, "/ExpRepPlus", experimental))
  
  QTL_files_add = dir(sim_folder_add, full.names = T, pattern = "mutations.txt")
  QTL_files_epi = dir(sim_folder_epi, full.names = T, pattern = "mutations.txt")
  
  epistatic_pairs = read_delim(file.path(sim_folder_epi, "_EpiQTL_list.txt"), delim = " ")
  non_epistatic_pairs = read_delim(file.path(sim_folder_epi, "_EpiQTL_list.txt"), delim = " ")
  
  mutation_table_add = ldply(QTL_files_add, read_mutation)
  mutation_table_epi = ldply(QTL_files_epi, read_mutation)
  
  tb_add = table(mutation_table_add$pos, mutation_table_add$gen)
  if(length(which(tb_add > 1)) > 0){
    bad_pos = as.numeric(rownames(tb_add)[which(tb_add[,1] > 1)])
    mutation_table_add = mutation_table_add[mutation_table_add$pos!=bad_pos ,]
    epistatic_pairs = epistatic_pairs[epistatic_pairs$pos1 != bad_pos & epistatic_pairs$pos2 != bad_pos, ] 
  }
  
  tb_epi = table(mutation_table_epi$pos, mutation_table_epi$gen)
  if(length(which(tb_epi > 1)) > 0){
    bad_pos = as.numeric(rownames(tb_epi)[which(tb_epi[,1] > 1)])
    mutation_table_epi = mutation_table_epi[mutation_table_epi$pos!=bad_pos ,]
    epistatic_pairs = epistatic_pairs[epistatic_pairs$pos1 != bad_pos & epistatic_pairs$pos2 != bad_pos, ] 
  }
  
  mutation_wide_add = mutation_table_add %>%
    select(gen, pos, freq) %>%
    spread(pos, value = freq, drop = FALSE, fill = NA) %>%
    select(-gen)
  
  mutation_wide_epi = mutation_table_epi %>%
    select(gen, pos, freq) %>%
    spread(pos, value = freq, drop = FALSE, fill = NA) %>%
    select(-gen)
  
  corr_mat_add = cor(mutation_wide_add)
  corr_mat_epi = cor(mutation_wide_epi)
  
  mean_corr_add = mean(diag(corr_mat_add[as.character(epistatic_pairs$pos1), 
                         as.character(epistatic_pairs$pos2)]), na.rm = TRUE)
  
  mean_corr_epi = mean(diag(corr_mat_epi[as.character(epistatic_pairs$pos1), 
                         as.character(epistatic_pairs$pos2)]), na.rm = TRUE)
  
  c(add = mean_corr_add, 
    epi = mean_corr_epi)
}

x = Map(measure_CorrPairs, 1:50, sample(1:3, 50, T))
mean_corrs = data.frame(do.call(rbind, x)) %>%
  reshape2::melt()

FreqCorr_boxplot = ggplot(mean_corrs, aes(variable, value, group = variable, fill = variable)) + 
  geom_boxplot() + labs(y = "Allele frequency correlation", x = "Scenario") + 
  scale_fill_discrete(labels = c("Average correlation between\nadditive QTLs", 
                                 "Average correlation between\nepistatic QTLs"), name = "") +
  scale_x_discrete(labels = c("Additive Scenario", "Epistatic Scenario")) +
  theme_cowplot() + theme(legend.position = "bottom")
save_plot(file= "HS_simulation_data/plots/FreqCorrelation_boxplot_epistatic_Additive.png", FreqCorr_boxplot, base_height = 8, base_asp = 1.2)

png("HS_simulation_data/plots/LD_truncation_selection_epistasis.png", width=15, height=15, units="in", res=300, pointsize=20)

dev.off()

