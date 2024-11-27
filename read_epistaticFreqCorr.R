#if(!require(pegas)){install.packages("pegas"); library(pegas)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
#if(!require(animation)){install.packages("animation"); library(animation)}
if(!require(ggrepel)){install.packages("ggrepel"); library(ggrepel)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
if(!require(ggthemes)){install.packages("ggthemes"); library(ggthemes)}

if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(yamdar)){remotes::install_github("diogro/yamda-r", subdir = "package")
; library(yamdar)}

file = QTL_files_add[[1]]
read_mutation = function(file){
  read_delim(file, 
           col_names = c("nothing", "out", "gen", "tracked", 
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
  non_epistatic_pairs = read_delim(file.path(sim_folder_add, "_EpiQTL_list.txt"), delim = " ")
  
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

scale = 1.4
axis.font.size = 9*scale
FreqCorr_boxplot = ggplot(mean_corrs, aes(variable, value, group = variable)) + 
   geom_boxplot(outlier.size = 0.5, size = 0.5 ) + 
   labs(y = "Allele frequency correlation\n between unlinked QTL pairs", x = "Scenario") + 
  scale_fill_manual(values=c("darkorange2", "purple"),
                    labels = c("Average correlation between\nadditive QTLs", 
                               "Average correlation between\nepistatic QTLs"), name = "") +
  scale_x_discrete(labels = c("Additive", "Epistatic")) +
  theme_classic() + theme(legend.position = "none") +
     theme(axis.title.x = element_text(size = axis.font.size*1.2), 
           axis.title.y = element_text(size = axis.font.size*1.),
           axis.text = element_text(size = axis.font.size))
save_plot(file= "figures/FreqCorrelation_boxplot_epistatic_Additive.png", 
          FreqCorr_boxplot,  base_width = scale*5  , base_height = scale*2.7)

axis.font.size = 13
FreqCorr_boxplot_poster = FreqCorr_boxplot + theme(legend.position = "none") +
     theme(axis.title.x = element_text(size = axis.font.size*1.2), 
           axis.title.y = element_text(size = axis.font.size*1.),
           axis.text = element_text(size = axis.font.size))

library(ggproto)

LD_boxplot = rio::import("figures/LD_boxplots.rds")
LD_boxplot_poster = LD_boxplot +
     theme(axis.title.x = element_text(size = axis.font.size), 
           axis.title.y = element_text(size = axis.font.size),
           axis.text = element_text(size = axis.font.size*0.8),
           legend.text = element_text(size = axis.font.size*0.8))

library(patchwork)

poster_panels = FreqCorr_boxplot_poster + LD_boxplot_poster + 
  plot_layout(widths = c(2, 3))
save_plot(file= "figures/FreqCorrelation_boxplot_epistatic_Additive_poster.png", 
          poster_panels,  base_width =  28/2.54 , base_height = 10/2.54)
