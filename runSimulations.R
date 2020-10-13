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
if(!require(yamdar)){ if (!requireNamespace("remotes", quietly = TRUE)) 
  install.packages("remotes")
  remotes::install_github("diogro/yamda-r", subdir = "package")
  ; library(yamdar)}
