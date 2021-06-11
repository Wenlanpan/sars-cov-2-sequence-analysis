## Genomic characteristics analysis
# 1. Top mutation analysis
# 2. Co-mutation analysis
# 3. Comparative analysis
# 4. Evolutionary analysis
# Phylogenetic analysis
## 1. Phylogenetic tree analysis
# 2. Clade temporal shift annalysis

# set the working directory to the data file
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/2021spring/FYP/sars-cov-2-sequence-analysis")

# required packages
library(tidyverse)
library(lubridate)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(ggvenn)
library(segmented)
library(ggeasy)

#  prepare your datasets download from GISAID and Nextclade
source('dataset_pre.R')
all_information("la", "data/GISAID/LA_0611_2025.metadata.tsv", 
                "data/Nextclade/LA_0611_nextclade.tsv")
all_information("global", "data/GISAID/Global_0611_4001.metadata.tsv", 
                "data/Nextclade/Global_0611_nextclade.tsv")
all_information("kings", "data/GISAID/Kings_0611_569.metadata.tsv", 
                "data/Nextclade/Kings_0611_nextclade.tsv")
all_information("cook", "data/GISAID/Cook_0611_1083.metadata.tsv", 
                "data/Nextclade/Cook_0611_nextclade.tsv")

# Read in prepared datasets
la <- read_tsv("data/la2025.tsv") 
global <- read_tsv("data/global3984.tsv")
kings <- read_tsv("data/kings569.tsv") 
cook <- read_tsv("data/cook1083.tsv") 

la_mut <- read_tsv("data/la2025mut.tsv")
global_mut <- read_tsv("data/global3984mut.tsv") 
kings_mut <- read_tsv("data/kings569mut.tsv") 
cook_mut <- read_tsv("data/cook1083mut.tsv")

### Genomic characteristics analysis
source('genomic_char_analysis.R')
genomic_char_la <- genomic_char_analysis(la, la_mut)

## 1. Top mutation analysis
source('top_mut_analysis.R')
top_mut_la <- top_mut_analysis(la, la_mut) # table 1
top_mut_fig_la <- top_mut_analysis_fig(la, la_mut) # default label top 10
# top_mut_analysis_fig(la, la_mut, 5) label top 5

top_mut_fig_nt_la <- top_mut_analysis_fig_nt(la, la_mut) # default label top 10
# error: frequence is different between position and substitution 
top_mut_fig_nt_gl <- top_mut_analysis_fig_nt(global, global_mut) # label nucleotide

# figure 1
top_mut_fig1_la <- top_mut_analysis_fig_la(la, la_mut)

## 2. Co-mutation analysis
source('co-mut_analysis.R')
co_mut_la <- co_mut_analysis(la, la_mut) # default label top 10
# co_mut_analysis(la, la_mut, 5) label top 5

#- co-exist: A23403G, C14408T, C3037T, and C241T 
coexist1 <- la_mut %>% filter(substitutions %in% 
                               c("A23403G", "C14408T", "C3037T", "C241T")) %>%
  group_by(strain) %>% summarize(count = n()) %>% filter(count == 4) # 1934
#- co-exist: C1059T and G25563T
coexist2 <- la_mut %>% filter(substitutions %in% c("C1059T", "G25563T")) %>%
  group_by(strain) %>% summarize(count = n()) %>% filter(count == 2) # 1070
#- co-exist: C26681T and T22917G
coexist3 <- la_mut %>% filter(substitutions %in% c("C26681T", "T22917G")) %>%
  group_by(strain) %>% summarize(count = n()) %>% filter(count == 2) # 502

## 3. Comparative analysis
source('comp_analysis.R')
# figure 3
comp_la <- comp_analysis ("Los Angeles", la, la_mut, global, global_mut)
# comp_analysis ("Los Angeles", la, la_mut, global, global_mut, 5)

# table 2
comp_tab_la <- comp_analysis_tab("Los Angeles", la, la_mut, kings, kings_mut,
                              cook, cook_mut, global, global_mut)
#- chi square
a <- 1983
c <- 3757
b <- 2025
d <- 3984
x <- matrix(c(a, c, b-a, d-c), nrow = 2, ncol = 2)
chisq.test(x)

# figure 4
venn_diag_la <- venn_diag("Los Angeles (California)", la_mut, 
                          "Kings (New York)", kings_mut, 
                          "Cook (Illinois)", cook_mut)
#- the percentage of mutations from only a region
uniq_per_la <- unique_percent (la, la_mut, kings_mut, cook_mut) #0.135158
uniq_per_kings <- unique_percent (kings, kings_mut, la_mut, cook_mut) #0.092641
uniq_per_cook <- unique_percent (cook, cook_mut, la_mut, kings_mut) #0.06873

## 4. Evolutionary analysis
source('evolut_analysis.R')
lm_loess_la <- lm_loess("Los Angeles", la)
lm_loess_gl <- lm_loess("Global", global)
lm_loess_k <- lm_loess("Kings", kings)
lm_loess_c <- lm_loess("Cook", cook)

# figure 5
ggarrange(lm_loess_la, lm_loess_gl, lm_loess_k, lm_loess_c,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

# table 3
seg_la <- seg_results(la)

lm1 <- lm(totalSubstitutions ~ day, la)
bic1 <- selgmented(lm1, type = "bic")$selection.psi # BIC

lm2 <- lm(totalSubstitutions ~ day, global) # 2 breakpoints
bic2 <- selgmented(lm2, type = "bic")$selection.psi # BIC

lm3 <- lm(totalSubstitutions ~ day, kings)
bic3 <- selgmented(lm3, type = "bic")$selection.psi # BIC
lm4 <- lm(totalSubstitutions ~ day, cook)
bic4 <- selgmented(lm4, type = "bic")$selection.psi # BIC

co_list <- list(la, kings, cook)
co_seg_result <- lapply(co_list, seg_results)

result1 <- c("Los Angeles", bic1, co_seg_result[[1]])
result3 <- c("Kings", bic3, co_seg_result[[2]])
result4 <- c("Cook", bic4, co_seg_result[[3]])
row_name <- c("Region", "BIC-0", "BIC-1", "BIC-2", "P-value of Score Test", 
              "Breakpoint", "Breakpoint 95% CI.l", "Breakpoint 95% CI.u", 
              "Slope 1", "Slope 1 95% CI.l", "Slope 1 95% CI.u", 
              "Slope 2", "Slope 2 95% CI.l", "Slope 2 95% CI.u")
all_seg_results <- as_tibble(cbind(row_name, result1, result3, result4))

row_name2 <- c("Region", "BIC-0", "BIC-1", "BIC-2", "P-value of Score Test", 
               "Breakpoint 1", "Breakpoint 95% CI.l", "Breakpoint 95% CI.u", 
               "Breakpoint 2", "Breakpoint 2 95% CI.l", "Breakpoint 2 95% CI.u", 
               "Slope 1", "Slope 1 95% CI.l", "Slope 1 95% CI.u", 
               "Slope 2", "Slope 2 95% CI.l", "Slope 2 95% CI.u",
               "Slope 3", "Slope 3 95% CI.l", "Slope 3 95% CI.u")
co_seg_result_2 <- seg_results_2(global)
result2 <- c("Global", bic2, co_seg_result_2)
all_seg_results_2 <- as_tibble(cbind(row_name2, result2))

# figure 6
par(mfrow = c(2, 2))
seg_results_fig(la)
title("A Los Angeles", adj = 0)
seg_results_fig_2(global) # 2 breakpoints
title("B Global", adj = 0)
seg_results_fig(kings)
title("C Kings", adj = 0)
seg_results_fig(cook)
title("D Cook", adj = 0)

### Phylogenetic analysis
source('phylo_analysis.R')
## 1. Phylogenetic tree analysis
clade_la <- clade_dist("Los Angeles", la)
clade_gl <- clade_dist("Global", global)
clade_k <- clade_dist("Kings", kings)
clade_c <- clade_dist("Cook", cook)
#- figure 9 arrange on one page
ggarrange(clade_la, clade_gl, clade_k, clade_c,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

## 2. Clade temporal shift analysis
clade_change_num_la <- clade_change_num("Los Angeles", la)
clade_change_per_la <- clade_change_per("Los Angeles", la)
clade_change_num_gl <- clade_change_num("Global", global)
clade_change_per_gl <- clade_change_per("Global", global)
clade_change_num_k <- clade_change_num("Kings", kings)
clade_change_per_k <- clade_change_per("Kings", kings)
clade_change_num_c <- clade_change_num("Cook", cook)
clade_change_per_c <- clade_change_per("Cook", cook)
#- combine 2 plots
clade_shift_la <- clade_change_num_la + 
  annotation_custom(ggplotGrob(clade_change_per_la), xmin = "2020-01", 
                    ymin = 190, xmax = "2020-10")
clade_shift_gl <- clade_change_num_gl + 
  annotation_custom(ggplotGrob(clade_change_per_gl), xmin = "2019-12", 
                    ymin = 350, xmax = "2021-01")
clade_shift_kings <- clade_change_num_k + 
  annotation_custom(ggplotGrob(clade_change_per_k), xmin = "2020-05", 
                    ymin = 70, xmax = "2020-11")
clade_shift_cook <- clade_change_num_c + 
  annotation_custom(ggplotGrob(clade_change_per_c), xmin = "2020-03", 
                    ymin = 110, xmax = "2020-12")

#- arrange on one page
ggarrange(clade_shift_la, clade_shift_gl, clade_shift_kings, clade_shift_cook,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

#- percent
gl_21B <- global %>% filter(clade == "21B (Kappa)") %>% 
  group_by(region_exposure) %>% summarise(count = n())

la_21B <- la %>% filter(clade == "21B (Kappa)")

