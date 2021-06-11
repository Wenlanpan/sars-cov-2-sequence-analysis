# sars-cov-2-sequence-analysis
 Master project: Analysis of Genomic Characteristics and Phylogeny of SARS-CoV-2 Introduced in Los Angeles
 
 Wenlan Pan, Department of Biostatistics, UCLA
 
## data
1. `GISAID`: Data downloaded from GISAID
2. `Nextclade`: Results from Nextclade
3. Save combined datasets using `dataset_pre.R`

## Analysis: analysis.R
1. Call different functions for analysis

## Data preparation: dataset_pre.R
1. Combined datasets of `data/GISAID` and `data/Nextclade`
2. Save in `data`

## Genomic characteristics analysis: genomic_char_analysis.R
1. `genomic_char_analysis`: a function producing a dataframe containing basic genomic characteristics 

### Top mutation analysis: top_mut_analysis.R
1. `top_mut_analysis`: a function producing a dataframe containing nucleotide mutation, count, frequency, date of first and last sample sorted by frequency from highest to lowest
2. `top_mut_analysis_fig`: a function producing a figure of the frequency of point mutations across the genomes in a region (unknown mutation type, label position)
3. `top_mut_analysis_fig_nt`: a function producing a figure of the frequency of point mutations across the genomes in a region (unknown mutation type, label nucleotide mutation)
4. `top_mut_analysis_fig_la`: a function producing a figure of the frequency of point mutations across the genomes in LA (known mutation types, label nucleotide mutation)

### Co-muation analysis: co-mut_analysis.R
1. `co_mut_analysis`: a function producing a scaled concurrence heat map of the top most common mutations in a region

### Comparative analysis: comp_analysis.R
1. `comp_analysis`: a function producing a lollipop figure of top mutations in a region compare with US & global datasets
2. `mut_freq_ca`: a function producing a dataframe containing nucleotide mutation, count, frequnecy
3. `comp_analysis_tab`: a function producing a table comparising the top frequent mutations in a region with the Kings, Cook, US, and Global datasets
4. `venn_diag`: a function producing a venn diagram of th mutually common and exclusive substitutions among genomes from three regions
5. `unique_percent`: a function calculating the number and percentage of mutation from only a region among three regions

### Evolutionary analysis: evolut_analysis.R
1. `lm_loess`: a function producing a figure of linear and loess regression between time of emergence and number of mutations per sample in a region
2. `seg_results`: a function producing the results of the segemented regression with one breakpoint between time of emergence and number of mutations per sample in a region (p-value of the score test for number of breakpoint, breakpoint est and 95\%CI, slope est and 95\% CI) - BIC results can be added in `analysis.R`
3. `seg_results_2`: a function producing the results of the segemented regression with two breakpoints between time of emergence and number of mutations per sample in a region (p-value of the score test for number of breakpoint, breakpoint est and 95\%CI, slope est and 95\% CI) - BIC results can be added in `analysis.R`
4. `seg_results_fig`: a function producing a figure of the segmented regression with one breakpoint between time of emergence and number of mutations per sample in a region (label slope and breakpoint est)
5. `seg_results_fig_2`: a function producing a figure of the segmented regression with two breakpoints between time of emergence and number of mutations per sample in a region (label slope and breakpoint est)

## Phylogenetic analysis: phylo_analysis.R
1. `clade_dist`: a function producing clade distribution in a region
2. `clade_change_num`: a function producing temporal change (number) in the predominant clades in a region
3. `clade_change_per`: a function producing temporal change (percentage) in the predominant clades in a region
