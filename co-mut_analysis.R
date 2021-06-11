# Scaled concurrence heat map of the ten most common mutations in a region
co_mut_analysis <- function(region, region_mut, top = 10){
  mut_top <- region_mut %>% group_by(substitutions) %>% summarise(num = n()) %>%
    mutate(freq = num/nrow(region)) %>% arrange(desc(freq)) %>% top_n(top, freq)
  df <- filter(region_mut, substitutions %in% mut_top$substitutions) 
  mat <- crossprod(table(df$strain, df$substitutions))
  p <- pheatmap(mat, display_numbers = TRUE, scale = "row", cluster_cols = FALSE, 
                border_color = NA) # scale by row
  return(p)
}