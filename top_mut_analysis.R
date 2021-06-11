# table 1
top_mut_analysis <- function(region, region_mut){
  df <- region_mut %>% group_by(substitutions) %>% 
    summarise(num = n(), first_date = min(date), last_date = max(date)) %>% 
    mutate(freq = round((num/nrow(region)) * 100, 2)) %>% arrange(desc(freq))
  return(df)
}

# Lollipop: mutations that are prevalent in a region (unknown mutation type)
top_mut_analysis_fig <- function(region, region_mut, top = 10){
  df <- region_mut %>% group_by(position) %>% summarise(num = n()) %>% 
    mutate(freq = num/nrow(region))
  mut_top <- df %>% arrange(desc(freq)) %>% top_n(top, freq)
  
  p <- ggplot(df, aes(x = position, y = freq)) +
    geom_segment(aes(x = position, xend = position, y = 0, yend = freq),
                 color = ifelse(df$position %in% mut_top$position, 
                                alpha("#4E84C4", 0.3), alpha("#999999", 0.3)),
                 size = ifelse(df$position %in% mut_top$position, 1.3, 0.7)) +
    geom_point(color = ifelse(df$position %in% mut_top$position, 
                              alpha("#4E84C4", 0.3), "#999999"),
               size = ifelse(df$position %in% mut_top$position, 5, 2),
               alpha = 0.7) + theme_classic() + 
    labs(x = "Base pair position", y = "Frequency") +
    geom_text(data = df %>% arrange(desc(freq)) %>% top_n(top, freq), 
              aes(label = position), vjust = -1, check_overlap = T) 
  return(p)
}

# Lollipop: mutations that are prevalent in a region 
# unknown mutation type, label nucleotide mutation
# error: frequence is different between position and substitution 
# i.e. 22917 rank is different from T22917G
top_mut_analysis_fig_nt <- function(region, region_mut, top = 10){
  df <- region_mut %>% group_by(position) %>% summarise(num = n()) %>% 
    mutate(freq = num/nrow(region))
  df2 <- region_mut %>% group_by(substitutions) %>% summarise(num = n()) %>% 
    mutate(freq = num/nrow(region)) %>% arrange(desc(freq)) %>% top_n(top, freq)
  mut_top <- df %>% arrange(desc(freq)) %>% top_n(top, freq)
  
  p <- ggplot(df, aes(x = position, y = freq)) +
    geom_segment(aes(x = position, xend = position, y = 0, yend = freq),
                 color = ifelse(df$position %in% mut_top$position, 
                                alpha("#4E84C4", 0.3), alpha("#999999", 0.3)),
                 size = ifelse(df$position %in% mut_top$position, 1.3, 0.7)) +
    geom_point(color = ifelse(df$position %in% mut_top$position, 
                              alpha("#4E84C4", 0.3), "#999999"),
               size = ifelse(df$position %in% mut_top$position, 5, 2),
               alpha = 0.7) + theme_classic() + 
    labs(x = "Base pair position", y = "Frequency") +
    geom_text(data = df %>% arrange(desc(freq)) %>% top_n(top, freq), 
              aes(label = df2$substitutions), vjust = -1, check_overlap = T) 
  return(p)
}

# figure 1 lollipop: mutations that are prevalent in LA (known mutation type)
# change label_NT
top_mut_analysis_fig_la <- function(region, region_mut){
  df <- region_mut %>% group_by(position) %>% summarise(num = n()) %>% 
    mutate(freq = num/nrow(region))
  mut_top <- df %>% arrange(desc(freq)) %>% top_n(10, freq)
  # annotate NT substitution
  label_NT <- c("A23403G", "C14408T", "C3037T", "C241T", "G25563T", "C1059T",
                "C28887T", "T22917G", "G17014T", "C26681T")
  
  p <- ggplot(df, aes(x = position, y = freq)) +
    geom_segment(aes(x = position, xend = position, y = 0, yend = freq),
                 color = ifelse(df$position == "241", "#D16103", 
                                ifelse(df$position %in% c("3037", "26681"), "#52854C", 
                                       ifelse(df$position %in% mut_top$position,
                                              "#4E84C4", "#999999"))),
                 size = ifelse(df$position %in% mut_top$position, 1.3, 0.7)) +
    geom_point(color = ifelse(df$position == "241", "#D16103", 
                              ifelse(df$position %in% c("3037", "26681"), "#52854C", 
                                     ifelse(df$position %in% mut_top$position, 
                                            "#4E84C4", "#999999"))),
               size = ifelse(df$position %in% mut_top$position, 5, 2)) + 
    theme_classic() + labs(x = "Base pair position", y = "Frequency") +
    geom_text(data = df %>% arrange(desc(freq)) %>% top_n(10, freq), 
              aes(label = label_NT), vjust = -1, check_overlap = T)
  return(p)
}
