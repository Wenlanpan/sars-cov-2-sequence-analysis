# figure 9 clade distribution in a region
clade_dist <- function(name, region){
  region_clade <- region %>% group_by(clade) %>% summarise(count=n()) %>% 
    mutate(freq = count/sum(count))
  region_label <- paste(region_clade$clade, "(", region_clade$count, ",",
                 round(region_clade$freq * 100, 1), "%)", sep = '')
  clade_num <- nrow(region_clade)
  if (clade_num  > 11){
    p <- ggplot(region_clade, aes(1, count, fill = clade)) +
      geom_bar(stat = "identity", position = "stack", width = 1) +
      coord_polar(theta = "y") + labs(title = name, fill = "Lineages") +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(hjust = 0.5)) + 
      scale_fill_manual(labels = region_label, 
                        values = colorRampPalette(brewer.pal(11, "Spectral"))(clade_num)) +
      theme_void() + easy_center_title()
  } else {
    p <- ggplot(region_clade, aes(1, count, fill = clade)) +
      geom_bar(stat = "identity", position = "stack", width = 1) +
      coord_polar(theta = "y") + labs(title = name, fill = "Lineages") +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(hjust = 0.5)) + 
      scale_fill_brewer(labels = region_label, palette = "Spectral") +
      theme_void() + easy_center_title()
  }
  return(p)
}

# figure 10 1 temporal change (number) in the predominant clades in a region
clade_change_num <- function(name, region){
  df <- region %>% group_by(month, clade) %>% summarize(count = n())
  region_clade <- region %>% group_by(clade) %>% summarise(count=n()) %>% 
    mutate(freq = count/sum(count))
  clade_num <- nrow(region_clade)
  if (clade_num > 11){
    p <- ggplot(df, aes(x = month, y = count, fill = clade)) +
      geom_bar(stat = "identity") + theme_classic() +
      labs(x = "Date (Month)", y = "Number of Genomes by Lineages ",
           title = name, fill = "Lineage") +
      theme(axis.text.x=element_text(angle=30, vjust=.8, hjust=.8)) +
      scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(clade_num ))
  } else {
    p <- ggplot(df, aes(x = month, y = count, fill = clade)) +
      geom_bar(stat = "identity") + theme_classic() +
      labs(x = "Date (Month)", y = "Number of Genomes by Lineages ",
           title = name, fill = "Lineage") +
      theme(axis.text.x=element_text(angle=30, vjust=.8, hjust=.8)) +
      scale_fill_brewer(palette = "Spectral")
  }
  return(p)
}

# figure 10 2 temporal change (percentage) in the predominant clades in a region
clade_change_per <- function(name, region){
  value_region <- c(rep(1, nrow(region)))
  df <- region %>% select(clade, month) %>% mutate(value = value_region)
  region_clade <- region %>% group_by(clade) %>% summarise(count=n()) %>% 
    mutate(freq = count/sum(count))
  clade_num <- nrow(region_clade)
  if (clade_num > 11){
    p <- ggplot(df, aes(fill = clade, y = value, x = month)) +
      geom_bar(position = "fill", stat = "identity") + theme_classic() +
      labs(y = "Percentage") +
      theme(axis.text.x= element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            legend.position = "none") +
      scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(clade_num))
  } else {
    p <- ggplot(df, aes(fill = clade, y = value, x = month)) +
      geom_bar(position = "fill", stat = "identity") + theme_classic() +
      labs(y = "Percentage") +
      theme(axis.text.x= element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            legend.position = "none") +
      scale_fill_brewer(palette = "Spectral")
  }
  return(p)
}