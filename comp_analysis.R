# figure 3 Lollipop: Top mutations in a region compare with US & global datasets
comp_analysis <- function(name, region, region_mut, global, global_mut, top = 10){
  region_df <- region_mut %>% group_by(substitutions) %>% 
    summarise(num = n()) %>% mutate(freq = num/nrow(region)) %>% 
    arrange(desc(freq)) %>% top_n(top, freq) %>% mutate(region = name)
  global_df <- global_mut %>% group_by(substitutions) %>% 
    summarise(num = n()) %>% mutate(freq = num/nrow(global)) %>% 
    arrange(desc(freq)) %>% top_n(top, freq) %>% mutate(region = "Global")
  us_num <- filter(global, country_exposure == "USA")
  us_df <- global_mut %>% filter(country_exposure == "USA") %>% 
    group_by(substitutions) %>% summarise(num = n()) %>% 
    mutate(freq = num/nrow(us_num)) %>%
    arrange(desc(freq)) %>% top_n(top, freq) %>% mutate(region = "United States")
  top_mutation <- rbind(global_df, us_df, region_df)
  top_mutation$region <- factor(top_mutation$region, 
                                levels = c("Global", "United States", name))
  
  p <- ggdotchart(top_mutation, x = "substitutions", y = "freq",
                  col = "region", add = "segments", rotate = TRUE,
                  dot.size = 8, label = (round(top_mutation$freq, 3) * 100),
                  palette = c("#D16103", "#52854C", "#4E84C4"),
                  font.label = list(size = 8, color = "white", vjust = 0.5),
                  ggtheme = theme_classic()) + labs(x = "Substitutions", 
                                                    y = "Frequency", 
                                                    color = "Region")
  return(p)
}

mut_freq_ca <- function(region, region_mut){
  df <- region_mut %>% group_by(substitutions) %>% summarise(num = n()) %>%
    mutate(freq = round((num/nrow(region)) * 100, 2))
  return(df)
}

# table 2 Top mutations in a region compare with Kings, Cook, US & global datasets
comp_analysis_tab <- function(name, region, region_mut, kings, kings_mut,
                              cook, cook_mut, global, global_mut, top = 10){
  region_df <- mut_freq_ca(region, region_mut) %>% arrange(desc(freq)) %>% 
    top_n(top, freq) %>% mutate(region = name)
  kings_df <- mut_freq_ca(kings, kings_mut) %>% 
    filter(substitutions %in% region_df$substitutions) %>% mutate(region = "Kings")
  cook_df <- mut_freq_ca(cook, cook_mut) %>% 
    filter(substitutions %in% region_df$substitutions) %>% mutate(region = "Cook")
  global_df <- mut_freq_ca(global, global_mut) %>%
    filter(substitutions %in% region_df$substitutions) %>% mutate(region = "Global")
  us_num <- global %>% filter(country_exposure == "USA")
  us_df <- global_mut %>% filter(country_exposure == "USA") %>% 
    group_by(substitutions) %>% summarise(num = n()) %>% 
    mutate(freq = round((num/nrow(us_num)) * 100, 2)) %>% 
    filter(substitutions %in% region_df$substitutions) %>% mutate(region = "US")
  dfs <- list(region_df, kings_df, cook_df, us_df, global_df)
  df <- dfs %>% reduce(left_join, by = c("substitutions"))
  return(df)
}

# figure 4 Venn diagram of three regions
venn_diag <- function(name1, region_mut1, name2, region_mut2, name3, region_mut3){
  region1_mut_uniq <- distinct(region_mut1, substitutions)
  region2_mut_uniq <- distinct(region_mut2, substitutions)
  region3_mut_uniq <- distinct(region_mut3, substitutions)
  df <- list(region2_mut_uniq$substitutions, region3_mut_uniq$substitutions,
             region1_mut_uniq$substitutions)
  names(df) <- c(name2, name3, name1)
  g <- ggvenn(df, show_percentage = T, stroke_color = "white",
         fill_color = c("#D16103","#52854C","#4E84C4"),
         set_name_color = c("#D16103","#52854C", "#4E84C4"))
  return(g)
}

# the percentage of mutations from only a region
unique_percent <- function(region1, region_mut1, region_mut2, region_mut3){
  region1_mut_uniq <- distinct(region_mut1, substitutions)
  region2_mut_uniq <- distinct(region_mut2, substitutions)
  region3_mut_uniq <- distinct(region_mut3, substitutions)
  only_region1 <- region1_mut_uniq %>% 
    anti_join(region2_mut_uniq, by = "substitutions") %>% 
    anti_join(region3_mut_uniq, by = "substitutions")
  region1_count <- region_mut1 %>% group_by(substitutions) %>% 
    summarise(num = n()) %>% right_join(only_region1, by = "substitutions")
  percent <- sum(region1_count$num)/sum(region1$totalSubstitutions)
  return(list(sum(region1_count$num), sum(region1$totalSubstitutions), percent))
}
