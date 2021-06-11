genomic_char_analysis <- function(region, region_mut){
  region_mut_pos_count <- region_mut %>% group_by(position) %>% 
    summarise(count = n(), freq = n()/nrow(region))
  region_mut_pos_ununi <- filter(region_mut_pos_count, count > 1)
  
  region_sub_count <- region_mut %>% group_by(substitutions) %>% 
    summarise(count = n(), freq = n()/nrow(region))
  
  genomic_char <- tibble(point_mut = nrow(region_mut))
  genomic_char <- genomic_char %>% mutate(
    distinct_point_mut = nrow(distinct(region_mut, substitutions)),
    private_mut = nrow(filter(region_sub_count, count == 1)),
    private_mut_p = round((private_mut / distinct_point_mut) * 100, 2),
    no_mut = nrow(filter(region, totalSubstitutions == 0)),
    no_mut_p = round((no_mut / nrow(region)) * 100, 2),
    mut = nrow(filter(region, totalSubstitutions > 0)),
    mut_p = round((mut / nrow(region)) * 100, 2),
    min_mut = min(region$totalSubstitutions),
    max_mut = max(region$totalSubstitutions),
    mean_mut = round(mean(region$totalSubstitutions), 2),
    sd_mut = round(sd(region$totalSubstitutions), 2),
    position_mut = nrow(distinct(region_mut, position)),
    mean_pos_mut = round(mean(region_mut_pos_count$count), 2),
    sd_pos_mut = round(sd(region_mut_pos_count$count), 2),
    pos_mut_ununi = nrow(region_mut_pos_ununi),
    pos_mut_ununi_p = round((pos_mut_ununi/position_mut) * 100, 2),
    mean_pos_mut_ununi = round(mean(region_mut_pos_ununi$count), 2),
    sd_pos_mut_ununi = round(sd(region_mut_pos_ununi$count), 2)
  )
  return(genomic_char)
}