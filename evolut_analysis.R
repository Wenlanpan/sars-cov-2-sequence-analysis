# figure 5 linear and loess regression between time of emergence and number of mutations 
# per sample in a region
lm_loess <- function(name, region){
  p <- ggplot(region, aes(x = day, y = totalSubstitutions)) +
    geom_point(col = c("#999999")) +
    geom_smooth(method = "loess", col = c("#4E84C4")) + 
    stat_smooth(method = "lm", formula = y ~ x, col = c("#D16103")) + 
    labs(title = name,
         x = paste0("Days from ", min(region$date)), y = "Number of Mutations per Sample") +
    theme_classic() 
  return(p)
}

# table 3 results of the segemented regression between time of emergence and number of mutations 
# per sample in a region
seg_results <- function(region){
  lm <- lm(totalSubstitutions ~ day, region) # linear regression
  segm <- segmented(lm) # segmented regression
  pscore <- pscore.test(segm)$p.value # <2.2e-16 reject -> one breakpoint better
  breakpoint <- t(confint(segm)) # 297.002, 283.366, 310.637 delta CI
  slope <- slope(segm)$day # slope
  result <- c(pscore, breakpoint, slope[1, 1], slope[1, 4:5], 
              slope[2, 1], slope[2, 4:5])
  return(result)
}

seg_results_2 <- function(region){
  lm <- lm(totalSubstitutions ~ day, region) # linear regression
  segm <- segmented(lm, npsi = 2) # segmented regression
  pscore <- pscore.test(segm, more.break = TRUE)$p.value
  # an additional breakpoint(s) for the variable
  breakpoint <- confint(segm)
  slope <- slope(segm)$day # slope
  result <- c(pscore, breakpoint[1, 1:3], breakpoint[2, 1:3], slope[1, 1], 
              slope[1, 4:5], slope[2, 1], slope[2, 4:5], slope[3, 1], slope[3, 4:5])
  return(result)
}

# figure 6 the segmented regression between time of emergence and number of mutations 
# per sample in a region
seg_results_fig <- function(region){
  lm <- lm(totalSubstitutions ~ day, region) # linear regression
  segm <- segmented(lm) # segmented regression
  breakpoint <- confint(segm)
  slope <- slope(segm)$day
  plot(region$day, region$totalSubstitutions, 
       xlab = paste0("Days from ", min(region$date)),
       ylab = "Number of Mutations per Sample", pch = 20, col = c("#999999"))
  plot(segm, add = TRUE, link = FALSE, lwd = 3, col = c("#D16103","#52854C"), 
       conf.level = 0.95, shade = TRUE)
  lines(segm, col = c("#4E84C4"), pch = 19, bottom = FALSE, lwd = 2) 
  # CI for breakpoint
  points(segm, col = c("#4E84C4"), link = FALSE, lwd = 3)
  legend(0, 45, c(paste0("Slope 1: ", round(slope[1, 1], 3)), 
                  paste0("Slope 2: ", round(slope[2, 1], 3)), 
                  paste0("Breakpoint: ", round(breakpoint[1, 1], 0))),
         col = c("#D16103","#52854C", "#4E84C4"), y.intersp = 0.5,
         lty = c(1, 1, NA), lwd = c(3, 3, 3), pch = c(NA, NA, 1), bty = "n")
}

seg_results_fig_2 <- function(region){
  lm <- lm(totalSubstitutions ~ day, region) # linear regression
  segm <- segmented(lm, npsi = 2) # segmented regression
  breakpoint <- confint(segm)
  slope <- slope(segm)$day
  plot(region$day, region$totalSubstitutions, 
       xlab = paste0("Days from ", min(region$date)),
       ylab = "Number of Mutations per Sample", pch = 20, col = c("#999999"))
  plot(segm, add = TRUE, link = FALSE, lwd = 3, col = c("#D16103","#52854C",
                                                        "#f88421"), 
       conf.level = 0.95, shade = TRUE)
  lines(segm, col = c("#4E84C4", "#006b7b"), pch = 19, bottom = FALSE, lwd = 2) 
  # CI for breakpoint
  points(segm, col = c("#4E84C4", "#006b7b"), link = FALSE, lwd = 3)
  legend(0, 60, c(paste0("Slope 1: ", round(slope[1, 1], 3)), 
                  paste0("Slope 2: ", round(slope[2, 1], 3)), 
                  paste0("Slope 3: ", round(slope[3, 1], 3)),
                  paste0("Breakpoint 1: ", round(breakpoint[1, 1], 0)),
                  paste0("Breakpoint 1: ", round(breakpoint[2, 1], 0))),
         col = c("#D16103","#52854C", "#f88421", "#4E84C4", "#006b7b"), 
         y.intersp = 0.5, lty = c(1, 1, 1, NA, NA), lwd = c(3, 3, 3, 3, 3), 
         pch = c(NA, NA, NA, 1, 1), bty = "n")
}
