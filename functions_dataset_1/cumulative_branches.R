cumulative_branches <- function(total) {
  
  
  
  # bin formation and collapsing events
  total <- total %>%
    mutate(
      "cBins" = cut(Collapse, breaks = seq(0, 24, 2), include.lowest = TRUE),
      "Bins" = cut(Formation, breaks = seq(0, 24, 2), include.lowest = TRUE)
    )
  
  
  
  # count newly formed per bin
  temp <- total %>%
    #filter(Collapse < 24) %>%
    group_by(Genotype, Location, Branchtype, CellSum, Bins) %>%
    summarise("Formed" = n()) %>%
    ungroup()
  
  # count collapsed per bin
  temp1 <- total %>%
    filter(Collapse < 24) %>%
    group_by(Genotype, Location, Branchtype, CellSum, cBins) %>%
    summarise("Collapsed" = n()) %>%
    ungroup()
  
  
  # merge to a binned count 
  binned_branches <- full_join(temp, temp1,
                               by = c(
                                 "Genotype" = "Genotype",
                                 "Location" = "Location",
                                 "Branchtype" = "Branchtype",
                                 "Bins" = "cBins",               
                                 "CellSum" = "CellSum"
                               )
  ) %>% 
    # add groups for non-existent bins
    complete(Bins, nesting(Genotype, Location, Branchtype))
  
  
  
  Cumulative <- binned_branches %>%
    # in this case, NAs mean 0 collapses or formations in this bin
    mutate(Collapsed = replace_na(Collapsed, 0),
           Formed = replace_na(Formed, 0)) %>%
    
    # sorting before cumulative sum!
    arrange(Genotype, Location, Branchtype, Bins) 
  
  
  rm(temp, temp1)
  
  
  Cumulative <- Cumulative %>%
    group_by(Genotype, Location, Branchtype) %>%
    # accumulate 
    mutate("flux" = Formed - Collapsed,
           "accumulated_total" = cumsum(flux))
  
  
  Cumulative <- Cumulative %>%
    group_by(Genotype) %>%
    mutate("Accumulated_normalized" = accumulated_total / CellSum)
  
  
  return(Cumulative)
}
