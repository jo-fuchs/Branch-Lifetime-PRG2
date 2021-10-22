clean_data <- function() {

# Data cleaning, creating secondary measures

branches <- read_csv(file.path("data", "dataset_2", "Individual_branches.csv"))
cells <-  read_csv2(file.path("data", "dataset_2", "Cells.csv"))




## Cells ----

## Branches ----
total <- branches %>% 
  left_join(cells, 
            by = c("Movie" = "Movie_ID",
                   "Experiment" = "Imaging_ID"))

# Total cells per treatment
total <- total %>% 
  group_by(Treatment) %>% 
  mutate("CellSum" = sum(unique(Neurons))) %>% 
  ungroup()


# remove BDNF
total <- total %>% filter(Treatment %in% c("Control", "Netrin", "FGF"))



# remove unmeasured
# total <- total %>% filter(!is.na(Collapse))

# remove DIV1 data from experiment C, cut down to 30 hours (181 frames) 
 
 total <- total %>%
   mutate(Formation = if_else(Culture %in% c("C"), Formation - 145, Formation),
          Collapse = if_else(Culture %in% c("C"), Collapse - 145, Collapse),
          Collapse = ifelse(Collapse > 181, 181, Collapse)) %>%
   filter(Formation > 0) %>% 
   filter(Formation < 181)

 
# convert 10 min frames to hours (Frame / 60 * 10 min)
total <- total %>% 
  mutate(
    Formation = (Formation - 1) / 6,
    Collapse  = (Collapse - 1) / 6)



# Lifetime
total <- total %>% 
  mutate(Lifetime = Collapse - Formation)




# Indicator for censoring
total <- total %>% mutate(CompleteData = as.numeric(Collapse != max(Collapse)))

# Inverse probability weighting indicator
# (branches forming close to censoring timepoint have a higher chance for being censored)
total <- total %>% mutate(Prob = 1 / (max(Collapse) - Formation))


# Set factor levels
total <- total %>% 
  mutate(
    Branchtype = factor(Branchtype, 
                        levels = c("Filopodium", "Mixed", "Lamellipodium", "Splitting")),
    Location   = factor(Location,
                        levels = c("Axon", "Neurite", "Unclear")),
    Treatment  = factor(Treatment,
                        levels = c("Control", "Netrin", "FGF")),
    Culture    = factor(Culture,
                        levels = c("A", "B", "C"))
    
  ) #%>% filter(Culture != "A")



return(total)
}
