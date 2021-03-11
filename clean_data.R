clean_data <- function() {

# Data cleaning, creating secondary measures

branches <- read_csv2(file.path("data", "Individual_branches.csv"))
cells <-  read_csv2(file.path("data", "Cells.csv"))




## Cells ----
# Calculate total cell number
cells <- cells %>% 
  mutate("Neurons" = Cells_polarized + Cells_unpolarized + Cells_Other)


## Branches ----
total <- branches %>% 
  left_join(cells, 
            by = c("Imaging_ID", 
                   "Movie" = "Movie_ID",
                   "Original_Name",
                   "Culture",
                   "Genotype"))

# Total cells per treatment
total <- total %>% 
  group_by(Genotype) %>% 
  mutate("CellSum" = sum(unique(Neurons))) %>% 
  ungroup()


# remove DIV1 data from experiments E & F (imaging started earlier than in A-D)
total <- total %>% 
  mutate(Formation = if_else(Culture %in% c("E", "F"), Formation - 144, Formation),
         Collapse = if_else(Culture %in% c("E", "F"), Collapse - 144, Collapse)) %>% 
  filter(Formation > 0)


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
    Branchtype = fct_recode(Color, 
                                 "Filopodium" = "red", 
                                 "Mixed" = "yellow",
                                 "Lamellipodium" = "green", 
                                 "Splitting" = "blue"),
    Branchtype = factor(Branchtype, 
                        levels = c("Filopodium", "Mixed", "Lamellipodium", "Splitting")),
    Location   = factor(Location,
                        levels = c("Axon", "Neurite", "Unclear")),
    Genotype   = fct_recode(Genotype,
                        "Plppr3 +/+" = "WT", 
                        "Plppr3 -/-" = "KO"),
    
    Genotype   = factor(Genotype,
                        levels = c("Plppr3 +/+", "Plppr3 -/-")),
    Culture    = factor(Culture,
                        levels = c("A", "B", "C", "D", "E", "F"))
    
  )



return(total)
}
