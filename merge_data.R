merge_data <- function() {

  library(tidyverse)

## Individual branchpoints ----
# list all folders
dirs <- list.dirs(path = "data", full.names = TRUE, recursive = F)

# load individual branchpoints from each folder
all_cells <- tibble()

for (folder in dirs) {
  
  subdir <- file.path(folder, "lifetimes")
  
  dirlist <- list.files(path = subdir , pattern = ".csv$")
  individual<- vroom::vroom(file.path(subdir, dirlist), id = "Movie_ID")
  
  all_cells <- rbind(all_cells, individual)
  rm(individual)
}


## clean individual branch points
all_cells_clean <- all_cells %>%
  extract(Movie_ID, into = c("Imaging_ID", "Movie"), regex = "data/(.+)/lifetimes/(.+).csv") 




# Cell counts and Randomizations ----
randomization <- tibble()
cells <- tibble()

for (folder in dirs) {
  rand_ind <- read_csv2(file.path(folder, "Randomization.csv"))
  cells_ind <- read_csv2(file.path(folder, "Cells.csv"))
  
  rand_ind$Imaging_ID <- folder
  cells_ind$Imaging_ID <- folder
   
  randomization <- rbind(randomization, rand_ind) 
  cells <- rbind(cells, cells_ind) 
  
  rm(rand_ind, cells_ind)
}

# clean data
randomization <- randomization %>%
  extract(Imaging_ID, into = "Imaging_ID", regex = "data/(.+)")

cells <- cells %>% 
  extract(Imaging_ID, into = "Imaging_ID", regex = "data/(.+)")


## Unblinding ----
cells_final <- cells %>% 
  left_join(randomization, by = c("Imaging_ID", 
                                  "Movie_ID" = "Random_Name"))

final <- all_cells_clean %>% 
  left_join(randomization, by = c("Imaging_ID", 
                                  "Movie" = "Random_Name")) 



## save results
write_csv2(final, file.path("data", "Individual_branches.csv"))
write_csv2(cells_final, file.path("data", "Cells.csv"))

}
