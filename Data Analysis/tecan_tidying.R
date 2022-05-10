#import libraries
library(dplyr)
library(tidyr)
library(janitor)
library(readxl)

# formatting raw file
tecan_raw <- read_excel("~/Thesis/RRR_Bacteria/Data/Li_24hr_growth.xlsx", skip = 34, col_names = FALSE)

tecan_raw <- tecan_raw %>%
  slice(1:58) %>%
  t() 

tecan_data <- row_to_names(tecan_raw, row_number = 1) 
#outputs as a matrix rather than a df. refused to run select or rename or filters

tecan_data <- as.data.frame(tecan_data) #this ensures that we're working with a df

tecan_data <- tecan_data %>%
  rename(cyclen = "Cycle Nr.",
         time = "Time [s]",
         temp = "Temp. [°C]") %>%
  mutate(cycle = 1:n()) %>%
  select(cycle, everything(), -cyclen) 

# tidying
tecan_tidy <- tecan_data %>% 
  pivot_longer(               #changes df from wider to longer
    cols = 4:58,              #columns with data
    values_to = "a600",       #new column name for all data
    names_to = "well"         #these were the original column names for data
  ) %>%
  mutate(strain = case_when(
    grepl("1", well) ~ "CONTROL", #checks for string in value. All wells with "1" are control LB
    grepl("2", well) ~ "DH5A",
    grepl("3", well) ~ "DRAD",
    grepl("4", well) ~ "CM1",
    grepl("5", well) ~ "1FH2",
    grepl("6", well) ~ "1F5",
    grepl("7", well) ~ "1F1"
  )) %>%                      
  #at this time, all values are currently stored as characters just because that's how they were imported. 
  transform(time = as.numeric(time),
            temp = as.numeric(temp),
            a600 = as.numeric(a600))

# export df to csv file
write.csv(tecan_tidy,"~/Thesis/RRR_Bacteria/Data/tecan_tidy.csv", row.names = FALSE)


