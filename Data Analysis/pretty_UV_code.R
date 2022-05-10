#import libraries
library(tidyverse); theme_set(theme_light())
library(googlesheets4)
library(viridis)

#import data from google sheets
url <- "[redacted sheet URL]"
streaks <- read_sheet(ss = url, sheet = "Streaks")
srl_dils <- read_sheet(ss = url, sheet = "Serial_Dilutions")
nming_dils <- read_sheet(ss = url, sheet = "Normalizing_Dilutions")
nmed_dils <- read_sheet(ss = url, sheet = "Normalized_Dilutions")

#normalize the starting CFU/mL for each strain
normalizing <-nming_dils %>%
  group_by(Strain)%>%
  summarise(mean = mean(CFU_mL)) %>%
  arrange(mean) %>%
  mutate(ul_culture = mean[1]/mean*1000,
         ul_PBS = 1000-ul_culture)
normalizing

#calculation and visualization of percent survivals
#this was performed for normalized and unnormalized, only the former is included
n_sixtys_uv <- nmed_dils %>%
  filter(Time == 60) %>%
  group_by(Strain, Set) %>%
  summarise(CFU_60 = mean(CFU_mL))

n_zeros_uv <- nmed_dils %>%
  filter(Time == 0) %>%
  group_by(Strain, Set) %>%
  summarise(CFU_0 = mean(CFU_mL))

n_sz_uv <- left_join(n_zeros_uv, n_sixtys_uv)

n_sz_uv <- n_sz_uv %>%
  mutate(p_survival = CFU_60/CFU_0)

write.csv(n_sz_uv,here::here("data", "n_sz_uv.csv"), row.names = FALSE)

ggplot(n_sz_uv, aes(x = Strain, y = p_survival*100, color = Strain)) +
  geom_boxplot() +
  labs(title = "Percent Survival of Normalized Strains in UV irradiation",
       y = "Percent Survival")+
  scale_color_viridis(discrete=TRUE)  +
  theme(text=element_text(family="serif"))

#t-tests
n_dh5a <- n_sz_uv %>%
  filter(Strain == 	"DH5A")
n_drad <- n_sz_uv %>%
  filter(Strain == 	"DRAD")
n_f5 <- n_sz_uv %>%
  filter(Strain == 	"1F5")
n_fh2 <- n_sz_uv %>%
  filter(Strain == 	"1FH2")
n_cm1 <- n_sz_uv %>%
  filter(Strain == 	"CM1")

t.test(x = n_sz_uv$CFU_0, y = n_sz_uv$CFU_60) #all strains
t.test(x = n_dh5a$CFU_0, y = n_dh5a$CFU_60) #dh5a
t.test(x = n_drad$CFU_0, y = n_drad$CFU_60) #drad
t.test(x = n_f5$CFU_0, y = n_f5$CFU_60) #f5
t.test(x = n_cm1$CFU_0, y = n_cm1$CFU_60) #cm1
t.test(x = n_fh2$CFU_0, y = n_fh2$CFU_60) #fh2