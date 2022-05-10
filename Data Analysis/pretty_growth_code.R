#import libraries
library(tidyverse); theme_set(theme_light())
library(viridis)
library(zoo)
library(moderndive)

#obtain relevant data and calculate average values
tecan_tidy_path <- here::here("data", "tecan_tidy.csv")
tecan_tidy <- read.csv(tecan_tidy_path, stringsAsFactors = FALSE)
tecan_avg <- tecan_tidy %>%
  group_by(strain, time)%>%
  summarise(time_min = time/60,
            time_hr = time_min/60,
            avg_A600 = mean(a600),
            temp = log10(avg_A600*1000)) %>%
  ungroup() %>%
  group_by(strain, time_hr) %>%
  summarize(log_A600 = mean(temp)) %>%
  select(strain, time_hr, log_A600) %>%
  filter(strain != "1F1")

write.csv(tecan_avg,here::here("data", "tecan_avg.csv"), row.names = FALSE)

#visualized data, this is Growth Curves graph.
ggplot(tecan_avg, aes(x = time_hr, y = log_A600)) +
  geom_point(aes(shape = strain, color = strain)) +
  scale_color_viridis(discrete=TRUE)  +
  labs(title = "Growth Curves from 24 hour Tecan Run in LB",
       x = "Time (hr)",
       y = "Log(Absorbance*1000) at 600nm",
       color = "Strain",
       shape = "Strain") +
  theme(text=element_text(family="serif"))
#in order for shape and color to be on the same legend, the manual settings must be the same in name and label. shape values are for specific shapes

#make individual dataframes to run through function
dh5a <- tecan_avg %>%
  filter(strain == "DH5A")
f5 <- tecan_avg %>%
  filter(strain == "1F5")
fh2 <- tecan_avg %>%
  filter(strain == "1FH2")
cm1 <- tecan_avg %>%
  filter(strain == "CM1")
drad <- tecan_avg %>%
  filter(strain == "DRAD")
strain_dfs <- list(dh5a, drad, f5, fh2, cm1)

#functions to optimize linear regressions
lin_eq <- function (df) {            # this function returns a slope
  d <- as.data.frame(df)
  lin_mod <- lm(log_A600~time_hr, d)    # m <- lm(y~x, as.data.frame(d)). outputs the linear model as list of values
  slope_temp <- coef(lin_mod)[2]
  return(slope_temp)
} 

lm_optimize <- function(df) {                                # this function calculates many slopes returns values with most common slope
  d <- as.data.frame(df)
  applied <- rollapply(d, 3, lin_eq, by.column = F)          # 3 is the reading frame of data. runs function many times. 
  clustering <- kmeans(applied, 3)                           # 3 is sets of matching. higher number results in fewer points per set
  cluster_positions <- which(clustering$cluster == match(max(clustering$centers), clustering$centers))+1   # looks for position of values
  RES <- d[cluster_positions,]                               # collects the best cluster set for best linear model. outputs df. 
  return(RES)
}

#run dataframes through function
# all the values put into one df for ease of graphing
tecan_sets <- data.frame()
i <- 1
for (i in 1:length(strain_dfs)) {
  d <- as.data.frame(strain_dfs[[i]])
  temp <- as.data.frame(lm_optimize(d))
  tecan_sets <- bind_rows(tecan_sets, temp)
  i <- i + 1
} 

# all the lm variables consolidated for easy comparison
lms_sets <- data.frame("CONTROL", 0, 0, 0)
names(lms_sets) <- c("strain", "intercept", "slope", "adj_r_sq")
i <- 1
for (i in 1:length(strain_dfs)) {
  d <- as.data.frame(strain_dfs[[i]])
  temp <- as.data.frame(lm_optimize(d))
  lin_mod <- lm(log_A600~time_hr, temp)    
  lms_sets <- rbind(lms_sets, c(strain = d[[1,1]],
                                intercept = coef(lin_mod)[1],
                                slope = coef(lin_mod)[2],
                                adj_r_sq = summary(lin_mod)$adj.r.squared))
  i <- i + 1
} 
lms_sets <- lms_sets %>%
  transform(intercept = as.numeric(intercept),
            slope = as.numeric(slope),
            adj_r_sq = as.numeric(adj_r_sq))

write.csv(tecan_sets,here::here("data", "tecan_sets.csv"), row.names = FALSE)
write.csv(lms_sets,here::here("data", "lms_sets.csv"), row.names = FALSE)

#plot linear portions of growth curves
ggplot(tecan_sets, aes(x = time_hr, y = log_A600)) +
  geom_point(aes(shape = strain, color = strain)) +
  geom_smooth(method = lm, se = FALSE, fullrange=TRUE, lwd = 0.5, aes(color = strain))+
  scale_color_viridis(discrete=TRUE)  +
  labs(title = "Growth Curves from 24 hour Tecan Run in LB",
       x = "Time (hr)",
       y = "Log(Absorbance*1000) at 600nm",
       color = "Strain",
       shape = "Strain") +
  theme(text=element_text(family="serif"))

#calculate doubling times using slope of linear models
doubling <- lms_sets %>%
  mutate(#k = slope*2.303, #can shorten into just g
    # g = 0.693/k,
    slope_min = slope/60,
    g = .301/slope_min, #in seconds, slope is 
    doubling_hr = g/60) %>%
  select(strain, slope, doubling_hr, adj_r_sq)
doubling
lms_sets
write.csv(doubling,here::here("data", "doubling.csv"), row.names = FALSE)