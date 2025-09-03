# Download NOAA tidal info for Annapolis
library(tidyverse) 
library(lubridate)

source("nuiscance_flood/getTidalDatums.R")

wl_table <- download6minWlData(station_id = 8575512, startDate = "2000-08-01 00:00", 
                               endDate = "2025-07-31 23:59", datum = "MHHW")
write_csv(wl_table, "nuiscance_flood/Annapolist_2020_2025_mhhw.csv")

wl_table <- read_csv("nuiscance_flood/Annapolist_2020_2025_mhhw.csv")

classify_floods <- wl_table %>% 
  mutate(is_flood = waterLevel > 0.29,
         event_group = consecutive_id(is_flood)) %>% 
  filter(is_flood == T) %>% 
  mutate(flood_id = event_group/2) %>% 
  group_by(event_group) %>% 
  mutate(n_obs = n()) %>% 
  ungroup() %>% 
  filter(n_obs > 1) %>% 
  group_by(flood_id) %>% 
  summarise(dateTimeMin = min(dateTime),
            dateTimeMax = max(dateTime),
            dateTimeMed = median(dateTime),
            maxWaterLevel = max(waterLevel),
            meanWaterLevel = mean(waterLevel)) %>% 
  ungroup() %>% 
  mutate(flood_time_hours =  as.numeric(dateTimeMax-dateTimeMin, units = "hours"))

write_csv(classify_floods, "nuiscance_flood/Floods_Classified.csv")

ggplot(classify_floods, aes(x = dateTimeMed, y = maxWaterLevel)) +
  geom_segment(aes(yend = 0, color = flood_time_hours), lwd = 1.25) +
  scale_colour_continuous(type = "viridis") +
  ylab("Maximum Water Level (m above MHHW)") +
  xlab(NULL)
