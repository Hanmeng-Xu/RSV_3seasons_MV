library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(writexl)
library(reactable)
library(ggplot2)
library(stringr)
library(flextable)
library(zoo)
library(cowplot)
library(gtsummary) 
library(gt)
library(lubridate)
library(openxlsx)
library(rio)
library(patchwork)
library(ggpubr)
library(vroom)

# ------------------------------------------------------------------------------
# RSV activity among 1-4 yr (EPIC)
# ------------------------------------------------------------------------------
# data downloaded from: https://github.com/PopHIVE/Ingest/blob/main/data/epic_resp_infections/standard/weekly.csv.gz
df.epic <- vroom("data/external_rsv_activity/weekly.csv.gz") %>% filter(geography == "09") %>%
  filter(age == "1-4 Years") %>%
  mutate(collection_week = floor_date(time, "week"))

plot(df.epic$epic_pct_rsv, type = "l")



# df.epic <- vroom("data/external_rsv_activity/weekly.csv.gz") %>% filter(geography == "00") %>% # natioanl 
#   filter(age == "1-4 Years" | age == "<1 Years" ) %>%
#   mutate(collection_week = floor_date(time, "week"))
# 
# plt.epic.trend <- df.epic |>
#   filter(time > as.Date("2020-10-01")) |>
#   ggplot(aes(x = time, y = epic_pct_rsv, color = age, group = age)) +
#   geom_line() +
#   scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
#   labs(x = NULL, y = "RSV (% of ED encounters)") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# calculate the log form
df.epic <- df.epic %>% 
  mutate(log_positivity_ext = log(epic_pct_rsv + 0.01)) %>%
  mutate(log_positivity_ext_c = log_positivity_ext - mean(log_positivity_ext, na.rm = TRUE))

plot(df.epic$log_positivity_ext, type = "l")
plot(df.epic$log_positivity_ext_c, type = "l")

# save this to interfile folder 
# saveRDS(df.epic, "data/external_rsv_activity/df.rsv.posrate.ext.pophive.1_4Y.epic.rds")

df.epic %>% 
  filter(time >= as.Date("2023-10-1") & time <= as.Date("2026-3-1")) %>% 
  ggplot() +
  geom_line(aes(x = time, y = epic_pct_rsv))



# ------------------------------------------------------------------------------
# RSV activity in all age grups (EPIC)
# ------------------------------------------------------------------------------
# data downloaded from: https://github.com/PopHIVE/Ingest/blob/main/data/epic_resp_infections/standard/weekly.csv.gz
df.epic <- vroom("data/external_rsv_activity/weekly.csv.gz") %>% filter(geography == "09") %>%
  filter(age == "Total") %>%
  mutate(collection_week = floor_date(time, "week"))

plot(df.epic$epic_pct_rsv, type = "l")



# calculate the log form
df.epic <- df.epic %>% 
  mutate(log_positivity_ext = log(epic_pct_rsv + 0.01)) %>%
  mutate(log_positivity_ext_c = log_positivity_ext - mean(log_positivity_ext, na.rm = TRUE))

plot(df.epic$log_positivity_ext, type = "l")
plot(df.epic$log_positivity_ext_c, type = "l")

# save this to interfile folder 
# saveRDS(df.epic, "data/external_rsv_activity/df.rsv.posrate.ext.pophive.allages.epic.rds")

df.epic %>% 
  filter(time >= as.Date("2023-10-1") & time <= as.Date("2026-3-1")) %>% 
  ggplot() +
  geom_line(aes(x = time, y = epic_pct_rsv))



# ------------------------------------------------------------------------------
# RSV activity in all age grups (NSSP)
# ------------------------------------------------------------------------------
# data downloaded from: https://github.com/PopHIVE/Ingest/blob/main/data/nssp/standard/data.csv.gz
df.nssp <- vroom("../season2_RSV_risk/data/external/data_pophive_nssp.csv.gz") %>% filter(geography == "09") %>%
  mutate(collection_week = floor_date(time, "week"))

plot(df.nssp$percent_visits_rsv, type = "l")

# calculate the log form
df.nssp <- df.nssp %>% 
  mutate(log_positivity_ext = log(percent_visits_rsv + 0.01)) %>%
  mutate(log_positivity_ext_c = log_positivity_ext - mean(log_positivity_ext, na.rm = TRUE))

plot(df.nssp$log_positivity_ext, type = "l")
plot(df.nssp$log_positivity_ext_c, type = "l")


# save this to interfile folder 
# saveRDS(df.nssp, "data/external_rsv_activity/df.rsv.posrate.ext.pophive.allages.nssp.rds")

df.epic %>% 
  filter(time >= as.Date("2023-10-1") & time <= as.Date("2026-3-1")) %>% 
  ggplot() +
  geom_line(aes(x = time, y = epic_pct_rsv))

# ------------------------------------------------------------------------------
# compare the RSV activity from different sources
# ------------------------------------------------------------------------------
ra.nssp <- readRDS("data/external_rsv_activity/df.rsv.posrate.ext.pophive.allages.nssp.rds") %>% 
  dplyr::select(collection_week, log_positivity_ext_c, percent_visits_rsv) %>% 
  dplyr::rename(original_data = percent_visits_rsv)
ra.epic.allages <- readRDS("data/external_rsv_activity/df.rsv.posrate.ext.pophive.allages.epic.rds") %>%
  dplyr::select(collection_week, log_positivity_ext_c, epic_pct_rsv) %>%
  dplyr::rename(original_data = epic_pct_rsv)
ra.epic.1_4Y <- readRDS("data/external_rsv_activity/df.rsv.posrate.ext.pophive.1_4Y.epic.rds") %>%
  dplyr::select(collection_week, log_positivity_ext_c, epic_pct_rsv) %>%
  dplyr::rename(original_data = epic_pct_rsv)

ra <- rbind(ra.nssp %>% mutate(source = "NSSP"),
            ra.epic.allages %>% mutate(source = "EPIC (total)"),
            ra.epic.1_4Y %>% mutate(source = "EPIC (1-4Y)")
            ) %>%
  filter(collection_week >= as.Date("2023-9-1") & collection_week <= as.Date("2026-3-31"))

# compare log rsv activity 
ra %>% 
  ggplot() +
  geom_line(aes(x = collection_week, y = log_positivity_ext_c, color = source, group = source))

# compare original data 
ra %>% 
  ggplot() +
  geom_line(aes(x = collection_week, y = original_data, color = source, group = source))

