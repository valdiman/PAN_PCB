# Initial code to check PAN PCB data in water



# load library
library(dplyr)

# Read data ---------------------------------------------------------------
pan_water <- read.csv("Data/Water/data.csv")

# Blanks
pan_water_blank <- pan_water %>%
  filter(experiment == "blank")

pan_water_blank$PCB4_ng_g <- pan_water_blank$PCB4 / pan_water_blank$mass
pan_water_blank$PCB18_ng_g <- pan_water_blank$PCB18 / pan_water_blank$mass
pan_water_blank$PCB52_ng_g <- pan_water_blank$PCB52 / pan_water_blank$mass

# Replace Inf for NA for next steps
pan_water_blank$PCB4_ng_g[is.infinite(pan_water_blank$PCB4_ng_g)] <- NA
pan_water_blank$PCB18_ng_g[is.infinite(pan_water_blank$PCB18_ng_g)] <- NA
pan_water_blank$PCB52_ng_g[is.infinite(pan_water_blank$PCB52_ng_g)] <- NA

# Histograms
# PCB4
hist(pan_water_blank$PCB4_ng_g)
hist(log10(pan_water_blank$PCB4_ng_g))
# PCB18
hist(pan_water_blank$PCB18_ng_g)
hist(log10(pan_water_blank$PCB18_ng_g))
# PCB52
hist(pan_water_blank$PCB52_ng_g)
hist(log10(pan_water_blank$PCB52_ng_g))





