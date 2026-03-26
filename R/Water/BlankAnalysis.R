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

# Q-Q plots
# PCB4
qqnorm(pan_water_blank$PCB4_ng_g, main = "Concentration (ng/g)")
qqline(pan_water_blank$PCB4_ng_g)
qqnorm(log10(pan_water_blank$PCB4_ng_g), main = "Concentration (ng/g)/log10")
qqline(log10(pan_water_blank$PCB4_ng_g))
# PCB18
qqnorm(pan_water_blank$PCB18_ng_g, main = "Concentration (ng/g)")
qqline(pan_water_blank$PCB18_ng_g)
qqnorm(log10(pan_water_blank$PCB18_ng_g), main = "Concentration (ng/g)/log10")
qqline(log10(pan_water_blank$PCB18_ng_g))
# PCB52
qqnorm(pan_water_blank$PCB52_ng_g, main = "Concentration (ng/g)")
qqline(pan_water_blank$PCB52_ng_g)
qqnorm(log10(pan_water_blank$PCB52_ng_g), main = "Concentration (ng/g)/log10")
qqline(log10(pan_water_blank$PCB52_ng_g))

# Shapiro test
# PCB4
shapiro.test(pan_water_blank[, 10])$p.value
shapiro.test(log10(pan_water_blank[, 10]))$p.value
# PCB18
shapiro.test(pan_water_blank[, 11])$p.value
shapiro.test(log10(pan_water_blank[, 11]))$p.value
# PCB52
shapiro.test(pan_water_blank[, 12])$p.value
shapiro.test(log10(pan_water_blank[, 12]))$p.value

# Calculate LOQ of the log10
# Upper 95 CI% (=mean + 1.96*sd/(n)^0.5)
log_PCB4 <- log10(pan_water_blank$PCB4_ng_g)
n <- sum(!is.na(log_PCB4))
loq_log4 <- mean(log_PCB4, na.rm = TRUE) +
  1.96 * sd(log_PCB4, na.rm = TRUE) / sqrt(n)
loq_4 <- 10^loq_log4

log_PCB18 <- log10(pan_water_blank$PCB18_ng_g)
n <- sum(!is.na(log_PCB18))
loq_log18 <- mean(log_PCB18, na.rm = TRUE) +
  1.96 * sd(log_PCB18, na.rm = TRUE) / sqrt(n)
loq_18 <- 10^loq_log18

log_PCB52 <- log10(pan_water_blank$PCB52_ng_g)
n <- sum(!is.na(log_PCB52))
loq_log52 <- mean(log_PCB52, na.rm = TRUE) +
  1.96 * sd(log_PCB52, na.rm = TRUE) / sqrt(n)
loq_52 <- 10^loq_log52

# Retrieve samples
pan_water_samples <- pan_water %>%
  filter(experiment != "blank")

pan_water_samples$PCB4_ng_g <- pan_water_samples$PCB4 / pan_water_samples$mass
pan_water_samples$PCB18_ng_g <- pan_water_samples$PCB18 / pan_water_samples$mass
pan_water_samples$PCB52_ng_g <- pan_water_samples$PCB52 / pan_water_samples$mass

