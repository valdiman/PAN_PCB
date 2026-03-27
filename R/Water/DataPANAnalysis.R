# Initial code to check PAN PCB data in water


# Install Packages
{
  install.packages("dplyr")
  install.packages("ggplot2")
}

# Load Libraries
{
  library(dplyr)
  library(ggplot2)
}

# Read data ---------------------------------------------------------------
pan_water <- read.csv("Data/Water/dataPAN.csv")

# Blank Analysis ----------------------------------------------------------
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
log_PCB4_bl <- log10(pan_water_blank$PCB4_ng_g)
n <- sum(!is.na(log_PCB4_bl))
loq_log4_bl <- mean(log_PCB4_bl, na.rm = TRUE) +
  1.96 * sd(log_PCB4_bl, na.rm = TRUE) / sqrt(n)
loq_4 <- 10^loq_log4_bl

log_PCB18_bl <- log10(pan_water_blank$PCB18_ng_g)
n <- sum(!is.na(log_PCB18_bl))
loq_log18_bl <- mean(log_PCB18_bl, na.rm = TRUE) +
  1.96 * sd(log_PCB18_bl, na.rm = TRUE) / sqrt(n)
loq_18 <- 10^loq_log18_bl

log_PCB52_bl <- log10(pan_water_blank$PCB52_ng_g)
n <- sum(!is.na(log_PCB52_bl))
loq_log52_bl <- mean(log_PCB52_bl, na.rm = TRUE) +
  1.96 * sd(log_PCB52_bl, na.rm = TRUE) / sqrt(n)
loq_52 <- 10^loq_log52_bl

# Sample Analysis ---------------------------------------------------------
# Retrieve samples
pan_water_samples <- pan_water %>%
  filter(experiment != "blank")

pan_water_samples$PCB4_ng_g <- pan_water_samples$PCB4 / pan_water_samples$mass
pan_water_samples$PCB18_ng_g <- pan_water_samples$PCB18 / pan_water_samples$mass
pan_water_samples$PCB52_ng_g <- pan_water_samples$PCB52 / pan_water_samples$mass

# Comparison between samples and loq @ log10 scale
# PCB4
log_PCB4_sa <- log10(pan_water_samples$PCB4_ng_g)
pan_water_samples$PCB4_ng_g_loq_co <- ifelse(
  is.na(log_PCB4_sa),
  NA,
  ifelse(log_PCB4_sa > loq_log4_bl, 10^(log_PCB4_sa), NA)
)

# PCB18
log_PCB18_sa <- log10(pan_water_samples$PCB18_ng_g)
pan_water_samples$PCB18_ng_g_loq_co <- ifelse(
  is.na(log_PCB18_sa),
  NA,
  ifelse(log_PCB18_sa > loq_log18_bl, 10^(log_PCB18_sa), NA)
)

# PCB52
log_PCB52_sa <- log10(pan_water_samples$PCB52_ng_g)
pan_water_samples$PCB52_ng_g_loq_co <- ifelse(
  is.na(log_PCB52_sa),
  NA,
  ifelse(log_PCB52_sa > loq_log52_bl, 10^(log_PCB52_sa), NA)
)

# Plots
# Shaking
ggplot(
  pan_water_samples %>% filter(experiment == "shaking"),
  aes(x = time, y = PCB4_ng_g_loq_co, color = SolutionPCB)
) +
  geom_point() +
  facet_wrap(~ SolutionPCB)

ggplot(
  pan_water_samples %>% filter(experiment == "shaking"),
  aes(x = time, y = PCB18_ng_g_loq_co, color = SolutionPCB)
) +
  geom_point() +
  facet_wrap(~ SolutionPCB)

ggplot(
  pan_water_samples %>% filter(experiment == "shaking"),
  aes(x = time, y = PCB52_ng_g_loq_co, color = SolutionPCB)
) +
  geom_point() +
  facet_wrap(~ SolutionPCB)

# No-shaking
ggplot(
  pan_water_samples %>% filter(experiment == "no_shaking"),
  aes(x = time, y = PCB4_ng_g_loq_co, color = SolutionPCB)
) +
  geom_point() +
  facet_wrap(~ SolutionPCB)

ggplot(
  pan_water_samples %>% filter(experiment == "no_shaking"),
  aes(x = time, y = PCB18_ng_g_loq_co, color = SolutionPCB)
) +
  geom_point() +
  facet_wrap(~ SolutionPCB)

ggplot(
  pan_water_samples %>% filter(experiment == "no_shaking"),
  aes(x = time, y = PCB52_ng_g_loq_co, color = SolutionPCB)
) +
  geom_point() +
  facet_wrap(~ SolutionPCB)

# Export data
write.csv(pan_water_samples, file = "Output/Data/Water/PANConcentration.csv")



