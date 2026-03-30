# Initial code to check SPME PCB data in water


# Install Packages
{
  install.packages("dplyr")
  install.packages("ggplot2")
}

# Load Libraries
{
  library(dplyr)
  library(ggplot2)
  library(tidyr)
}

# Read data ---------------------------------------------------------------
spme_water <- read.csv("Data/Water/dataSPME.csv")

# Blank Analysis ----------------------------------------------------------
spme_water_blank <- spme_water %>%
  filter(experiment == "blank")

spme_water_blank$PCB4_ng_cm <- spme_water_blank$PCB4 / spme_water_blank$length
spme_water_blank$PCB18_ng_cm <- spme_water_blank$PCB18 / spme_water_blank$length
spme_water_blank$PCB52_ng_cm <- spme_water_blank$PCB52 / spme_water_blank$length

# Histograms
# PCB4
hist(spme_water_blank$PCB4_ng_cm)
hist(log10(spme_water_blank$PCB4_ng_cm))
# PCB18
hist(spme_water_blank$PCB18_ng_cm)
hist(log10(spme_water_blank$PCB18_ng_cm))
# PCB52
hist(spme_water_blank$PCB52_ng_cm)
hist(log10(spme_water_blank$PCB52_ng_cm))

# Q-Q plots
# PCB4
qqnorm(spme_water_blank$PCB4_ng_cm, main = "Concentration (ng/cm)")
qqline(spme_water_blank$PCB4_ng_cm)
qqnorm(log10(spme_water_blank$PCB4_ng_cm), main = "Concentration (ng/cm)/log10")
qqline(log10(spme_water_blank$PCB4_ng_cm))
# PCB18
qqnorm(spme_water_blank$PCB18_ng_cm, main = "Concentration (ng/cm)")
qqline(spme_water_blank$PCB18_ng_cm)
qqnorm(log10(spme_water_blank$PCB18_ng_cm), main = "Concentration (ng/cm)/log10")
qqline(log10(spme_water_blank$PCB18_ng_cm))
# PCB52
qqnorm(spme_water_blank$PCB52_ng_cm, main = "Concentration (ng/cm)")
qqline(spme_water_blank$PCB52_ng_cm)
qqnorm(log10(spme_water_blank$PCB52_ng_cm), main = "Concentration (ng/cm)/log10")
qqline(log10(spme_water_blank$PCB52_ng_cm))

# Shapiro test
# PCB4
shapiro.test(spme_water_blank[, 10])$p.value
shapiro.test(log10(spme_water_blank[, 10]))$p.value
# PCB18
shapiro.test(spme_water_blank[, 11])$p.value
shapiro.test(log10(spme_water_blank[, 11]))$p.value
# PCB52
shapiro.test(spme_water_blank[, 12])$p.value
shapiro.test(log10(spme_water_blank[, 12]))$p.value

# Calculate LOQ of the log10
# Upper 95 CI% (=mean + 1.96*sd/(n)^0.5)
log_PCB4_bl <- log10(spme_water_blank$PCB4_ng_cm)
n <- sum(!is.na(log_PCB4_bl))
loq_log4_bl <- mean(log_PCB4_bl, na.rm = TRUE) +
  1.96 * sd(log_PCB4_bl, na.rm = TRUE) / sqrt(n)
loq_4 <- 10^loq_log4_bl

log_PCB18_bl <- log10(spme_water_blank$PCB18_ng_cm)
n <- sum(!is.na(log_PCB18_bl))
loq_log18_bl <- mean(log_PCB18_bl, na.rm = TRUE) +
  1.96 * sd(log_PCB18_bl, na.rm = TRUE) / sqrt(n)
loq_18 <- 10^loq_log18_bl

log_PCB52_bl <- log10(spme_water_blank$PCB52_ng_cm)
n <- sum(!is.na(log_PCB52_bl))
loq_log52_bl <- mean(log_PCB52_bl, na.rm = TRUE) +
  1.96 * sd(log_PCB52_bl, na.rm = TRUE) / sqrt(n)
loq_52 <- 10^loq_log52_bl

# Sample Analysis ---------------------------------------------------------
# Retrieve samples
spme_water_samples <- spme_water %>%
  filter(experiment != "blank")

spme_water_samples$PCB4_ng_cm <- spme_water_samples$PCB4 / spme_water_samples$length
spme_water_samples$PCB18_ng_cm <- spme_water_samples$PCB18 / spme_water_samples$length
spme_water_samples$PCB52_ng_cm <- spme_water_samples$PCB52 / spme_water_samples$length

# Comparison between samples and loq @ log10 scale
# PCB4
log_PCB4_sa <- log10(spme_water_samples$PCB4_ng_cm)
spme_water_samples$PCB4_ng_cm_loq_co <- ifelse(
  is.na(log_PCB4_sa),
  NA,
  ifelse(log_PCB4_sa > loq_log4_bl, 10^(log_PCB4_sa), NA)
)

# PCB18
log_PCB18_sa <- log10(spme_water_samples$PCB18_ng_cm)
spme_water_samples$PCB18_ng_cm_loq_co <- ifelse(
  is.na(log_PCB18_sa),
  NA,
  ifelse(log_PCB18_sa > loq_log18_bl, 10^(log_PCB18_sa), NA)
)

# PCB52
log_PCB52_sa <- log10(spme_water_samples$PCB52_ng_cm)
spme_water_samples$PCB52_ng_cm_loq_co <- ifelse(
  is.na(log_PCB52_sa),
  NA,
  ifelse(log_PCB52_sa > loq_log52_bl, 10^(log_PCB52_sa), NA)
)

# Plots
# Shaking
ggplot(spme_water_samples,
  aes(x = time, y = PCB4_ng_cm_loq_co, color = SolutionPCB)) +
  geom_point() +
  facet_wrap(~ SolutionPCB)

ggplot(spme_water_samples,
       aes(x = time, y = PCB18_ng_cm_loq_co, color = SolutionPCB)) +
  geom_point() +
  facet_wrap(~ SolutionPCB)

ggplot(spme_water_samples,
  aes(x = time, y = PCB52_ng_cm_loq_co, color = SolutionPCB)) +
  geom_point() +
  facet_wrap(~ SolutionPCB)

# Estimate concentrations -------------------------------------------------
# Cw = _ng_cm_loq_co / (Kspme * Vspme) [ng/L]

Vspme <- 0.000000069 # Lspme/cmspme

Kspme.4 <- 18845.3 # Lw/Lspme
Kspme.18 <- 97795.57 # Lw/Lspme
Kspme.52 <- 200314.7 # Lw/Lspme

spme_water_samples$PCB4_ng_L_loq_co <- spme_water_samples$PCB4_ng_cm_loq_co / (Kspme.4 * Vspme)
spme_water_samples$PCB18_ng_L_loq_co <- spme_water_samples$PCB18_ng_cm_loq_co / (Kspme.18 * Vspme)
spme_water_samples$PCB52_ng_L_loq_co <- spme_water_samples$PCB52_ng_cm_loq_co / (Kspme.52 * Vspme)

# Need to remove sample: SPME_S05_WP_S. Extraction issue
spme_water_samples <- spme_water_samples[-2,]

# Summary concentrations
summary_spme <- spme_water_samples %>%
  select(SolutionPCB, ends_with("_ng_L_loq_co")) %>%
  pivot_longer(
    cols = ends_with("_ng_L_loq_co"),
    names_to = "PCB",
    values_to = "value"
  ) %>%
  mutate(PCB = gsub("_ng_L_loq_co", "", PCB)) %>%
  group_by(SolutionPCB, PCB) %>%
  summarise(
    mean_ng_L = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    rsd_percent = 100 * sd / mean_ng_L,
    n = sum(!is.na(value)),
    .groups = "drop"
  )

summary_spme$PCB <- factor(
  summary_spme$PCB,
  levels = c("PCB4", "PCB18", "PCB52")
)

summary_spme <- summary_spme %>%
  arrange(SolutionPCB, PCB)

# Export data
write.csv(summary_spme, file  = "Output/Data/Water/SolutionPCBConcentrations.csv")

