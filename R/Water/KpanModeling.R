# Initial code to model uptake PCB in PAN
# Shaking experiments

# Install Packages
{
  install.packages("dplyr")
  install.packages("ggplot2")
  install.packages("stringr")
  install.packages("tidyr")
  install.packages("minpack.lm")
}

# Load Libraries
{
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(minpack.lm)
}

# Read data ---------------------------------------------------------------
dat <- read.csv("Data/Water/data_Uptake_Shaking.csv")

# Remove extra rows
dat <- dat %>%
  filter(sample_id != "sample_id")

# Long format
dat_long <- dat %>%
  pivot_longer(
    cols = c(`PCB4`, `PCB18.30`, `PCB52`),
    names_to = "compound",
    values_to = "concentration")

# Change format and units
dat_long <- dat_long %>%
  mutate(
    concentration = as.numeric(concentration),
    mass_mg = as.numeric(mass_mg),
    time_hr = as.numeric(time_hr)
  ) %>%
  mutate(
    concentration_new = case_when(
      unit == "ng" ~ concentration / mass_mg * 1e6,  # ng/kg
      unit == "ng/100 mL" ~ concentration * 10       # ng/L
    ),
    unit_new = case_when(
      unit == "ng" ~ "ng/kg",
      unit == "ng/100 mL" ~ "ng/L"))

# Select PCB, 4, 18.30 or 52
pcb52 <- dat_long %>% # select PCB
  filter(
    compound == "PCB52", # select PCB
    !is.na(time_hr))

# PAN
pan <- dat_long %>%
  filter(
    compound == "PCB52", # select PCB
    unit_new == "ng/kg",
    !is.na(time_hr)
  ) %>%
  mutate(replicate = str_extract(sample_id, "R\\d")) %>%
  select(time_hr, replicate, Cpan = concentration_new)

water <- dat_long %>%
  filter(
    compound == "PCB52", # select PCB
    unit_new == "ng/L",
    !is.na(time_hr)
  ) %>%
  mutate(replicate = str_extract(sample_id, "R\\d")) %>%
  select(time_hr, replicate, Cw = concentration_new)

# Prepare modeling dataset ------------------------------------------------
model_dat <- left_join(
  pan,
  water,
  by = c("time_hr", "replicate")
) %>%
  mutate(
    Kpan_obs = Cpan / Cw)

# Quick inspection
ggplot(model_dat, aes(time_hr, log10(Kpan_obs))) +
  geom_point(size = 3) +
  theme_bw()

# Fit first-order uptake model --------------------------------------------
fitK <- nlsLM(
  Kpan_obs ~ Kpan * (1 - exp(-ke * time_hr)),
  data = model_dat,
  start = list(
    Kpan = max(model_dat$Kpan_obs),
    ke = 0.5))

summary(fitK)

# Extract parameter estimates ---------------------------------------------
coef_table <- summary(fitK)$coefficients

Kpan <- coef_table["Kpan", "Estimate"]
Kpan_SE <- coef_table["Kpan", "Std. Error"]

ke <- coef_table["ke", "Estimate"]
ke_SE <- coef_table["ke", "Std. Error"]

t90 <- log(10) / ke
t90_SE <- log(10) * ke_SE / ke^2

# Predict fitted curve ----------------------------------------------------
pred_dat <- data.frame(
  time_hr = seq(
    min(model_dat$time_hr),
    max(model_dat$time_hr),
    length.out = 1000))

pred_dat <- pred_dat %>%
  mutate(
    Kpan_pred = Kpan * (1 - exp(-ke * time_hr)))

# Plot fitted model -------------------------------------------------------
plot.k <- ggplot(model_dat, aes(time_hr, Kpan_obs)) +
  geom_point(size = 3, shape = 21) +
  geom_line(data = pred_dat, aes(time_hr, Kpan_pred), linewidth = 0.5,
            color = "black") +
  annotate("text", x = 1, y = 20000, label = "PCB 52", hjust = 0,
           size = 5) +
  scale_y_log10() +
  theme_classic() +
  labs(x = "t (h)", y = expression(K[PAN]~"(L kg"^{-1}*")")) +
  theme(
    axis.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right")

# See plot
plot.k

# Export plot
ggsave("Output/Plots/Water/PCB52_shaking.png",
       plot = plot.k, width = 6, height = 5, dpi = 500)

# Model diagnostics -------------------------------------------------------
model_dat <- model_dat %>%
  mutate(
    fitted = predict(fitK),
    residuals = residuals(fitK))

# Residuals vs time
ggplot(model_dat, aes(time_hr, residuals)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()

# Residuals vs fitted
ggplot(model_dat, aes(fitted, residuals)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()

# Goodness of fit ---------------------------------------------------------
RMSE <- sqrt(mean(model_dat$residuals^2))

RSS <- sum(model_dat$residuals^2)
TSS <- sum((model_dat$Kpan_obs - mean(model_dat$Kpan_obs))^2)

R2 <- 1 - RSS / TSS

fit_stats <- data.frame(
  compound = "PCB52",
  RMSE = RMSE,
  R2 = R2)

fit_stats

# Export parameter estimates ----------------------------------------------
results <- data.frame(
  compound = "PCB52",
  Kpan = Kpan,
  Kpan_SE = Kpan_SE,
  ke = ke,
  ke_SE = ke_SE,
  t90_hr = t90,
  t90_SE = t90_SE,
  RMSE = RMSE,
  R2 = R2)

write.csv(results, "Output/Data/Water/PCB52_kinetic_parameters_shaking.csv",
          row.names = FALSE)









