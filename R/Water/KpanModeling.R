# Initial code to model uptake PCB in PAN
# Shaking experiments

# Install Packages
{
  install.packages("dplyr")
  install.packages("ggplot2")
  install.packages("stringr")
  install.packages("tidyr")
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

# Remove exta rows
dat <- dat %>%
  filter(sample_id != "sample_id")

# Long format
dat_long <- dat %>%
  pivot_longer(
    cols = c(`PCB4`, `PCB18.30`, `PCB52`),
    names_to = "compound",
    values_to = "concentration"
  )

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
      unit == "ng/100 mL" ~ "ng/L"
    )
  )

# PCB4
pcb52 <- dat_long %>%
  filter(
    compound == "PCB52",
    !is.na(time_hr)
    )

# PAN
pan <- dat_long %>%
  filter(
    compound == "PCB52",
    unit_new == "ng/kg",
    !is.na(time_hr)
  ) %>%
  mutate(replicate = str_extract(sample_id, "R\\d")) %>%
  select(time_hr, replicate, Cpan = concentration_new)

water <- dat_long %>%
  filter(
    compound == "PCB52",
    unit_new == "ng/L",
    !is.na(time_hr)
  ) %>%
  mutate(replicate = str_extract(sample_id, "R\\d")) %>%
  select(time_hr, replicate, Cw = concentration_new)

model_dat <- left_join(
  pan,
  water,
  by = c("time_hr", "replicate")
)

model_dat <- model_dat %>%
  mutate(
    Kobs = Cpan / Cw
  )

summary(model_dat$Kobs)

ggplot(model_dat, aes(time_hr, log10(Kobs))) +
  geom_point(size = 3) +
  theme_bw()

fit <- nlsLM(
  Cpan ~ Kpan * Cw * (1 - exp(-ke * time_hr)),
  data = model_dat,
  start = list(
    Kpan = max(model_dat$Kobs),
    ke = 0.5
  )
)

summary(fit)

model_dat <- model_dat %>%
  mutate(
    Kobs = Cpan / Cw
  )

coef_fit <- coef(fit)

pred_dat <- data.frame(
  time_hr = seq(
    min(model_dat$time_hr),
    max(model_dat$time_hr),
    length.out = 1000
  )
)

pred_dat$Kpred <- coef_fit["Kpan"] *
  (1 - exp(-coef_fit["ke"] * pred_dat$time_hr))

ggplot(model_dat, aes(time_hr, Kobs)) +
  geom_point(size = 3) +
  geom_line(
    data = pred_dat,
    aes(time_hr, Kpred),
    linewidth = 1,
    color = "blue"
  ) +
  theme_bw() +
  scale_y_log10() +
  labs(
    x = "Time (h)",
    y = expression(K[PAN]~"(L kg"^{-1}*")")
  )


