# Code to estimate partition coefficients (Kpan)


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
pan_water <- read.csv("Output/Data/Water/PANConcentration.csv")
pcb_solution <- read.csv("Output/Data/Water/SolutionPCBConcentrations.csv")

pan_water <- pan_water %>% select(-X)
pcb_solution <- pcb_solution %>% select(-X)

pan_long <- pan_water %>%
  select(SolutionPCB, sample_id, experiment, time, ends_with("_ng_g_loq_co")) %>%
  pivot_longer(
    cols = ends_with("_ng_g_loq_co"),
    names_to = "PCB",
    values_to = "pan_ng_g"
  ) %>%
  mutate(PCB = gsub("_ng_g_loq_co", "", PCB))

Kpan <- pan_long %>%
  left_join(
    pcb_solution %>% select(SolutionPCB, PCB, mean_ng_L),
    by = c("SolutionPCB", "PCB")
  ) %>%
  mutate(Kpan_L_g = pan_ng_g / mean_ng_L)

# Plots
# PCB4
ggplot(
  Kpan %>% filter(PCB == "PCB4", experiment == "shaking"),
  aes(x = time, y = Kpan_L_g, color = SolutionPCB)) +
  geom_point() +
  facet_wrap(~ SolutionPCB) +
  labs(x = "Time", y = "Kpan (L/g) - PCB4")

# PCB18
ggplot(
  Kpan %>% filter(PCB == "PCB18", experiment == "shaking"),
  aes(x = time, y = Kpan_L_g, color = SolutionPCB)) +
  geom_point() +
  facet_wrap(~ SolutionPCB) +
  labs(x = "Time", y = "Kpan (L/g) - PCB18")

# PCB52
ggplot(
  Kpan %>% filter(PCB == "PCB52", experiment == "shaking"),
  aes(x = time, y = Kpan_L_g, color = SolutionPCB)) +
  geom_point() +
  facet_wrap(~ SolutionPCB) +
  labs(x = "Time", y = "Kpan (L/g) - PCB52")
