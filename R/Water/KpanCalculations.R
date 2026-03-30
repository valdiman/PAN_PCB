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

# < 30 days
ggplot(
  Kpan %>% filter(PCB == "PCB4", experiment == "shaking"),
  aes(x = time, y = Kpan_L_g, color = SolutionPCB)) +
  geom_point() +
  xlim(c(0, 30)) +
  facet_wrap(~ SolutionPCB) +
  labs(x = "Time", y = "Kpan (L/g) - PCB4")

# PCB18
ggplot(
  Kpan %>% filter(PCB == "PCB18", experiment == "shaking"),
  aes(x = time, y = Kpan_L_g, color = SolutionPCB)) +
  geom_point() +
  facet_wrap(~ SolutionPCB) +
  labs(x = "Time", y = "Kpan (L/g) - PCB18")

# < 30 days
ggplot(
  Kpan %>% filter(PCB == "PCB18", experiment == "shaking"),
  aes(x = time, y = Kpan_L_g, color = SolutionPCB)) +
  geom_point() +
  xlim(c(0, 30)) +
  facet_wrap(~ SolutionPCB) +
  labs(x = "Time", y = "Kpan (L/g) - PCB18")

# PCB52
ggplot(
  Kpan %>% filter(PCB == "PCB52", experiment == "shaking"),
  aes(x = time, y = Kpan_L_g, color = SolutionPCB)) +
  geom_point() +
  facet_wrap(~ SolutionPCB) +
  labs(x = "Time", y = "Kpan (L/g) - PCB52")

# < 30 days
ggplot(
  Kpan %>% filter(PCB == "PCB52", experiment == "shaking"),
  aes(x = time, y = Kpan_L_g, color = SolutionPCB)) +
  geom_point() +
  xlim(c(0, 30)) +
  facet_wrap(~ SolutionPCB) +
  labs(x = "Time", y = "Kpan (L/g) - PCB52")

# All together
# PCB4
ggplot(
  Kpan %>% filter(PCB == "PCB4", experiment == "shaking"),
  aes(x = time, y = log10(Kpan_L_g), color = SolutionPCB)) +
  geom_point() +
  labs(x = "Time", y = "log10 Kpan (L/g) - PCB4")

# PCB18
ggplot(
  Kpan %>% filter(PCB == "PCB18", experiment == "shaking"),
  aes(x = time, y = log10(Kpan_L_g), color = SolutionPCB)) +
  geom_point() +
  labs(x = "Time", y = "log10 Kpan (L/g) - PCB18")

# PCB52
ggplot(
  Kpan %>% filter(PCB == "PCB52", experiment == "shaking"),
  aes(x = time, y = log10(Kpan_L_g), color = SolutionPCB)) +
  geom_point() +
  labs(x = "Time", y = "log10 Kpan (L/g) - PCB52")

# Uptake model (Kpan) -----------------------------------------------------
# Remove s2, S4
Kpan2 <- Kpan %>%
  filter(!SolutionPCB %in% c("S2", "S4"))

# PCB4
m_PCB4 <- nls(
  Kpan_L_g ~ Kpan * (1 - exp(-time * ku)),
  data = subset(Kpan2, PCB == "PCB4"),
  start = list(Kpan = 10, ku = 0.01)
)

summary(m_PCB4)
d4 <- subset(Kpan2, PCB == "PCB4")

cols <- as.numeric(factor(d4$SolutionPCB))

plot(d4$time, d4$Kpan_L_g,
     col = cols,
     pch = 16,
     xlab = "Time",
     ylab = "Kpan_L_g (PCB 4)")

curve(
  coef(m_PCB4)["Kpan"] * (1 - exp(-x * coef(m_PCB4)["ku"])),
  add = TRUE,
  col = "red",
  lwd = 2
)

legend("topright",
       legend = levels(factor(d4$SolutionPCB)),
       col = 1:length(unique(d4$SolutionPCB)),
       pch = 16, cex = 0.8)

# PCB18
m_PCB18 <- nls(
  Kpan_L_g ~ Kpan * (1 - exp(-time * ku)),
  data = subset(Kpan2, PCB == "PCB18"),
  start = list(Kpan = 10, ku = 0.01)
)

summary(m_PCB18)
d18 <- subset(Kpan2, PCB == "PCB18")

cols <- as.numeric(factor(d18$SolutionPCB))

plot(d18$time, d18$Kpan_L_g,
     col = cols,
     pch = 16,
     xlab = "Time",
     ylab = "Kpan_L_g (PCB 18)")

curve(
  coef(m_PCB18)["Kpan"] * (1 - exp(-x * coef(m_PCB18)["ku"])),
  add = TRUE,
  col = "red",
  lwd = 2
)

legend("topright",
       legend = levels(factor(d18$SolutionPCB)),
       col = 1:length(unique(d18$SolutionPCB)),
       pch = 16, cex = 0.8)

# PCB52
m_PCB52 <- nls(
  Kpan_L_g ~ Kpan * (1 - exp(-time * ku)),
  data = subset(Kpan2, PCB == "PCB52"),
  start = list(Kpan = 10, ku = 0.01)
)

summary(m_PCB52)
d52 <- subset(Kpan2, PCB == "PCB52")

cols <- as.numeric(factor(d52$SolutionPCB))

plot(d52$time, d52$Kpan_L_g,
     col = cols,
     pch = 16,
     xlab = "Time",
     ylab = "Kpan_L_g (PCB 52)")

curve(
  coef(m_PCB52)["Kpan"] * (1 - exp(-x * coef(m_PCB52)["ku"])),
  add = TRUE,
  col = "red",
  lwd = 2
)

legend("topright",
       legend = levels(factor(d52$SolutionPCB)),
       col = 1:length(unique(d52$SolutionPCB)),
       pch = 16, cex = 0.8)

# Shaking
Kpan_shaking <- Kpan %>%
  filter(experiment == "shaking")
Kpan_shaking <- Kpan_shaking %>%
  filter(!SolutionPCB %in% c("S4"))

# PCB4
m_PCB4 <- nls(
  Kpan_L_g ~ Kpan * (1 - exp(-time * ku)),
  data = subset(Kpan_shaking, PCB == "PCB4"),
  start = list(Kpan = 10, ku = 0.01)
)

summary(m_PCB4)
d4 <- subset(Kpan_shaking, PCB == "PCB4")

cols <- as.numeric(factor(d4$SolutionPCB))

plot(d4$time, d4$Kpan_L_g,
     col = cols,
     pch = 16,
     xlab = "Time",
     ylab = "Kpan_L_g (PCB 4)")

curve(
  coef(m_PCB4)["Kpan"] * (1 - exp(-x * coef(m_PCB4)["ku"])),
  add = TRUE,
  col = "red",
  lwd = 2
)

legend("topright",
       legend = levels(factor(d4$SolutionPCB)),
       col = 1:length(unique(d4$SolutionPCB)),
       pch = 16, cex = 0.8)

# PCB18
m_PCB18 <- nls(
  Kpan_L_g ~ Kpan * (1 - exp(-time * ku)),
  data = subset(Kpan_shaking, PCB == "PCB18"),
  start = list(Kpan = 10, ku = 0.01)
)

summary(m_PCB18)
d18 <- subset(Kpan_shaking, PCB == "PCB18")

cols <- as.numeric(factor(d18$SolutionPCB))

plot(d18$time, d18$Kpan_L_g,
     col = cols,
     pch = 16,
     xlab = "Time",
     ylab = "Kpan_L_g (PCB 18)")

curve(
  coef(m_PCB18)["Kpan"] * (1 - exp(-x * coef(m_PCB18)["ku"])),
  add = TRUE,
  col = "red",
  lwd = 2
)

legend("topright",
       legend = levels(factor(d18$SolutionPCB)),
       col = 1:length(unique(d18$SolutionPCB)),
       pch = 16, cex = 0.8)

# PCB52
m_PCB52 <- nls(
  Kpan_L_g ~ Kpan * (1 - exp(-time * ku)),
  data = subset(Kpan_shaking, PCB == "PCB52"),
  start = list(Kpan = 10, ku = 0.01)
)

summary(m_PCB52)
d52 <- subset(Kpan_shaking, PCB == "PCB52")

cols <- as.numeric(factor(d52$SolutionPCB))

plot(d52$time, d52$Kpan_L_g,
     col = cols,
     pch = 16,
     xlab = "Time",
     ylab = "Kpan_L_g (PCB 52)")

curve(
  coef(m_PCB52)["Kpan"] * (1 - exp(-x * coef(m_PCB52)["ku"])),
  add = TRUE,
  col = "red",
  lwd = 2
)

legend("topright",
       legend = levels(factor(d52$SolutionPCB)),
       col = 1:length(unique(d52$SolutionPCB)),
       pch = 16, cex = 0.8)





