library(tidyverse)
x <- seq(400, 4000, by = 2)

params <- tibble(
  band_id = c("OH", "CH", "CO"),
  type = c("gaussian", "lorentzian", "voigtian"),
  x0 = c(3400, 2900, 1700),
  amp = c(1.0, 0.8, 1.2),
  sigma = c(100, NA, 50),
  gamma = c(NA, 80, 30)
)

# Uso básico
band_preview(x, params)

# Sem soma
band_preview(x, params, show_sum = FALSE)

# Com cores customizadas
cores <- c(OH = "blue", CH = "red", CO = "green", sum = "black")
band_preview(x, params, colors = cores)
