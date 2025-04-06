library(tidyverse)
library(tidyspec)

set_spec_wn("Wavenumber")

Region <- CoHAspec %>%
  tidyspec::spec_smooth_sga() %>%
  spec_blc_rollingBall(wn_col = "Wavenumber",
                       Wn_min = 1030,
                       Wn_max = 1285,
                       ws = 10, wm = 50) %>%
  spec_norm_01()

Region %>%
  tidyspec::spec_smartplot(geom = "line")

#Use spec_smartplotly to find the peaks centers

Region %>%
  spec_select(CoHA100) %>%
  spec_smartplotly()

# Region %>%
#   write_csv2(file = "region_coha100.csv")

v1 <- \(x) voigtian(x, x0 = 1071, sigma = 9, gamma = 5,amp = 28.5)
v2 <- \(x) voigtian(x, x0 = 1130, sigma = 9, gamma = 5,amp = 22.5)
v3 <- \(x) voigtian(x, x0 = 1165, sigma = 9, gamma = 5,amp = 3.5)
v4 <- \(x) voigtian(x, x0 = 1234, sigma = 9, gamma = 9,amp = 2)

Region %>%
  spec_select(CoHA100) %>%
  spec_smartplot() +
  geom_function(fun = v1) +
  geom_function(fun = v2) +
  geom_function(fun = v3) +
  geom_function(fun = v4)

#Initial parameters
params <- tibble(x0 = c(1070, 1130, 1172, 1230),
                 gamma = c(5,5,5,9),
                 sigma = c(9,9,9,9),
                 amp = c(28.5, 22.5, 3.5, 2))


# Função com 4 picos do tipo Voigt
four_peaks <- function(x, x0_1, x0_2, x0_3, x0_4,
                       amp_1, amp_2, amp_3, amp_4,
                       sigma_1, sigma_2, sigma_3, sigma_4,
                       gamma_1, gamma_2, gamma_3, gamma_4) {
  voigtian(x, x0 = x0_1, amp = amp_1, sigma = sigma_1, gamma = gamma_1) +
    voigtian(x, x0 = x0_2, amp = amp_2, sigma = sigma_2, gamma = gamma_2) +
    voigtian(x, x0 = x0_3, amp = amp_3, sigma = sigma_3, gamma = gamma_3) +
    voigtian(x, x0 = x0_4, amp = amp_4, sigma = sigma_4, gamma = gamma_4)
}

# Seleção de espectro da região
d1 <- Region %>%
  spec_select(CoHA100)

# Ajuste do modelo não linear (NLS)

Wavenumber <- d1$Wavenumber
y <- d1$CoHA100
model <- nls(
  y ~ four_peaks(x = Wavenumber,
                 x0_1 = x0_1,x0_2 = x0_2,x0_3 = x0_3,x0_4 = x0_4,
                       amp_1, amp_2, amp_3, amp_4,
                       sigma_1, sigma_2, sigma_3, sigma_4,
                       gamma_1, gamma_2, gamma_3, gamma_4),
  data = d1,
  start = list(
    x0_1 = 1070, x0_2 = 1130, x0_3 = 1172, x0_4 = 1230,
    amp_1 = 28.5, amp_2 = 22.5, amp_3 = 3.5, amp_4 = 2,
    sigma_1 = 9, sigma_2 = 9, sigma_3 = 9, sigma_4 = 9,
    gamma_1 = 5, gamma_2 = 5, gamma_3 = 5, gamma_4 = 9
  )
)
