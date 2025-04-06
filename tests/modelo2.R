library(minpack.lm)
model_function <- function(x, x0_1, gamma_1, sigma_1, amp_1,
                           x0_2, gamma_2, sigma_2, amp_2,
                           x0_3, gamma_3, sigma_3, amp_3,
                           x0_4, gamma_4, sigma_4, amp_4) {
  voigtian(x, x0_1, sigma_1, gamma_1, amp_1) +
    voigtian(x, x0_2, sigma_2, gamma_2, amp_2) +
    voigtian(x, x0_3, sigma_3, gamma_3, amp_3) +
    voigtian(x, x0_4, sigma_4, gamma_4, amp_4)
}

dados <- tibble(
  x = d1$Wavenumber,
  y = model_function(
    x,
    params$x0[1], params$gamma[1], params$sigma[1], params$amp[1],
    params$x0[2], params$gamma[2], params$sigma[2], params$amp[2],
    params$x0[3], params$gamma[3], params$sigma[3], params$amp[3],
    params$x0[4], params$gamma[4], params$sigma[4], params$amp[4]
  ) + rnorm(length(x), sd = 0.5)  # Adicionando algum ruído
)
start_params <- list(
  x0_1 = 1070, gamma_1 = 5, sigma_1 = 9, amp_1 = 28.5,
  x0_2 = 1130, gamma_2 = 5, sigma_2 = 9, amp_2 = 22.5,
  x0_3 = 1165, gamma_3 = 5, sigma_3 = 9, amp_3 = 3.5,
  x0_4 = 1230, gamma_4 = 9, sigma_4 = 9, amp_4 = 2
)

# Ajuste não linear
fit <- nlsLM(
  CoHA100 ~ model_function(
    Wavenumber,
    x0_1, gamma_1, sigma_1, amp_1,
    x0_2, gamma_2, sigma_2, amp_2,
    x0_3, gamma_3, sigma_3, amp_3,
    x0_4, gamma_4, sigma_4, amp_4
  ),
  data = d1,
  start = start_params,
  control = list(maxiter = 1000)
)

summary(fit)

broom::tidy(fit) %>%
  select(term, estimate) %>%
  mutate(peak = case_when(
    str_detect(term,"_1") ~ "v1",
    str_detect(term,"_2") ~ "v2",
    str_detect(term,"_3") ~ "v3",
    str_detect(term,"_4") ~ "v4"
    )) %>%
  mutate(term = str_split_i(term,"_",1)) %>%
  pivot_wider(names_from = term, values_from = estimate)

# Visualizar o ajuste
d1$predicted <- predict(fit)

ggplot(d1, aes(x = Wavenumber)) +
  geom_point(aes(y = CoHA100), color = "black", alpha = 0.5) +
  geom_line(aes(y = predicted), color = "red", linewidth = 1) +
  labs(title = "Ajuste Voigt com 4 picos",
       x = "Wavenumber", y = "Intensidade")
