# Helpers individuais
band_gaussian <- function(band_id = NULL, x0, sigma, amp) {
  tibble(
    band_id = band_id %||% NA_character_,
    type = "gaussian",
    x0 = x0,
    sigma = sigma,
    gamma = NA_real_,
    wid = NA_real_,
    eta = NA_real_,
    amp = amp
  )
}

band_lorentzian <- function(band_id = NULL, x0, gamma, amp) {
  tibble(
    band_id = band_id %||% NA_character_,
    type = "lorentzian",
    x0 = x0,
    sigma = NA_real_,
    gamma = gamma,
    wid = NA_real_,
    eta = NA_real_,
    amp = amp
  )
}

band_voigtian <- function(band_id = NULL, x0, sigma, gamma, amp) {
  tibble(
    band_id = band_id %||% NA_character_,
    type = "voigtian",
    x0 = x0,
    sigma = sigma,
    gamma = gamma,
    wid = NA_real_,
    eta = NA_real_,
    amp = amp
  )
}

band_pseudo_voigtian <- function(band_id = NULL, x0, wid, eta, amp) {
  tibble(
    band_id = band_id %||% NA_character_,
    type = "pseudo_voigtian",
    x0 = x0,
    sigma = NA_real_,
    gamma = NA_real_,
    wid = wid,
    eta = eta,
    amp = amp
  )
}

# Função principal para combinar
band_params <- function(...) {
  bands <- list(...)

  # Combinar todos os tibbles
  result <- bind_rows(bands)

  # Gerar band_ids automáticos se necessário
  if (any(is.na(result$band_id))) {
    missing_ids <- is.na(result$band_id)
    result$band_id[missing_ids] <- paste0("band", seq_len(sum(missing_ids)))
  }

  # Reordenar colunas para melhor visualização
  result <- result %>%
    select(band_id, type, x0, amp, everything())

  return(result)
}
