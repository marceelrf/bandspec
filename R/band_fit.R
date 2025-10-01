#' Fit spectral bands using nonlinear regression
#'
#' @param .data A data.frame or tibble with spectral data
#' @param spectrum Character. Name of the column containing the spectrum to fit
#' @param params A tibble created with band_params() containing initial parameters
#' @param wn_col Character. Name of the wavenumber column (uses default if NULL)
#' @param method Character. Fitting method (currently only "nlsLM" supported)
#' @param control List. Control parameters passed to nlsLM
#' @param ... Additional arguments passed to nlsLM
#'
#' @return An object of class 'band_fit' containing fitted parameters and statistics
#' @export
band_fit <- function(.data,
                     spectrum,
                     params,
                     wn_col = NULL,
                     method = "nlsLM",
                     control = list(),
                     ...) {

  # Validações padrão tidyspec
  if (missing(.data)) {
    stop("Argument '.data' is missing. Please provide a data.frame or tibble with spectral data.")
  }
  if (!inherits(.data, c("data.frame", "tbl", "tbl_df"))) {
    stop("Argument '.data' must be a data.frame or tibble.")
  }
  if (is.null(wn_col)) {
    wn_col <- get0(".wn_col_default", envir = tidyspec_env,
                   ifnotfound = NULL)
    if (is.null(wn_col)) {
      stop("wn_col not specified and no default defined with set_spec_wn().")
    } else {
      warn_missing_param_once("wn_col", wn_col)
    }
  }
  if (!wn_col %in% colnames(.data)) {
    stop(paste0("Column '", wn_col, "' not found in the input data."))
  }

  # Validações específicas
  if (missing(spectrum)) {
    stop("Argument 'spectrum' is missing. Please specify which column to fit.")
  }
  if (!spectrum %in% colnames(.data)) {
    stop(paste0("Column '", spectrum, "' not found in the input data."))
  }
  if (missing(params) || !inherits(params, c("data.frame", "tbl", "tbl_df"))) {
    stop("Argument 'params' must be a tibble created with band_params().")
  }
  if (method != "nlsLM") {
    stop("Currently only 'nlsLM' method is supported.")
  }

  # Extrair dados
  x <- .data[[wn_col]]
  y_obs <- .data[[spectrum]]

  # Criar função modelo composta
  model_function <- create_composite_model(params)

  # Criar lista de parâmetros iniciais
  start_params <- create_start_list(params)

  # Criar fórmula
  param_names <- names(start_params)
  formula_str <- paste0("y_obs ~ model_function(x, ",
                        paste(param_names, collapse = ", "), ")")
  formula_obj <- as.formula(formula_str)

  # Fazer o fit
  fit_result <- tryCatch({
    minpack.lm::nlsLM(
      formula = formula_obj,
      start = start_params,
      control = control,
      ...
    )
  }, error = function(e) {
    stop("Fitting failed: ", e$message,
         "\nTry adjusting initial parameters or control settings.")
  })

  # Extrair parâmetros otimizados
  fitted_params <- update_params_with_fit(params, fit_result)

  # Calcular predições
  y_pred <- predict(fit_result)

  # Calcular componentes individuais
  individual_bands <- calculate_individual_bands(x, fitted_params)

  # Calcular resíduos
  residuals <- y_obs - y_pred

  # Calcular métricas
  metrics <- calculate_fit_metrics(y_obs, y_pred, fit_result)

  # Calcular estatísticas das bandas
  band_stats <- calculate_band_statistics(x, fitted_params)

  # Criar objeto band_fit
  result <- structure(
    list(
      data = tibble(
        x = x,
        y_observed = y_obs,
        y_fitted = y_pred,
        residuals = residuals
      ),
      params_initial = params,
      params_fitted = fitted_params,
      individual_bands = individual_bands,
      band_stats = band_stats,
      metrics = metrics,
      fit_object = fit_result,
      spectrum_name = spectrum,
      wn_col = wn_col,
      method = method
    ),
    class = "band_fit"
  )

  return(result)
}

# Função auxiliar: criar modelo composto
create_composite_model <- function(params) {
  function(x, ...) {
    args <- list(...)
    n_bands <- nrow(params)
    result <- numeric(length(x))

    for (i in seq_len(n_bands)) {
      band_id <- params$band_id[i]
      type <- params$type[i]

      # Extrair parâmetros desta banda
      x0 <- args[[paste0(band_id, "_x0")]]
      amp <- args[[paste0(band_id, "_amp")]]

      if (type == "gaussian") {
        sigma <- args[[paste0(band_id, "_sigma")]]
        result <- result + gaussian(x, x0, sigma, amp)

      } else if (type == "lorentzian") {
        gamma <- args[[paste0(band_id, "_gamma")]]
        result <- result + lorentzian(x, x0, gamma, amp)

      } else if (type == "voigtian") {
        sigma <- args[[paste0(band_id, "_sigma")]]
        gamma <- args[[paste0(band_id, "_gamma")]]
        result <- result + voigtian(x, x0, sigma, gamma, amp)

      } else if (type == "pseudo_voigtian") {
        wid <- args[[paste0(band_id, "_wid")]]
        eta <- args[[paste0(band_id, "_eta")]]
        result <- result + pseudo_voigtian(x, amp, x0, wid, eta)
      }
    }

    return(result)
  }
}

# Função auxiliar: criar lista de parâmetros iniciais
create_start_list <- function(params) {
  start_list <- list()

  for (i in seq_len(nrow(params))) {
    band_id <- params$band_id[i]
    type <- params$type[i]

    start_list[[paste0(band_id, "_x0")]] <- params$x0[i]
    start_list[[paste0(band_id, "_amp")]] <- params$amp[i]

    if (type == "gaussian") {
      start_list[[paste0(band_id, "_sigma")]] <- params$sigma[i]

    } else if (type == "lorentzian") {
      start_list[[paste0(band_id, "_gamma")]] <- params$gamma[i]

    } else if (type == "voigtian") {
      start_list[[paste0(band_id, "_sigma")]] <- params$sigma[i]
      start_list[[paste0(band_id, "_gamma")]] <- params$gamma[i]

    } else if (type == "pseudo_voigtian") {
      start_list[[paste0(band_id, "_wid")]] <- params$wid[i]
      start_list[[paste0(band_id, "_eta")]] <- params$eta[i]
    }
  }

  return(start_list)
}

# Função auxiliar: atualizar params com valores fitted
update_params_with_fit <- function(params, fit) {
  fitted_params <- params
  coefs <- coef(fit)

  for (i in seq_len(nrow(params))) {
    band_id <- params$band_id[i]
    type <- params$type[i]

    fitted_params$x0[i] <- coefs[[paste0(band_id, "_x0")]]
    fitted_params$amp[i] <- coefs[[paste0(band_id, "_amp")]]

    if (type == "gaussian") {
      fitted_params$sigma[i] <- coefs[[paste0(band_id, "_sigma")]]

    } else if (type == "lorentzian") {
      fitted_params$gamma[i] <- coefs[[paste0(band_id, "_gamma")]]

    } else if (type == "voigtian") {
      fitted_params$sigma[i] <- coefs[[paste0(band_id, "_sigma")]]
      fitted_params$gamma[i] <- coefs[[paste0(band_id, "_gamma")]]

    } else if (type == "pseudo_voigtian") {
      fitted_params$wid[i] <- coefs[[paste0(band_id, "_wid")]]
      fitted_params$eta[i] <- coefs[[paste0(band_id, "_eta")]]
    }
  }

  return(fitted_params)
}

# Função auxiliar: calcular bandas individuais
calculate_individual_bands <- function(x, params) {
  bands <- tibble(x = x)

  for (i in seq_len(nrow(params))) {
    band_id <- params$band_id[i]
    type <- params$type[i]

    if (type == "gaussian") {
      y <- gaussian(x, params$x0[i], params$sigma[i], params$amp[i])
    } else if (type == "lorentzian") {
      y <- lorentzian(x, params$x0[i], params$gamma[i], params$amp[i])
    } else if (type == "voigtian") {
      y <- voigtian(x, params$x0[i], params$sigma[i], params$gamma[i], params$amp[i])
    } else if (type == "pseudo_voigtian") {
      y <- pseudo_voigtian(x, params$amp[i], params$x0[i], params$wid[i], params$eta[i])
    }

    bands[[band_id]] <- y
  }

  return(bands)
}

# Função auxiliar: calcular métricas do fit
calculate_fit_metrics <- function(y_obs, y_pred, fit) {
  n <- length(y_obs)
  p <- length(coef(fit))

  # Resíduos
  residuals <- y_obs - y_pred

  # SSE, SST, SSR
  sse <- sum(residuals^2)
  sst <- sum((y_obs - mean(y_obs))^2)
  ssr <- sst - sse

  # R²
  r_squared <- 1 - (sse / sst)

  # R² ajustado
  r_squared_adj <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

  # RMSE
  rmse <- sqrt(sse / n)

  # MAE
  mae <- mean(abs(residuals))

  # AIC e BIC
  logLik_val <- logLik(fit)
  aic <- AIC(fit)
  bic <- BIC(fit)

  tibble(
    r_squared = r_squared,
    r_squared_adj = r_squared_adj,
    rmse = rmse,
    mae = mae,
    sse = sse,
    aic = aic,
    bic = bic,
    n_obs = n,
    n_params = p
  )
}

# Função auxiliar: calcular estatísticas das bandas
calculate_band_statistics <- function(x, params) {
  stats_list <- list()

  for (i in seq_len(nrow(params))) {
    band_id <- params$band_id[i]
    type <- params$type[i]

    # Calcular banda individual
    if (type == "gaussian") {
      y <- gaussian(x, params$x0[i], params$sigma[i], params$amp[i])
      fwhm <- 2 * sqrt(2 * log(2)) * params$sigma[i]

    } else if (type == "lorentzian") {
      y <- lorentzian(x, params$x0[i], params$gamma[i], params$amp[i])
      fwhm <- 2 * params$gamma[i]

    } else if (type == "voigtian") {
      y <- voigtian(x, params$x0[i], params$sigma[i], params$gamma[i], params$amp[i])
      # Aproximação de FWHM para Voigt
      fg <- 2 * sqrt(2 * log(2)) * params$sigma[i]
      fl <- 2 * params$gamma[i]
      fwhm <- 0.5346 * fl + sqrt(0.2166 * fl^2 + fg^2)

    } else if (type == "pseudo_voigtian") {
      y <- pseudo_voigtian(x, params$amp[i], params$x0[i], params$wid[i], params$eta[i])
      fwhm <- params$wid[i]
    }

    # Altura do pico
    height <- max(y)

    # Posição do pico
    peak_position <- params$x0[i]

    # Área (integral usando regra do trapézio)
    dx <- diff(x)
    area <- sum((y[-1] + y[-length(y)]) / 2 * dx)

    stats_list[[i]] <- tibble(
      band_id = band_id,
      type = type,
      peak_position = peak_position,
      height = height,
      fwhm = fwhm,
      area = area
    )
  }

  bind_rows(stats_list) %>%
    mutate(
      height_relative = height / sum(height) * 100,
      area_relative = area / sum(area) * 100
    )
}

# Método print
#' @export
print.band_fit <- function(x, ...) {
  cat("Band Fitting Results\n")
  cat("====================\n\n")
  cat("Spectrum:", x$spectrum_name, "\n")
  cat("Method:", x$method, "\n")
  cat("Number of bands:", nrow(x$params_fitted), "\n\n")

  cat("Fit Quality:\n")
  cat("  R² =", sprintf("%.4f", x$metrics$r_squared), "\n")
  cat("  RMSE =", sprintf("%.4f", x$metrics$rmse), "\n\n")

  cat("Use summary() for detailed results\n")
  invisible(x)
}

# Método summary
#' @export
summary.band_fit <- function(object, ...) {
  cat("Band Fitting Summary\n")
  cat("====================\n\n")

  cat("Fitted Parameters:\n")
  print(object$params_fitted, n = Inf)
  cat("\n")

  cat("Fit Metrics:\n")
  print(object$metrics)
  cat("\n")

  cat("Band Statistics:\n")
  print(object$band_stats, n = Inf)

  invisible(object)
}

# Função para extrair parâmetros
#' @export
band_get_params <- function(fit) {
  if (!inherits(fit, "band_fit")) {
    stop("Object must be of class 'band_fit'")
  }
  fit$params_fitted
}

# Função para extrair métricas
#' @export
band_get_metrics <- function(fit) {
  if (!inherits(fit, "band_fit")) {
    stop("Object must be of class 'band_fit'")
  }
  fit$metrics
}

# Função para extrair estatísticas das bandas
#' @export
band_get_stats <- function(fit) {
  if (!inherits(fit, "band_fit")) {
    stop("Object must be of class 'band_fit'")
  }
  fit$band_stats
}

# Função para extrair resíduos
#' @export
band_get_residuals <- function(fit) {
  if (!inherits(fit, "band_fit")) {
    stop("Object must be of class 'band_fit'")
  }
  fit$data %>% select(x, residuals)
}
