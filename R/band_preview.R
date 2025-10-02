#' Preview spectral bands
#'
#' @param .data A numeric vector of wavenumbers OR a data.frame/tibble with spectral data
#' @param spectrum Character. Name of the column containing the spectrum to preview (only for tibble)
#' @param params A tibble created with band_params() containing band parameters
#' @param wn_col Character. Name of the wavenumber column (only for tibble)
#' @param show_sum Logical. Show the sum of all bands? Default is TRUE
#' @param colors Character vector. Colors for bands. NULL uses ggplot2 default colors
#'
#' @return A ggplot object
#' @export
band_preview <- function(.data, spectrum = NULL, params, wn_col = NULL, show_sum = TRUE, colors = NULL) {

  # Determinar se .data é um vetor ou tibble
  is_vector <- is.numeric(.data) && is.null(dim(.data))

  if (is_vector) {
    # É um vetor
    x <- .data
    has_spectrum <- FALSE

  } else {
    # É um tibble - validações padrão tidyspec
    if (missing(.data)) {
      stop("Argument '.data' is missing. Please provide a numeric vector or data.frame/tibble.")
    }
    if (!inherits(.data, c("data.frame", "tbl", "tbl_df"))) {
      stop("Argument '.data' must be a numeric vector or a data.frame/tibble.")
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

    # Verificar se spectrum foi fornecido
    has_spectrum <- !is.null(spectrum)
    if (has_spectrum && !spectrum %in% colnames(.data)) {
      stop(paste0("Column '", spectrum, "' not found in the input data."))
    }

    # Extrair vetor x do tibble
    x <- .data[[wn_col]]

    # Extrair y_observed se fornecido
    if (has_spectrum) {
      y_obs <- .data[[spectrum]]
    }
  }

  # Validar params
  if (missing(params) || !inherits(params, c("data.frame", "tbl", "tbl_df"))) {
    stop("Argument 'params' must be a tibble created with band_params().")
  }

  required_cols <- c("type", "x0", "amp")
  missing_cols <- setdiff(required_cols, names(params))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in params: ", paste(missing_cols, collapse = ", "))
  }

  # Criar band_id se não existir
  if (!"band_id" %in% names(params)) {
    params$band_id <- paste0("band", seq_len(nrow(params)))
  }

  # Validar parâmetros por tipo de banda
  for (i in seq_len(nrow(params))) {
    type <- params$type[i]
    band_id <- params$band_id[i]

    if (type == "gaussian") {
      if (is.na(params$sigma[i]) || is.null(params$sigma[i])) {
        stop("Band '", band_id, "' is gaussian but sigma is missing")
      }
    } else if (type == "lorentzian") {
      if (is.na(params$gamma[i]) || is.null(params$gamma[i])) {
        stop("Band '", band_id, "' is lorentzian but gamma is missing")
      }
    } else if (type == "voigtian") {
      if (is.na(params$sigma[i]) || is.null(params$sigma[i])) {
        stop("Band '", band_id, "' is voigtian but sigma is missing")
      }
      if (is.na(params$gamma[i]) || is.null(params$gamma[i])) {
        stop("Band '", band_id, "' is voigtian but gamma is missing")
      }
    } else if (type == "pseudo_voigtian") {
      if (is.na(params$wid[i]) || is.null(params$wid[i])) {
        stop("Band '", band_id, "' is pseudo_voigtian but wid is missing")
      }
      if (is.na(params$eta[i]) || is.null(params$eta[i])) {
        stop("Band '", band_id, "' is pseudo_voigtian but eta is missing")
      }
    } else {
      stop("Unknown band type: ", type, " for band '", band_id, "'")
    }
  }

  # Calcular as bandas
  bands_data <- tibble(x = x)

  for (i in seq_len(nrow(params))) {
    band_id <- params$band_id[i]
    type <- params$type[i]

    if (type == "gaussian") {
      y <- gaussian(x, x0 = params$x0[i], sigma = params$sigma[i], amp = params$amp[i])
    } else if (type == "lorentzian") {
      y <- lorentzian(x, x0 = params$x0[i], gamma = params$gamma[i], amp = params$amp[i])
    } else if (type == "voigtian") {
      y <- voigtian(x, x0 = params$x0[i], sigma = params$sigma[i],
                    gamma = params$gamma[i], amp = params$amp[i])
    } else if (type == "pseudo_voigtian") {
      y <- pseudo_voigtian(x, x0 = params$x0[i], wid = params$wid[i],
                           eta = params$eta[i], amp = params$amp[i])
    }

    bands_data[[band_id]] <- y
  }

  # Calcular soma se solicitado
  if (show_sum) {
    sum_values <- rowSums(bands_data[, -1, drop = FALSE])
    bands_data$Sum <- sum_values
  }

  # Preparar dados para ggplot (formato long)
  bands_long <- bands_data %>%
    tidyr::pivot_longer(cols = -x, names_to = "band_id", values_to = "intensity")

  # Criar o plot base
  p <- ggplot()

  # Adicionar dados observados como pontos (se fornecidos)
  if (has_spectrum) {
    obs_data <- tibble(x = x, intensity = y_obs)
    p <- p +
      geom_point(data = obs_data, aes(x = x, y = intensity, color = spectrum),
                 size = 1.5, alpha = 0.6)
  }

  # Adicionar bandas como linhas
  p <- p +
    geom_line(data = bands_long, aes(x = x, y = intensity, color = band_id, group = band_id),
              linewidth = 1)

  # Configurar cores
  if (!is.null(colors)) {
    # Adicionar a cor do spectrum observado se necessário
    if (has_spectrum) {
      if (!spectrum %in% names(colors)) {
        colors <- c(setNames("black", spectrum), colors)
      }
    }
    p <- p + scale_color_manual(values = colors)
  }

  # Estilo
  p <- p +
    labs(x = "Wavenumber (cm⁻¹)", y = "Intensity", color = "Legend") +
    theme_minimal() +
    theme(legend.position = "right")

  return(p)
}
