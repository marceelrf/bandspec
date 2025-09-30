band_preview <- function(x, params, show_sum = TRUE, colors = NULL) {

  # Validações básicas
  if (!is.data.frame(params)) {
    stop("params must be a data.frame or tibble")
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
    }

    bands_data[[band_id]] <- y
  }

  # Calcular soma se solicitado
  if (show_sum) {
    bands_data$sum <- rowSums(bands_data[, -1, drop = FALSE])
  }

  # Preparar dados para ggplot (formato long)
  bands_long <- bands_data %>%
    tidyr::pivot_longer(cols = -x, names_to = "band_id", values_to = "intensity")

  # Configurar cores
  if (is.null(colors)) {
    # Cor única padrão
    plot_colors <- "steelblue"
  } else if (length(colors) == 1) {
    # Uma cor para todas
    plot_colors <- colors
  } else {
    # Vetor nomeado de cores
    if (is.null(names(colors))) {
      warning("colors vector is not named, will use in order of bands")
    }
    plot_colors <- colors
  }

  # Criar o plot
  p <- ggplot(bands_long, aes(x = x, y = intensity, color = band_id, group = band_id)) +
    geom_line(linewidth = 1) +
    labs(x = "Wavenumber", y = "Intensity", color = "Band") +
    theme_minimal() +
    theme(legend.position = "right")

  # Aplicar cores
  if (length(plot_colors) == 1) {
    p <- p + scale_color_manual(values = rep(plot_colors, nrow(params) + ifelse(show_sum, 1, 0)))
  } else {
    p <- p + scale_color_manual(values = plot_colors)
  }

  # Destacar soma se presente
  if (show_sum) {
    p <- p +
      geom_line(data = filter(bands_long, band_id == "sum"),
                aes(x = x, y = intensity),
                color = "black", linewidth = 1.5, linetype = "solid")
  }

  return(p)
}
