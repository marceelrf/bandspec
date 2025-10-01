# band_preview <- function(x, params, show_sum = TRUE, colors = NULL) {
#
#   # Validações básicas
#   if (!is.data.frame(params)) {
#     stop("params must be a data.frame or tibble")
#   }
#
#   required_cols <- c("type", "x0", "amp")
#   missing_cols <- setdiff(required_cols, names(params))
#   if (length(missing_cols) > 0) {
#     stop("Missing required columns in params: ", paste(missing_cols, collapse = ", "))
#   }
#
#   # Criar band_id se não existir
#   if (!"band_id" %in% names(params)) {
#     params$band_id <- paste0("band", seq_len(nrow(params)))
#   }
#
#   # Validar parâmetros por tipo de banda
#   for (i in seq_len(nrow(params))) {
#     type <- params$type[i]
#     band_id <- params$band_id[i]
#
#     if (type == "gaussian") {
#       if (is.na(params$sigma[i]) || is.null(params$sigma[i])) {
#         stop("Band '", band_id, "' is gaussian but sigma is missing")
#       }
#     } else if (type == "lorentzian") {
#       if (is.na(params$gamma[i]) || is.null(params$gamma[i])) {
#         stop("Band '", band_id, "' is lorentzian but gamma is missing")
#       }
#     } else if (type == "voigtian") {
#       if (is.na(params$sigma[i]) || is.null(params$sigma[i])) {
#         stop("Band '", band_id, "' is voigtian but sigma is missing")
#       }
#       if (is.na(params$gamma[i]) || is.null(params$gamma[i])) {
#         stop("Band '", band_id, "' is voigtian but gamma is missing")
#       }
#     } else if (type == "pseudo_voigtian") {
#       if (is.na(params$wid[i]) || is.null(params$wid[i])) {
#         stop("Band '", band_id, "' is pseudo_voigtian but wid is missing")
#       }
#       if (is.na(params$eta[i]) || is.null(params$eta[i])) {
#         stop("Band '", band_id, "' is pseudo_voigtian but eta is missing")
#       }
#     } else {
#       stop("Unknown band type: ", type, " for band '", band_id, "'")
#     }
#   }
#
#   # Calcular as bandas
#   bands_data <- tibble(x = x)
#
#   for (i in seq_len(nrow(params))) {
#     band_id <- params$band_id[i]
#     type <- params$type[i]
#
#     if (type == "gaussian") {
#       y <- gaussian(x, x0 = params$x0[i], sigma = params$sigma[i], amp = params$amp[i])
#     } else if (type == "lorentzian") {
#       y <- lorentzian(x, x0 = params$x0[i], gamma = params$gamma[i], amp = params$amp[i])
#     } else if (type == "voigtian") {
#       y <- voigtian(x, x0 = params$x0[i], sigma = params$sigma[i],
#                     gamma = params$gamma[i], amp = params$amp[i])
#     } else if (type == "pseudo_voigtian") {
#       y <- pseudo_voigtian(x, x0 = params$x0[i], wid = params$wid[i],
#                            eta = params$eta[i], amp = params$amp[i])
#     }
#
#     bands_data[[band_id]] <- y
#   }
#
#   # Calcular soma se solicitado
#   if (show_sum) {
#     bands_data$sum <- rowSums(bands_data[, -1, drop = FALSE])
#   }
#
#   # Preparar dados para ggplot (formato long)
#   bands_long <- bands_data %>%
#     tidyr::pivot_longer(cols = -x, names_to = "band_id", values_to = "intensity")
#
#   # Configurar cores
#   if (is.null(colors)) {
#     # Usar paleta de cores padrão do ggplot2
#     n_bands <- nrow(params)
#     if (show_sum) n_bands <- n_bands + 1
#     plot_colors <- NULL  # deixar ggplot escolher
#   } else if (length(colors) == 1) {
#     # Uma cor para todas
#     plot_colors <- rep(colors, nrow(params) + ifelse(show_sum, 1, 0))
#     names(plot_colors) <- c(params$band_id, if(show_sum) "sum")
#   } else {
#     # Vetor de cores (nomeado ou não)
#     plot_colors <- colors
#   }
#
#   # Criar o plot
#   p <- ggplot(bands_long, aes(x = x, y = intensity, color = band_id, group = band_id)) +
#     geom_line(linewidth = 1) +
#     labs(x = "Wavenumber (cm⁻¹)", y = "Intensity", color = "Band") +
#     theme_minimal() +
#     theme(legend.position = "right")
#
#   # Aplicar cores se fornecidas
#   if (!is.null(plot_colors)) {
#     p <- p + scale_color_manual(values = plot_colors)
#   }
#
#   # Destacar soma se presente
#   if (show_sum) {
#     sum_data <- filter(bands_long, band_id == "sum")
#     p <- p +
#       geom_line(data = sum_data,
#                 aes(x = x, y = intensity),
#                 color = "black", linewidth = 1.2, linetype = "solid",
#                 inherit.aes = FALSE)
#   }
#
#   return(p)
# }

band_preview <- function(x, params, wn_col = "x",
                         show_sum = TRUE, colors = NULL,
                         exp_color = "grey40", exp_alpha = 0.6, exp_size = 1.5) {

  # valida params
  if (!is.data.frame(params)) stop("params must be a data.frame or tibble")
  required_cols <- c("type", "x0", "amp")
  missing_cols <- setdiff(required_cols, names(params))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in params: ", paste(missing_cols, collapse = ", "))
  }
  if (!"band_id" %in% names(params)) {
    params$band_id <- paste0("band", seq_len(nrow(params)))
  }

  # checar se x é tibble do tidyspec (com colunas wn_col e y)
  experimental_data <- NULL
  if (is.data.frame(x)) {
    if (!(wn_col %in% names(x))) stop("Coluna ", wn_col, " não encontrada no tibble")
    if (!("y" %in% names(x))) stop("O tibble deve ter uma coluna 'y' com os valores experimentais")

    x_vals <- x[[wn_col]]
    experimental_data <- tibble::tibble(x = x_vals, y = x$y)

  } else if (is.numeric(x)) {
    x_vals <- x
  } else {
    stop("x deve ser vetor numérico ou tibble do tidyspec")
  }

  # calcular bandas
  bands_data <- tibble::tibble(x = x_vals)
  for (i in seq_len(nrow(params))) {
    band_id <- params$band_id[i]
    type <- params$type[i]
    if (type == "gaussian") {
      y <- gaussian(x_vals, x0 = params$x0[i], sigma = params$sigma[i], amp = params$amp[i])
    } else if (type == "lorentzian") {
      y <- lorentzian(x_vals, x0 = params$x0[i], gamma = params$gamma[i], amp = params$amp[i])
    } else if (type == "voigtian") {
      y <- voigtian(x_vals, x0 = params$x0[i], sigma = params$sigma[i],
                    gamma = params$gamma[i], amp = params$amp[i])
    } else if (type == "pseudo_voigtian") {
      y <- pseudo_voigtian(x_vals, x0 = params$x0[i], wid = params$wid[i],
                           eta = params$eta[i], amp = params$amp[i])
    }
    bands_data[[band_id]] <- y
  }

  if (show_sum) {
    bands_data$sum <- rowSums(bands_data[, -1, drop = FALSE])
  }

  bands_long <- tidyr::pivot_longer(bands_data, cols = -x,
                                    names_to = "band_id", values_to = "intensity")

  # base plot com bandas
  p <- ggplot(bands_long, aes(x = x, y = intensity, color = band_id)) +
    geom_line(linewidth = 1) +
    labs(x = wn_col, y = "Intensity", color = "Band") +
    theme_minimal()

  # adicionar pontos experimentais
  if (!is.null(experimental_data)) {
    p <- p +
      geom_point(data = experimental_data,
                 aes(x = x, y = y),
                 color = exp_color,
                 alpha = exp_alpha,
                 size = exp_size,
                 inherit.aes = FALSE)
  }

  # destaque soma
  if (show_sum) {
    sum_data <- dplyr::filter(bands_long, band_id == "sum")
    p <- p +
      geom_line(data = sum_data, aes(x = x, y = intensity),
                color = "black", linewidth = 1.2)
  }

  return(p)
}

