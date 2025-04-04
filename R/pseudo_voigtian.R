#' Pseudo-Voigt Peak Profile Function
#'
#' Computes a Pseudo-Voigt profile, which is an approximation of the Voigt profile
#' (convolution of Gaussian and Lorentzian profiles) using a linear combination
#' of the two distributions.
#'
#' @param x Numeric vector, the x-values (e.g., wavelength, energy, or frequency values)
#' @param amp Numeric, the amplitude (peak height) at the center position
#' @param x0 Numeric, the center position of the peak
#' @param wid Numeric, the full width at half maximum (FWHM) of the peak. Must be positive.
#' @param eta Numeric between 0 and 1, the mixing parameter:
#'   - 0 = pure Gaussian profile
#'   - 1 = pure Lorentzian profile
#'   - Values between represent mixtures
#'
#' @return Numeric vector of the same length as x containing the computed profile values
#'
#' @details
#' The Pseudo-Voigt profile is given by:
#' \deqn{PseudoVoigt(x) = A[\eta L(x) + (1-\eta)G(x)]}
#' where:
#' \itemize{
#'   \item \eqn{A} is the amplitude (amp)
#'   \item \eqn{L(x)} is the Lorentzian component
#'   \item \eqn{G(x)} is the Gaussian component
#'   \item \eqn{\eta} (eta) is the mixing parameter
#' }
#'
#' The function uses `RcppFaddeeva` for efficient computation of the base profiles,
#' with automatic conversion between FWHM and the native width parameters of each distribution.
#'
#' @examples
#' x <- seq(-10, 10, length.out = 500)
#' # Pure Gaussian profile
#' y_g <- pseudo_voigtian(x, amp = 1, x0 = 0, wid = 2, eta = 0)
#' # Pure Lorentzian profile
#' y_l <- pseudo_voigtian(x, amp = 1, x0 = 0, wid = 2, eta = 1)
#' # Mixed profile (eta = 0.5)
#' y_pv <- pseudo_voigtian(x, amp = 1, x0 = 0, wid = 2, eta = 0.5)
#'
#' plot(x, y_g, type = "l", col = "blue", ylim = c(0, 1))
#' lines(x, y_l, col = "red")
#' lines(x, y_pv, col = "purple", lwd = 2)
#' legend("topright", legend = c("Gaussian (η=0)", "Lorentzian (η=1)", "Pseudo-Voigt (η=0.5)"),
#'        col = c("blue", "red", "purple"), lty = 1)
#'
#' @seealso
#' \code{\link[RcppFaddeeva]{Gauss}} and \code{\link[RcppFaddeeva]{Lorentz}} for the base functions,
#' \code{\link[stats]{nls}} for curve fitting applications
#'
#' @export
pseudo_voigtian <- function(x, amp, x0, wid, eta) {
  # Verificação dos parâmetros
  if(eta < 0 | eta > 1) stop("eta must be between 0 and 1")
  if(any(wid <= 0)) stop("wid must be positive")

  # Convertendo wid para os parâmetros específicos de cada função
  sigma_g <- wid / (2*sqrt(2*log(2)))  # Converte FWHM para sigma
  gamma_l <- wid / 2                    # Converte FWHM para gamma

  # Obtendo os perfis (verifique a documentação do RcppFaddeeva para confirmar)
  g <- RcppFaddeeva::Gauss(x = x, x0 = x0, sigma = sigma_g)
  l <- RcppFaddeeva::Lorentz(x = x, x0 = x0, gamma = gamma_l)

  # Combinação linear
  amp * (eta * l + (1 - eta) * g)
}
