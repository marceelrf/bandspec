#' Voigt Profile Function
#'
#' Computes the Voigt profile, which is the convolution of a Gaussian profile and
#' a Lorentzian profile. This is commonly used to model spectral lines in spectroscopy,
#' X-ray diffraction patterns, and other physical phenomena where both Gaussian
#' (instrumental) and Lorentzian (natural) broadening effects are present.
#'
#' @param x Numeric vector of x-values where the function will be evaluated
#'          (e.g., wavelength, energy, or frequency values)
#' @param x0 Numeric, the center position of the peak (location parameter)
#' @param sigma Numeric, the standard deviation of the Gaussian component
#'              (related to Gaussian broadening). Must be positive.
#' @param gamma Numeric, the half-width at half-maximum (HWHM) of the Lorentzian
#'              component (related to natural broadening). Must be positive.
#' @param amp Numeric, the amplitude (peak height) at the center position
#'
#' @return Numeric vector of the same length as x containing the computed Voigt
#'         profile values
#'
#' @details
#' The Voigt profile is given by:
#' \deqn{V(x; x_0, \sigma, \gamma) = A \cdot \frac{1}{\sigma\sqrt{2\pi}} \int_{-\infty}^{\infty} \frac{e^{-(x'-x_0)^2/(2\sigma^2)}}{\gamma^2 + (x - x')^2} dx'}
#' where:
#' \itemize{
#'   \item \eqn{A} is the amplitude (amp)
#'   \item \eqn{x_0} is the center position (x0)
#'   \item \eqn{\sigma} controls the Gaussian broadening
#'   \item \eqn{\gamma} controls the Lorentzian broadening
#' }
#'
#' This implementation uses the efficient computation from the \code{RcppFaddeeva} package.
#' The full width at half maximum (FWHM) of the Voigt profile can be approximated by:
#' \deqn{FWHM_V \approx 0.5346 FWHM_L + \sqrt{0.2166 FWHM_L^2 + FWHM_G^2}}
#' where \eqn{FWHM_L = 2\gamma} and \eqn{FWHM_G = 2\sigma\sqrt{2\ln 2}}.
#'
#' @examples
#' x <- seq(-10, 10, length.out = 500)
#' # Voigt profile with dominant Gaussian broadening
#' y_v1 <- voigtian(x, x0 = 0, sigma = 1, gamma = 0.2, amp = 1)
#' # Voigt profile with dominant Lorentzian broadening
#' y_v2 <- voigtian(x, x0 = 0, sigma = 0.2, gamma = 1, amp = 1)
#'
#' plot(x, y_v1, type = "l", col = "blue",
#'      main = "Voigt Profile Examples", ylab = "Intensity")
#' lines(x, y_v2, col = "red")
#' legend("topright",
#'        legend = c(expression(paste(sigma, "=1.0, ", gamma, "=0.2")),
#'                    expression(paste(sigma, "=0.2, ", gamma, "=1.0"))),
#'        col = c("blue", "red"), lty = 1)
#'
#' @seealso
#' \code{\link[RcppFaddeeva]{Voigt}} for the underlying implementation,
#' \code{\link{pseudo_voigtian}} for a faster approximation,
#' \code{\link[stats]{nls}} for curve fitting applications
#'
#' @references
#' The algorithm is based on:
#' Zaghloul, M. R. (2007). "On the calculation of the Voigt line profile".
#' Monthly Notices of the Royal Astronomical Society.
#'
#' @export

voigtian <- function(x, x0, sigma, gamma, amp){

  amp * RcppFaddeeva::Voigt(x = x, x0 = x0, sigma = sigma, gamma = gamma)
}
