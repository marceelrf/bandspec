#' Lorentzian Profile Function
#'
#' Computes a Lorentzian (Cauchy) distribution profile, commonly used to model
#' natural line shapes in spectroscopy and other physical phenomena with
#' homogeneous broadening.
#'
#' @param x Numeric vector of x-values where the function will be evaluated
#' @param x0 Numeric, the center position of the peak (location parameter)
#' @param gamma Numeric, the half-width at half-maximum (HWHM) of the distribution.
#'              Must be positive.
#' @param amp Numeric, the amplitude (peak height) at the center position
#'
#' @return Numeric vector of the same length as x containing the computed Lorentzian
#'         profile values
#'
#' @details
#' The Lorentzian profile is given by:
#' \deqn{L(x; x_0, \gamma) = A \cdot \frac{\gamma^2}{(x - x_0)^2 + \gamma^2}}
#' where:
#' \itemize{
#'   \item \eqn{A} is the amplitude (amp)
#'   \item \eqn{x_0} is the center position
#'   \item \eqn{\gamma} is the HWHM (half-width at half-maximum)
#' }
#' The full width at half maximum (FWHM) is \eqn{2\gamma}.
#'
#' @examples
#' x <- seq(-5, 5, length.out = 200)
#' y <- lorentzian(x, x0 = 0, gamma = 1, amp = 1)
#' plot(x, y, type = "l", main = "Lorentzian Profile", ylab = "Intensity")
#'
#' @seealso
#' \code{\link[RcppFaddeeva]{Lorentz}} for the underlying implementation,
#' \code{\link{gaussian}} for Gaussian profile,
#' \code{\link{voigtian}} for Voigt profile
#'
#' @export
lorentzian <- function(x, x0, gamma, amp){

  amp * RcppFaddeeva::Lorentz(x = x, x0 = x0, gamma = gamma)
}
