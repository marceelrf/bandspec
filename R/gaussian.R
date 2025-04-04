#' Gaussian Profile Function
#'
#' Computes a Gaussian (normal) distribution profile, commonly used to model
#' instrumental broadening effects in spectroscopy and other measurement systems.
#'
#' @param x Numeric vector of x-values where the function will be evaluated
#' @param x0 Numeric, the center position of the peak (mean)
#' @param sigma Numeric, the standard deviation of the distribution. Must be positive.
#' @param amp Numeric, the amplitude (peak height) at the center position
#'
#' @return Numeric vector of the same length as x containing the computed Gaussian
#'         profile values
#'
#' @details
#' The Gaussian profile is given by:
#' \deqn{G(x; x_0, \sigma) = A \cdot e^{-\frac{(x - x_0)^2}{2\sigma^2}}}
#' where:
#' \itemize{
#'   \item \eqn{A} is the amplitude (amp)
#'   \item \eqn{x_0} is the mean/center position
#'   \item \eqn{\sigma} is the standard deviation
#' }
#' The full width at half maximum (FWHM) is \eqn{2\sigma\sqrt{2\ln 2} ≈ 2.355\sigma.
#'
#' @examples
#' x <- seq(-5, 5, length.out = 200)
#' y <- gaussian(x, x0 = 0, sigma = 1, amp = 1)
#' plot(x, y, type = "l", main = "Gaussian Profile", ylab = "Intensity")
#'
#' @seealso
#' \code{\link[RcppFaddeeva]{Gauss}} for the underlying implementation,
#' \code{\link{lorentzian}} for Lorentzian profile,
#' \code{\link{voigtian}} for Voigt profile
#'
#' @export
gaussian <- function(x, x0, sigma, amp){

  amp * RcppFaddeeva::Gauss(x = x, x0 = x0, sigma = sigma)
}
