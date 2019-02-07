#' @title Lamb wave dispersion function for symmetric modes.
#' @description A utility function that calculates the Lamb wave dispersion function for symmetric modes
#'              according Wilcox's PhD thesis eq. 2.15, p. 28.
#' @param x The wavenumber
#' @param par The parameters thickness, longitudinal velocity, transvers velocity and angular frequency
#'            as a \code{list} object.
#' @return The value of the Lamb wave dispersion function for symmetric modes according
#'         Wilcox's PhD thesis eq. 2.15, p. 28.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname delta_S
#' @export
delta_S <- function(x, par) {
  # extract variables
  k <- x
  # extract parameters
  d <- par[[1]]
  v_l <- par[[2]]
  v_t <- par[[3]]
  omega <- par[[4]]
  # compute auxillary parameters
  k_l <- as.complex(omega/v_l) # longitudinal wave number (Wilcox's PhD thesis eq. 2.10, p. 27)
  k_t <- as.complex(omega/v_t) # transverse wave number (Wilcox's PhD thesis eq. 2.10, p. 27)
  q_l <- sqrt(k^2 - k_l^2) # Wilcox's PhD thesis eq. 2.12, p. 27
  q_t <- sqrt(k^2 - k_t^2) # Wilcox's PhD thesis eq. 2.12, p. 27
  # Wilcox's PhD thesis eq. 2.15, p. 28
  delta_S <- (k^2 + q_t^2)^2*cosh(q_l*d/2)*sinh(q_t*d/2) -
    4*k^2*q_t*q_l*cosh(q_t*d/2)*sinh(q_l*d/2)
}
