#' @title Generating excitability function data.
#' @description A function that generates excitability function data.
#' @param dispersionData The dispersion data dataset, Default: 'titaniumDD'
#' @param plot If \code{TRUE}, plots the excitability function data, Default: FALSE
#' @return A dataset with excitability function data embedded into the \code{pongaR} library.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @seealso
#'  \code{\link[numDeriv]{grad}}
#' @rdname makeEF
#' @export
#' @importFrom numDeriv grad
makeEF <- function(dispersionData = "titaniumDD",
                   plot = FALSE) {

  longVelocity <- tranVelocity <- freq <- thickness <- phasVelocity <- NULL

  ## sanity check
  if (is.character(dispersionData)) {
    data(list = dispersionData, envir = environment())
  } else {
    stop("Argument 'dispersionData' must be a character string.")
  }
  if (!is.logical(plot)) {
    stop("Argument 'plot' must be an object of type logical")
  }

  thickness <- eval(parse(text = paste0(dispersionData, "$", "thickness")))
  density <- eval(parse(text = paste0(dispersionData, "$", "density")))
  v_l <- eval(parse(text = paste0(dispersionData, "$", "longVelocity")))
  v_t <- eval(parse(text = paste0(dispersionData, "$", "tranVelocity")))
  freq <- eval(parse(text = paste0(dispersionData, "$", "freq")))
  phasVelocity <- eval(parse(text = paste0(dispersionData, "$", "phasVelocity")))

  omega <- 2*pi*freq*v_t/thickness # rad/s
  v_ph <- phasVelocity*v_t         # m/s
  wavenumber <- omega/v_ph         # rad/m

  k_l <- as.complex(omega/v_l) # longitudinal wave number (Wilcox's PhD thesis eq. 2.10, p. 27)
  k_t <- as.complex(omega/v_t) # transverse wave number (Wilcox's PhD thesis eq. 2.10, p. 27)

  q_l <- sqrt((wavenumber^2 - k_l^2)) # Wilcox's PhD thesis eq. 2.12, p. 27
  q_t <- sqrt((wavenumber^2 - k_t^2)) # Wilcox's PhD thesis eq. 2.12, p. 27

  nModes <- ncol(wavenumber)

  # determinant derivative data
  ddelta <- matrix(data = NA, nrow = nrow(freq), ncol = nModes)

  # symmetric determinant
  for (mode in seq(1, nModes, 2)) {
    ddelta[!is.na(wavenumber[, mode]), mode] <- numDeriv::grad(func = delta_S,
                                                               x = wavenumber[!is.na(wavenumber[, mode]), mode],
                                                               method = "Richardson",
                                                               par = list(thickness,
                                                                          v_l,
                                                                          v_t,
                                                                          omega[!is.na(wavenumber[, mode]), mode]))
  }

  # anti-symmetric determinant
  for (mode in seq(2, nModes, 2)) {
    ddelta[!is.na(wavenumber[, mode]), mode] <- numDeriv::grad(func = delta_A,
                                                               x = wavenumber[!is.na(wavenumber[, mode]), mode],
                                                               method = "Richardson",
                                                               par = list(thickness,
                                                                          v_l,
                                                                          v_t,
                                                                          omega[!is.na(wavenumber[, mode]), mode]))
  }

  # make dimensionless
  ddelta <- ddelta*thickness^3
  wavenumber <- wavenumber*thickness
  k <- wavenumber
  q_l <- q_l*thickness
  q_t <- q_t*thickness

  # Wilcox's PhD thesis eq. 5.18, p. 135
  EF_c <- matrix(data = NA, nrow = nrow(freq), ncol = nModes)
  EF_c[, seq(1, nModes, 2)] <- k[, seq(1, nModes, 2)]*q_l[, seq(1, nModes, 2)]*((q_t[, seq(1, nModes, 2)]^2 - k[, seq(1, nModes, 2)]^2)*sinh(q_t[, seq(1, nModes, 2)]/2)*sinh(q_l[, seq(1, nModes, 2)]/2))/4/ddelta[, seq(1, nModes, 2)]
  EF_c[, seq(2, nModes, 2)] <- k[, seq(2, nModes, 2)]*q_l[, seq(2, nModes, 2)]*((q_t[, seq(2, nModes, 2)]^2 - k[, seq(2, nModes, 2)]^2)*cosh(q_t[, seq(2, nModes, 2)]/2)*cosh(q_l[, seq(2, nModes, 2)]/2))/4/ddelta[, seq(2, nModes, 2)]

  # patches for Wilcox
  EF_c[Re(EF_c) < 0] <- NA # no negative values
  # EF_c[Re(EF_c) > 1] <- NA # cut-off frequency first 2 modes
  # EF_c[freq[, 3] < 0.9, 3] <- NA # cut-off frequency third mode

  # Wilcox's PhD thesis eq. 5.20, p. 136
  EF_s <- 2*EF_c/(1+1i)/k

  # Wilcox's PhD thesis fig. 5.1a, p. 146
  if (plot) {
    plot(freq[, 1],
         Re(EF_s[, 1]),
         type = "l",
         xlab = "Frequency [-]",
         ylab = "Excitability [-]",
         ylim = c(0, 0.5),
         col = "green")
    lines(freq[, 2:dim(freq)[2]],
          Re(EF_s[, 2:dim(EF_s)[2]]),
          col = "green")
    # grid()
  }

  # Wilcox's PhD thesis fig. 5.2, p. 147
  if (plot) {
    plot(freq[, 1],
         Re(EF_c[, 1]),
         type = "l",
         xlab = "Frequency [-]",
         ylab = "Excitability [-]",
         ylim = c(0, 1),
         col = "green")
    lines(freq[, 2:dim(freq)[2]],
          Re(EF_c[, 2:dim(EF_c)[2]]),
          col = "green")
    # grid()
  }

  # make excitability function dataset
  dataName <- paste0(strsplit(dispersionData, split = "DD", fixed = T)[[1]][1],
                     "EF")
  filePath <- paste0("data/", dataName, ".rda")
  data <- list(density = density,
               longVelocity = v_l,
               tranVelocity = v_t,
               freq = freq,
               wavenumber = wavenumber,
               EF_s = EF_s,
               EF_c = EF_c)
  assign(dataName, data)
  save(list = dataName, file = filePath)

  # filePath <- paste("data/", substr(dispersionData, start = 1, stop = nchar(dispersionData) - 2), "EF.rda", sep = "")
  # save(density, longVelocity, tranVelocity, freq, wavenumber, EF_s, EF_c, file = filePath)
}
