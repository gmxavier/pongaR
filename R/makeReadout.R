#' @title Generating readout impulse response.
#' @description A function that generates readout impulse response.
#' @param sensorName The sensor dataset name, Default: 'r15iuc'
#' @param filterRange The band-pass filter range, in Hz, Default: c(NULL, NULL)
#' @param amplifierGain The amplfier gain in dB, Default: 0
#' @param plot If \code{TRUE}, plots the readout impulse response, Default: FALSE
#' @return The readout impulse response waveform as a \code{readout} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname makeReadout
#' @export
makeReadout <- function(sensorName = "r15iuc",
                        filterRange = c(NULL, NULL),
                        amplifierGain = 0,
                        plot = FALSE) {
  FR <- NULL

  ## sanity check
  if (is.character(sensorName)) {
    data(list = sensorName, envir = environment())
    FR <- eval(parse(text = sensorName))
  } else {
    stop("Argument 'sensorName' must be a character string.")
  }
  if (!is.null(filterRange)) {
    if (!is.vector(x = filterRange, mode = "numeric") | any(filterRange < 0)) {
      stop("Argument 'filterRange' must be a non negative numeric vector.")
    }
    if (diff(filterRange) < 0) {
      stop("High frequency must be greater than low frequency.")
    }
  }
  if (!is.numeric(amplifierGain) || (amplifierGain < 0)) {
    stop("Argument 'amplifierGain' must be a non negative number.")
  }
  if (!is.logical(plot)) {
    stop("Argument 'plot' must be an object of type logical")
  }

  ## define variables
  # samp.rate <- 2*FR$freq[length(FR$freq)]
  samp.rate <- 2*FR[nrow(FR), 1]

  ## calculate filter frequency response
  if (!is.null(filterRange)) {
    n <- 100 # filter order
    W <- filterRange/(samp.rate/2) # filter critical frequencies
    f <- c(0, W[1], W[1], W[2], W[2], 1)
    m <- c(0, 0, 1, 1, 0, 0)
    filterFR <- freqz(fir2(100, f, m))
    # plot(f, m, type = "b", ylab = "magnitude", xlab = "Frequency")
    # lines(filterFR$f / pi, abs(filterFR$h), col = "blue")
    # plot(f, 20*log10(m+1e-5), type = "b", ylab = "dB", xlab = "Frequency")
    # lines(filterFR$f / pi, 20*log10(abs(filterFR$h)), col = "blue")
    filterFR$f <- (samp.rate/2)*(filterFR$f / pi)
    filterFR$h <- abs(filterFR$h)
    # plot(1e-6*filterFR$f, filterFR$h, type = "l", ylab = "Magnitude / a.u.", xlab = "Frequency / MHz")
    filterFR <- approx(x = filterFR$f,
                       y = filterFR$h,
                       xout = FR$freq,
                       rule = 2) # add d.c. frequency and amplitude
    filterFR <- data.frame(freq = filterFR$x, amp = filterFR$y)
  }

  ## calculate readout frequency response
  if (!is.null(filterRange)) {
    readoutFR <- FR
    readoutFR[, 2] <- FR[, 2]*filterFR$amp*10^(amplifierGain/20)
  } else {
    readoutFR <- FR
    readoutFR[, 2] <- FR[, 2]*10^(amplifierGain/20)
  }

  ## plot readout FR if plot is TRUE
  if (plot) {
    plot(1e-6*readoutFR[, 1],
         20*log10(readoutFR[, 2]),
         type = "l",
         xlab = "Frequency [MHz]",
         ylab = "Magnitude [dB]")
    grid()
  }

  ## make output
  object <- list(sensorName = sensorName,
                 filterRange = filterRange,
                 amplifierGain = amplifierGain,
                 samp.rate = samp.rate,
                 readoutFR = readoutFR)
  class(object) <- "readout"
  return(object)
}
