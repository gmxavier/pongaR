#' @title Generating sensor frequency response data.
#' @description A function that generates sensor frequency response data.
#' @param sensorCalib The .csv file with the digitized  calibration  curve, Default: 'r15iuc.csv'
#' @param xunit Character indicating which unit is used in for the frequency axis, Default: c("Hz", "kHz", "MHz")
#' @param yunit Character indicating which unit is used in for the magnitude axis, Default: c("V/ubar", "V/(m/s)", "V/m")
#' @param dBref A dB reference value, Default: 1
#' @param n Number of frequency response points, Default: 4097
#' @param plot If \code{TRUE}, plots the sensor frequency response data, Default: FALSE
#' @return A dataset with sensor frequency response data embedded into the \code{pongaR} library.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname makeFR
#' @export
makeFR <- function(sensorCalib = "r15iuc.csv",
                   xunit = c("Hz", "kHz", "MHz"),
                   yunit = c("V/ubar", "V/(m/s)", "V/m"),
                   dBref = 1,
                   n = 4097,
                   plot = FALSE) {
  ## sanity check
  if (sensorCalib == "r15iuc.csv") {
    try(sensorCalib <- system.file("extdata",
                                   "r15iuc.csv",
                                   package = "pongaR",
                                   mustWork = T),
        silent = T)
  }
  if (is.character(sensorCalib)) {
    if (rev(strsplit(sensorCalib, "([/.])")[[1]])[1] != "csv") {
      stop("Sensor calibration data must be supplied as a comma separated value file (*.csv).")
    }
    file <- file(sensorCalib, "r")
    on.exit(close(file))
  } else {
    stop("Argument 'sensorCalib' must be a character string containing file names or paths.")
  }
  if (!isOpen(file)) {
    open(file, "r")
    on.exit(close(file))
  }
  xunit <- match.arg(xunit)
  if (xunit != "Hz" & xunit != "kHz" & xunit != "MHz") {
    stop("Argument 'xunit' must be 'Hz', 'kHz' or 'MHz'.")
  }
  yunit <- match.arg(yunit)
  if (yunit != "V/ubar" & yunit != "V/(m/s)" & yunit != "V/m") {
    stop("Argument 'yunit' must be 'V/ubar', 'V/(m/s)' or 'V/m'.")
  }
  if (!is.numeric(n) | (n < 0)) {
      stop("Argument 'dBref' must be equal or greater than 0.")
  }
  if (!is.numeric(n) | (n < 4097)) {
    stop("Argument 'n' must be equal or greater than 4097.")
  }
  if (!is.logical(plot)) {
    stop("Argument 'plot' must be an object of type logical")
  }

  ## read file
  df <- read.csv(file = file,
                 header = FALSE,
                 comment.char = "#",
                 stringsAsFactors = FALSE,
                 skip = 6)
  df <- df[, -3]

  ## change frequency to Hz
  if (xunit == "kHz") {
    df$V1 <- df$V1*1e3
  } else if (xunit == "MHz") {
    df$V1 <- df$V1*1e6
  }

  ## change amplitude scale from log to linear
  if (dBref != 0) {
    df$V2 <- 10^(df$V2/20)*dBref
  }

  ## compute sensor frequency response (frequency in Hz, amplitude in dB, dBref = 1 V/ubar)
  upperfreq <- df$V1[length(df$V1)]
  df$V1[length(df$V1)] <- round(upperfreq/10^floor(log10(upperfreq)),
                                digits = 1)*10^floor(log10(upperfreq)) # round last frequency
  FR <- approx(x = c(0, df$V1),
               y = c(df$V2[1], df$V2),
               n = n,
               rule = 2) # add d.c. frequency and amplitude
  FR <- data.frame(freq = FR$x, amp = FR$y)
  names(FR) <- c("freq_Hz", paste("amp_", yunit, sep = ""))

  ## plot sensor FR if plot is TRUE
  if (plot) {
      plot(1e-6*FR$freq,
           20*log10(FR$amp),
           type = "l",
           xlab = "Frequency [MHz]",
           ylab = ifelse(dBref,
                         paste("Magnitude [dB ref ", dBref, yunit, "]", sep = ""),
                         paste("Magnitude [", yunit, "]", sep = "")))
      grid()
      plot(2*pi*FR$freq,
           FR$amp,
           type = "l",
           xlab = expression(paste(omega, " [rad ", "s"^{-1}, "]")),
           ylab = expression(paste("G"["C"],"(", omega, ") [V ", mu, "bar"^{-1}, "]")))
      grid()
  }

  ## make sensor FR dataset
  dataName <- strsplit(basename(sensorCalib), split = ".csv", fixed = T)[[1]][1]
  filePath <- paste0("data/", dataName, ".rda")
  assign(dataName, FR)
  save(list = dataName, file = filePath)

  # filePath <- paste("data/", rev(strsplit(sensorCalib, "([/.])")[[1]])[2], ".rda", sep = "")
  # save(FR, file = filePath)
}
