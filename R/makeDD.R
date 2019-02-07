#' @title Generating dispersion data.
#' @description A function that generates dispersion data.
#' @param phaseData.file A text file with exported dispersion data from Disperse , Default: "titanium.txt"
#' @param thickness The plate thickness in m, Default: 0.001
#' @param density The material density in kg/m3, Default: 4460
#' @param longVelocity The material longitudinal velocity in m/s, Default: 6060
#' @param tranVelocity The material transverse velocity in m/s, Default: 3230
#' @param plot If \code{TRUE}, plots the dispersion data, Default: FALSE
#' @return A dataset with dispersion data embedded into the \code{pongaR} library.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname makeDD
#' @export
makeDD <- function(phaseData.file = "titanium.txt",
                   thickness = 1e-3,
                   density = 4460,
                   longVelocity = 6060,
                   tranVelocity = 3230,
                   plot = FALSE) {

  ## sanity check
  if (phaseData.file == "titanium.txt") {
    try(phaseData.file <- system.file("extdata",
                                      "titanium.txt",
                                      package = "pongaR",
                                      mustWork = T),
        silent = T)
  }
  if (is.character(phaseData.file)) {
    if (file_ext(basename(phaseData.file)) != "txt") {
      stop("Phase data must be supplied as a text file (*.txt) exported from Disperse software.")
    }
    file <- file(phaseData.file, "r")
    on.exit(close(file))
  } else {
    stop("Argument 'phaseData.file' must be a character string containing file names or paths.")
  }
  if (!isOpen(file)) {
    open(file, "r")
    on.exit(close(file))
  }
  if (!is.numeric(thickness) | (thickness <= 0)) {
    stop("Argument 'thickness' must be a positive number.")
  }
  if (!is.numeric(density) | (density <= 0)) {
    stop("Argument 'density' must be a positive number.")
  }
  if (!is.numeric(longVelocity) | (longVelocity <= 0)) {
    stop("Argument 'longVelocity' must be a positive number.")
  }
  if (!is.numeric(tranVelocity) | (tranVelocity <= 0)) {
    stop("Argument 'tranVelocity' must be a positive number.")
  }
  if (!is.logical(plot)) {
    stop("Argument 'plot' must be an object of type logical")
  }

  ## read file
  df <- read.delim(file = file,
                   header = FALSE,
                   stringsAsFactors = FALSE,
                   skip = 2)
  nameIndex <- paste(c("S", "A"),
                     rep(0:(dim(df)[2]/(2*2) - 1), each = 2),
                     sep = "")

  ## get frequency and phase velocity from all modes
  freq <- df[, seq(1, dim(df)[2], 2)] # [MHz]
  v_ph <- df[, seq(2, dim(df)[2], 2)] # [m/ms]
  colnames(freq) <- freq[1, ]
  colnames(v_ph) <- freq[1, ]
  freq <- freq[-(1:2), nameIndex]
  v_ph <- v_ph[-(1:2), nameIndex]
  # freq <- freq[-(1:2), ]
  # v_ph <- v_ph[-(1:2), ]
  # freq <- as.numeric(as.matrix(freq))
  # v_ph <- as.numeric(as.matrix(v_ph))
  # freq <- matrix(data = freq, ncol = dim(df)[2]/2)*1e6 # [Hz]
  # v_ph <- matrix(data = v_ph, ncol = dim(df)[2]/2)*1e3 # [m/s]
  freq <- matrix(data = as.numeric(as.matrix(freq)),
                 ncol = ncol(freq),
                 dimnames = list(1:dim(freq)[1],
                                 colnames(freq)))*1e6 # [Hz]
  v_ph <- matrix(data = as.numeric(as.matrix(v_ph)),
                 ncol = ncol(freq),
                 dimnames = list(1:dim(v_ph)[1],
                                 colnames(v_ph)))*1e3 # [m/s]
  # dimnames(freq) <- list(1:dim(freq)[1],
  #                        nameIndex)
  # dimnames(v_ph) <- list(1:dim(v_ph)[1],
  #                        nameIndex)
  # dimnames(freq) <- list(1:dim(freq)[1],
  #                        colnames(freq))
  # dimnames(v_ph) <- list(1:dim(v_ph)[1],
  #                        colnames(v_ph))

  ## limit maximum phase velocity to 10.000 m/s (see Wilcox's PhD thesis p. 143, first paragraph)
  v_ph[v_ph > 10e3] <- NA

  ## make dimesionless
  freq <- freq*thickness/tranVelocity # [-]
  phasVelocity <- v_ph/tranVelocity # [-]

  ## make material dataset
  dataName <- paste0(strsplit(basename(phaseData.file), ".", fixed = T)[[1]][1],
                     "DD")
  filePath <- paste0("data/", dataName, ".rda")
  data <- list(thickness = thickness,
               density = density,
               longVelocity = longVelocity,
               tranVelocity = tranVelocity,
               freq = freq,
               phasVelocity = phasVelocity)
  assign(dataName, data)
  save(list = dataName, file = filePath)

  ## Wilcox's PhD thesis fig. 2.3, p. 47
  if (plot) {
      plot(freq[, 1],
           phasVelocity[, 1],
           type = "l",
           ylim = c(0, 5),
           xlab = "Frequency [-]",
           ylab = "Phase velocity [-]")
      for (k in 2:dim(freq)[2]) {
          lines(freq[, k], phasVelocity[, k])
      }
      grid()
  }
}
