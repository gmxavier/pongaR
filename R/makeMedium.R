#' @title Generating medium impulse response.
#' @description A function that generates medium impulse response.
#' @param plate Plate specification, Default: makePlate()
#' @param plateTheory Plate theory model, Default: c("lamb", "mindlin", "kirchhoff")
#' @param excitabilityData Excitability function data, Default: 'titaniumEF'
#' @param modes Maximum number of modes \code{c(m,n)}, Default: NULL
#' @param sensorDiameter Sensor diameter (in m), Default: 0
#' @param sensorPosition Sensor position coordinates \code{c(x,y)} in m, Default: c(0.254, 0.19)
#' @param sourcePosition Sensor position coordinates \code{c(x0,y0)} in m, Default: c(0.204, 0.19)
#' @param sourceType Source type, Default: c("impulse", "step")
#' @param duration Duration of the \code{Wave} in \code{xunit}, Default: 0.0001024
#' @param xunit Character indicating which unit is used in argument \code{duration}.
#'              If \code{xunit = "time"}, the unit is time in s,
#'              otherwise the number of samples, Default: c("time", "samples")
#' @param samp.rate Sampling rate of the \code{Wave}, in Hz, Default: 1e+07
#' @param bit Resolution of the \code{Wave} and rescaling unit, Default: 32. This may be
#'            1 for rescaling to numeric values in [-1,1],
#'            8 (i.e. 8-bit) for rescaling to integers in [0, 254],
#'            16 (i.e. 16-bit) for rescaling to integers in [-32767, 32767],
#'            24 (i.e. 24-bit) for rescaling to integers in [-8388607, 8388607],
#'            32 (i.e. 32-bit) for rescaling either to integers in [-2147483647, 2147483647],
#'            64 (i.e. 64-bit) for rescaling to numeric values in [-1, 1] (FLOAT_IEEE Wave format), and
#'            0 for not rescaling at all.
#' @return The medium impulse response waveform as a \code{medium} object.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname makeMedium
#' @export
makeMedium <- function(plate = makePlate(),
                       plateTheory = c("lamb", "mindlin", "kirchhoff"),
                       excitabilityData = "titaniumEF",
                       modes = NULL,
                       sensorDiameter = 0,
                       sensorPosition = c(0.254, 0.190),
                       sourcePosition = c(0.204, 0.190),
                       sourceType = c("impulse", "step"),
                       duration = 102.4e-6,
                       xunit = c("time", "samples"),
                       samp.rate = 10e6,
                       bit = 32) {

    ## this function must reproduce fig. 3.11 (Prosser's PhD thesis)
    # test_makeMedium <- function() {
    #   mediumData <- makeMedium(sensorPosition = c(0.254, 0.254),
    #                            sourceType = "step",
    #                            plateTheory = "mindlin",
    #                            modes = c(300, 300),
    #                            duration = 204.8e-6,
    #                            samp.rate = 5e6)
    #   time <- seq(from = 1/5e6,
    #               to = 204.8e-6,
    #               by = 1/5e6)
    #   plot(1e6*time,
    #        mediumData$waveform@left,
    #        type = "l",
    #        ylim = c(-2e-9, 2e-9),
    #        c(0,100),
    #        ylab = "Amplitude [m]",
    #        xlab = expression(paste("Time [", mu, "s]")))
    #   grid()
    # }

    longVelocity <- tranVelocity <- EF_c <- freq <- wavenumber <- NULL

    ## sanity check
    if ((class(plate) != "plate")) {
        stop("Argument 'plate' must be an object of class plate.")
    }
    plateTheory <- match.arg(plateTheory)
    if (plateTheory != "mindlin" & plateTheory != "kirchhoff" & plateTheory != "lamb") {
        stop("Argument 'plateTheory' must be 'mindlin', 'kirchhoff' or 'lamb'.")
    }
    if (plateTheory == "lamb") {
        if (is.character(excitabilityData)) {
            data(list = excitabilityData, envir = environment())
        } else {
            stop("Argument 'excitabilityData' must be a character string.")
        }
    } else {
        if (!is.vector(x = modes, mode = "numeric") | any(modes <= 0)) {
            stop("Argument 'modes' must be a positive numeric vector.")
        }
    }
    if (!is.numeric(sensorDiameter) | (sensorDiameter < 0)) {
        stop("Argument 'sensorDiameter' must be a non negative number.")
    }
    if (!is.vector(x = sensorPosition, mode = "numeric") | any(sensorPosition < 0)) {
        stop("Argument 'sensorPosition' must be a non negative numeric vector.")
    } else {
        if (any((c(plate$length, plate$width) >= sensorPosition) != TRUE)) {
            stop("Receiver is out of plate. Check receiver position coordinates.")
        }
        if (!is.null(sensorDiameter) & any((sensorPosition >= sensorDiameter/2) != TRUE)) {
            stop("Receiver is out of plate. Check receiver diameter.")
        }
    }
    if (!is.vector(x = sourcePosition, mode = "numeric") | any(sourcePosition < 0)) {
        stop("Argument 'sourcePosition' must be a non negative numeric vector.")
    } else {
        if (any((c(plate$length, plate$width) >= sourcePosition) != TRUE)) {
            stop("Source is out of plate. Check source position coordinates.")
        }
    }
    sourceType <- match.arg(sourceType)
    if (sourceType != "step" & sourceType != "impulse") {
        stop("Argument 'sourceType' must be 'step' or 'impulse'.")
    }
    if (!is.numeric(duration) | (duration <= 0)) {
        stop("Argument 'duration' must be a positive number.")
    }
    xunit <- match.arg(xunit)
    if (xunit != "samples" & xunit != "time") {
        stop("Argument 'xunit' must be 'samples' or 'time'.")
    }
    if (!is.numeric(samp.rate) | (samp.rate <= 0)) {
        stop("Argument 'samp.rate' must be a positive number.")
    }

    ## read parameters
    m <- modes[1]
    n <- modes[2]
    qsi <- sourcePosition[1]
    zeta <- sourcePosition[2]
    x0 <- sensorPosition[1]
    y0 <- sensorPosition[2]

    ## compute response
    if (plateTheory == "mindlin") {
        ## compute parameters
        alpha_n <- seq(1:n)*pi/plate$length
        gamma_m <- seq(1:m)*pi/plate$width
        E2_mn <- matrix(data = plate$transverseShearmodulus44*gamma_m^2,
                        nrow = m,
                        ncol = n)
        E2_mn <- E2_mn + matrix(data = plate$transverseShearmodulus55*alpha_n^2,
                                nrow = m,
                                ncol = n,
                                byrow = TRUE)
        E2_mn <- E2_mn/(plate$density*plate$thickness)
        omega2_mn <- matrix(data = plate$flexuralRigidity22*gamma_m^4,
                            nrow = m,
                            ncol = n)
        omega2_mn <- omega2_mn + matrix(data = plate$flexuralRigidity11*alpha_n^4,
                                        nrow = m,
                                        ncol = n,
                                        byrow = TRUE)
        omega2_mn <- omega2_mn + 2*(plate$flexuralRigidity12 + 2*plate$flexuralRigidity66)*gamma_m^2%*%t(alpha_n^2)
        omega2_mn <- omega2_mn/(plate$density*plate$thickness)
        beta2_mn <- omega2_mn/(1 + omega2_mn/E2_mn)
        beta_mn <- sqrt(beta2_mn)
        # compute aperture effect
        if (sensorDiameter == 0) {
            aux <- sin(gamma_m*y0)%*%t(sin(alpha_n*x0)) # response without aperture effect
        } else {
            aux <- matrix(data = 0, nrow = m, ncol = n)
            lower <- x0 - (sensorDiameter/2)
            upper <- x0 + (sensorDiameter/2)
            receiverRadius2 <- (sensorDiameter/2)^2
            for (i in 1:m) {
                for (j in 1:n) {
                    aux[i, j] <- (1/gamma_m[i])*integrate(apertureFunction,
                                                          lower = lower,
                                                          upper = upper,
                                                          gamma_m[i],
                                                          alpha_n[j],
                                                          x0,
                                                          y0,
                                                          receiverRadius2,
                                                          stop.on.error = FALSE)$value
                }
            }
            aux <- (1/pi/(sensorDiameter/2)^2)*aux # response with aperture effect
        }
        aux <- aux*sin(gamma_m*zeta)%*%t(sin(alpha_n*qsi))
        aux <- (4/(plate$density*plate$thickness*plate$length*plate$width))*aux/beta2_mn
        ## compute displacement
        # sample <- 1
        if (xunit == "samples") {
            time <- seq(from = 0,
                        to = (duration - 1),
                        by = 1)/samp.rate
            displacement <- rep(0, duration)
        } else {
            time <- seq(from = 0,
                        to = (duration*samp.rate - 1),
                        by = 1)/samp.rate
            displacement <- rep(0, duration*samp.rate)
        }
        if (sourceType == "step") {
            aux1 <- plate$lossFactor*beta_mn
            aux2 <- beta_mn*sqrt(1 - plate$lossFactor^2)
            theta <- atan(plate$lossFactor/sqrt(1 - plate$lossFactor^2))
            # while (sample <= length(displacement)) {
            #     displacement[sample] <- sum(aux - (aux/sqrt(1 - plate$lossFactor^2)*exp(-aux1*time[sample])*cos(aux2*time[sample] - theta)))
            #     sample <- sample + 1
            # }
            cl <- makeCluster(getOption("cl.cores", (detectCores() - 2))) # cluster startup
            displacement <- parSapply(cl,
                                      time,
                                      function(time, aux, plate, aux1, aux2, theta) sum(aux - (aux/sqrt(1 - plate$lossFactor^2)*exp(-aux1*time)*cos(aux2*time - theta))),
                                      aux = aux,
                                      plate = plate,
                                      aux1 = aux1,
                                      aux2 = aux2,
                                      theta = theta) # parallel calculation
            stopCluster(cl) # cluster shutdown
            displacement <- Wave(left = displacement, samp.rate = samp.rate, bit = bit)
            object <- list(plate = plate,
                           excitabilityData = excitabilityData,
                           sensorDiameter = sensorDiameter,
                           sensorPosition = sensorPosition,
                           sourcePosition = sourcePosition,
                           sourceType = sourceType,
                           modes = modes,
                           plateTheory = plateTheory,
                           duration = duration,
                           xunit = xunit,
                           samp.rate = samp.rate,
                           bit = bit,
                           frequencyRange = c(min(beta_mn), max(beta_mn)/2/pi),
                           waveform = displacement)
            class(object) <- "medium"
            return(object)
        }
        if (sourceType == "impulse") {
            aux <- aux/sqrt(1 - plate$lossFactor^2)
            aux1 <- plate$lossFactor*beta_mn
            aux2 <- beta_mn*sqrt(1 - plate$lossFactor^2)
            # while (sample <= length(displacement)) {
            #     displacement[sample] <- sum(aux*(beta_mn*exp(-aux1*time[sample])*sin(aux2*time[sample])))
            #     sample <- sample + 1
            # }
            # displacement <- sapply(time,
            #                        function(time, aux, beta_mn, aux1, aux2) sum(aux*(beta_mn*exp(-aux1*time)*sin(aux2*time))),
            #                        aux = aux,
            #                        beta_mn = beta_mn,
            #                        aux1 = aux1,
            #                        aux2 = aux2)
            cl <- makeCluster(getOption("cl.cores", (detectCores() - 1))) # cluster startup
            cl <- makeForkCluster(nnodes = getOption("mc.cores", (detectCores() - 1)))
            displacement <- parSapply(cl,
                                      time,
                                      function(time, aux, beta_mn, aux1, aux2) sum(aux*(beta_mn*exp(-aux1*time)*sin(aux2*time))),
                                      aux = aux,
                                      beta_mn = beta_mn,
                                      aux1 = aux1,
                                      aux2 = aux2) # parallel calculation
            stopCluster(cl) # cluster shutdown
            displacement <- Wave(left = displacement, samp.rate = samp.rate, bit = bit)
            object <- list(plate = plate,
                           excitabilityData = excitabilityData,
                           sensorDiameter = sensorDiameter,
                           sensorPosition = sensorPosition,
                           sourcePosition = sourcePosition,
                           sourceType = sourceType,
                           modes = modes,
                           plateTheory = plateTheory,
                           duration = duration,
                           xunit = xunit,
                           samp.rate = samp.rate,
                           bit = bit,
                           frequencyRange = c(min(beta_mn), max(beta_mn)/2/pi),
                           waveform = displacement)
            class(object) <- "medium"
            return(object)
        }
    }
    if (plateTheory == "kirchhoff") {
        ## compute parameters
        alpha_n <- seq(1:n)*pi/plate$length
        gamma_m <- seq(1:m)*pi/plate$width
        omega2_mn <- matrix(data = plate$flexuralRigidity22*gamma_m^4,
                            nrow = m,
                            ncol = n)
        omega2_mn <- omega2_mn + matrix(data = plate$flexuralRigidity11*alpha_n^4,
                                        nrow = m,
                                        ncol = n,
                                        byrow = TRUE)
        omega2_mn <- omega2_mn + 2*(plate$flexuralRigidity12 + 2*plate$flexuralRigidity66)*gamma_m^2%*%t(alpha_n^2)
        omega2_mn <- omega2_mn/(plate$density*plate$thickness)
        beta2_mn <- omega2_mn
        beta_mn <- sqrt(beta2_mn)
        # compute aperture effect
        if (sensorDiameter == 0) {
            aux <- sin(gamma_m*y0)%*%t(sin(alpha_n*x0)) # response without aperture effect
        } else {
            aux <- matrix(data = 0, nrow = m, ncol = n)
            lower <- x0 - (sensorDiameter/2)
            upper <- x0 + (sensorDiameter/2)
            receiverRadius2 <- (sensorDiameter/2)^2
            for (i in 1:m) {
                for (j in 1:n) {
                    aux[i, j] <- (1/gamma_m[i])*integrate(apertureFunction,
                                                          lower = lower,
                                                          upper = upper,
                                                          gamma_m[i],
                                                          alpha_n[j],
                                                          x0,
                                                          y0,
                                                          receiverRadius2,
                                                          stop.on.error = FALSE)$value
                }
            }
            aux <- (1/pi/(sensorDiameter/2)^2)*aux # response with aperture effect
        }
        aux <- aux*sin(gamma_m*zeta)%*%t(sin(alpha_n*qsi))
        aux <- (4/(plate$density*plate$thickness*plate$length*plate$width))*aux/beta2_mn

        ## compute displacement
        if (xunit == "samples") {
            time <- seq(from = 0,
                        to = (duration - 1),
                        by = 1)/samp.rate
            displacement <- rep(0, duration)
        } else {
            time <- seq(from = 0,
                        to = (duration*samp.rate - 1),
                        by = 1)/samp.rate
            displacement <- rep(0, duration*samp.rate)
        }
        if (sourceType == "step") {
            aux1 <- plate$lossFactor*beta_mn
            aux2 <- beta_mn*sqrt(1 - plate$lossFactor^2)
            theta <- atan(plate$lossFactor/sqrt(1 - plate$lossFactor^2))
            # sample <- 1
            # while (sample <= length(displacement)) {
            #     displacement[sample] <- sum(aux - (aux/sqrt(1 - plate$lossFactor^2)*exp(-aux1*time[sample])*cos(aux2*time[sample] - theta)))
            #     sample <- sample + 1
            # }
            cl <- makeCluster(getOption("cl.cores", (detectCores() - 1))) # cluster startup
            displacement <- parSapply(cl,
                                      time,
                                      function(time, aux, plate, aux1, aux2, theta) sum(aux - (aux/sqrt(1 - plate$lossFactor^2)*exp(-aux1*time)*cos(aux2*time - theta))),
                                      aux = aux,
                                      plate = plate,
                                      aux1 = aux1,
                                      aux2 = aux2,
                                      theta = theta) # parallel calculation
            stopCluster(cl) # cluster shutdown
            displacement <- Wave(left = displacement, samp.rate = samp.rate, bit = bit)
            object <- list(plate = plate,
                           excitabilityData = excitabilityData,
                           sensorDiameter = sensorDiameter,
                           sensorPosition = sensorPosition,
                           sourcePosition = sourcePosition,
                           sourceType = sourceType,
                           modes = modes,
                           plateTheory = plateTheory,
                           duration = duration,
                           xunit = xunit,
                           samp.rate = samp.rate,
                           bit = bit,
                           frequencyRange = c(min(beta_mn), max(beta_mn)/2/pi),
                           waveform = displacement)
            class(object) <- "medium"
            return(object)
        }
        if (sourceType == "impulse") {
            aux <- aux/sqrt(1 - plate$lossFactor^2)
            aux1 <- plate$lossFactor*beta_mn
            aux2 <- beta_mn*sqrt(1 - plate$lossFactor^2)
            # while (sample <= length(displacement)) {
            #     displacement[sample] <- sum(aux*(beta_mn*exp(-aux1*time[sample])*sin(aux2*time[sample])))
            #     sample <- sample + 1
            # }
            cl <- makeCluster(getOption("cl.cores", (detectCores() - 1))) # cluster startup
            displacement <- parSapply(cl,
                                      time,
                                      function(time, aux, beta_mn, aux1, aux2) sum(aux*(beta_mn*exp(-aux1*time)*sin(aux2*time))),
                                      aux = aux,
                                      beta_mn = beta_mn,
                                      aux1 = aux1,
                                      aux2 = aux2) # parallel calculation
            stopCluster(cl) # cluster shutdown
            displacement <- Wave(left = displacement, samp.rate = samp.rate, bit = bit)
            object <- list(plate = plate,
                           excitabilityData = excitabilityData,
                           sensorDiameter = sensorDiameter,
                           sensorPosition = sensorPosition,
                           sourcePosition = sourcePosition,
                           sourceType = sourceType,
                           modes = modes,
                           plateTheory = plateTheory,
                           duration = duration,
                           xunit = xunit,
                           samp.rate = samp.rate,
                           bit = bit,
                           frequencyRange = c(min(beta_mn), max(beta_mn))/(2*pi),
                           waveform = displacement)
            class(object) <- "medium"
            return(object)
        }
    }
    if (plateTheory == "lamb") {
        if (sourceType == "step") {
            stop("This source type is not available for Lamb theory.")
        }
        if (sourceType == "impulse") {

            density <- eval(parse(text = paste0(excitabilityData, "$", "density")))
            longVelocity <- eval(parse(text = paste0(excitabilityData, "$", "longVelocity")))
            tranVelocity <- eval(parse(text = paste0(excitabilityData, "$", "tranVelocity")))
            freq <- eval(parse(text = paste0(excitabilityData, "$", "freq")))
            EF_c <- eval(parse(text = paste0(excitabilityData, "$", "EF_c")))

            ## compute parameters
            thickness <- plate$thickness # [m]
            distance <- sqrt(sum((sensorPosition - sourcePosition)^2)) # [m]
            distance <- distance/thickness # [dimensionless]
            lameMu <- vel2lame(density = density, longVelocity = longVelocity, tranVelocity = tranVelocity)[[2]]
            ## make time vector [dimensionless]
            timeLength <- duration*samp.rate
            time <- seq(from = 0,
                        to = (timeLength - 1),
                        by = 1)/samp.rate
            time <- time*tranVelocity/thickness
            ## make input signal FT (input signal is a Dirac's delta)
            inputFT <- seq(0, samp.rate/2, length.out = timeLength)*thickness/tranVelocity #[-]
            inputFT <- cbind(inputFT, 1+0*1i)
            ## assign excitability function [dimensionless]
            EF <- EF_c
            ## plate TF collector for each mode
            plateTF <- matrix(data = NA, nrow = NROW(inputFT), ncol = NCOL(freq))
            ## Hankel transform collector for each mode
            Hfirst0 <- matrix(data = NA, nrow = NROW(inputFT), ncol = NCOL(freq))
            ## compute plate TF for each mode
            for (mode in seq(1, NCOL(freq))) {
                aux1 <- approx(x = freq[, mode], y = EF[, mode], xout = inputFT[, 1])
                aux2 <- approx(x = freq[, mode], y = wavenumber[, mode], xout = inputFT[, 1])
                Hfirst0[!is.na(aux2$y), mode] <- BesselJ(aux2$y[!is.na(aux2$y)]*distance, 0) + 1i*BesselY(aux2$y[!is.na(aux2$y)]*distance, 0)
                plateTF[, mode] <- aux1$y*Hfirst0[, mode]
            }
            ## compute final plate TF (sum of all modes)
            plateTF <- matrix(data = rowSums(plateTF[, ], na.rm = TRUE), ncol = 1)
            ## displacement per force collector
            displacement <- rep(0, timeLength)
            ## calculate dimensionless displacement per force by inverse FT due to an unit impulse dimensionless force
            for (n in 1:NROW(plateTF)) {
                u_z <- -inputFT[n, 2]*plateTF[n, 1]*exp(-1i*2*pi*inputFT[n, 1]*time)
                displacement <- displacement + u_z
            }
            ## dimensionless displacement per force corrections
            # displacement <- displacement - displacement[1] # remove offset
            displacement <- displacement*(1/(timeLength/2)) # due inverse FT
            ## displacement per force in m/N ( displacement*thickness/(lameMu*thickness^2) )
            displacement <- displacement/(lameMu*thickness)
            displacement <- Wave(left = Re(displacement),
                                 samp.rate = samp.rate,
                                 bit = bit)
            object <- list(plate = plate,
                           excitabilityData = excitabilityData,
                           sensorDiameter = sensorDiameter,
                           sensorPosition = sensorPosition,
                           sourcePosition = sourcePosition,
                           sourceType = sourceType,
                           modes = modes,
                           plateTheory = plateTheory,
                           duration = duration,
                           xunit = xunit,
                           samp.rate = samp.rate,
                           bit = bit,
                           frequencyRange = c(0, samp.rate/2),
                           waveform = displacement)
            class(object) <- "medium"
            return(object)
        }
    }
}
