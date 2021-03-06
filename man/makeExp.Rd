% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeExp.R
\name{makeExp}
\alias{makeExp}
\title{Generating experiment data.}
\usage{
makeExp(impactFrequency = 2000, expDuration = 0.001,
  sphereMaterial = "glassMcLaskey",
  diameterDist = "rnorm(n, mean = 1e-3, sd = 10e-6)",
  velocityDist = "rnorm(n, mean = 1, sd = 1e-3)",
  plateMaterial = "titaniumDisperse", plateGeometry = c(0.508, 0.381,
  0.005), lossFactor = 0, shearCorrectionfactor = NA,
  plateTheory = "lamb", excitabilityData = "titaniumEF",
  modes = NULL, sensorDiameter = 0, sensorPosition = c(round(0.508/2,
  3), round(0.381/2, 3)), sourceType = "impulse",
  sourceModel = "hertz", riseTime = 1.5e-06, forceAmplitude = 1,
  centerFreq = 2e+05, cyclesNum = 5, sensorName = "r15iuc",
  filterRange = c(50000, 1e+06), amplifierGain = 80, xunit = "time",
  samp.rate = 1e+07, bit = 32)
}
\arguments{
\item{impactFrequency}{Impact frequency, in impacts/s, Default: 2000}

\item{expDuration}{Experiment duration, in s, Default: 0.001}

\item{sphereMaterial}{Solid sphere material name, Default: 'glassMcLaskey'}

\item{diameterDist}{Distribution of the solid sphere diameter, in m, Default: 'rnorm(n, mean = 1e-3, sd = 10e-6)'}

\item{velocityDist}{Distribution of the impact velocity, in m/s, Default: 'rnorm(n, mean = 1, sd = 1e-3)'}

\item{plateMaterial}{Plate material name, Default: 'titaniumDisperse'}

\item{plateGeometry}{Plate geometry, in m, Default: c(0.508, 0.381, 0.005)}

\item{lossFactor}{Material loss factor, Default: 0}

\item{shearCorrectionfactor}{Material shear correction factor, Default: NA}

\item{plateTheory}{Plate theory model name, Default: 'lamb'}

\item{excitabilityData}{Excitability data name, Default: 'titaniumEF'}

\item{modes}{Maximum number of modes \code{c(m,n)}, Default: NULL}

\item{sensorDiameter}{Sensor diameter (in m), Default: 0}

\item{sensorPosition}{Sensor position coordinates \code{c(x,y)} in m, Default: c(round(0.508/2, 3), round(0.381/2, 3))}

\item{sourceType}{Source type, Default: 'impulse'}

\item{sourceModel}{The source model, Default: 'hertz'}

\item{riseTime}{The plb rise time, in s, Default: 1.5e-06}

\item{forceAmplitude}{The plb and tone burst force amplitude, in N, Default: 1}

\item{centerFreq}{The tone burst center frequency, in Hz, Default: 2e+05}

\item{cyclesNum}{The tone burst number of cycles, Default: 5}

\item{sensorName}{The sensor dataset name, Default: 'r15iuc'}

\item{filterRange}{The band-pass filter range, in Hz, Default: c(50000, 1e+06)}

\item{amplifierGain}{The amplfier gain in dB, Default: 80}

\item{xunit}{Character indicating which unit is used in argument \code{duration}.
If \code{xunit = "time"}, the unit is time in s,
otherwise the number of samples, Default: 'time'}

\item{samp.rate}{Sampling rate of the \code{Wave}, in Hz, Default: 1e+07}

\item{bit}{Resolution of the \code{Wave} and rescaling unit, Default: 32. This may be
1 for rescaling to numeric values in [-1,1],
8 (i.e. 8-bit) for rescaling to integers in [0, 254],
16 (i.e. 16-bit) for rescaling to integers in [-32767, 32767],
24 (i.e. 24-bit) for rescaling to integers in [-8388607, 8388607],
32 (i.e. 32-bit) for rescaling either to integers in [-2147483647, 2147483647],
64 (i.e. 64-bit) for rescaling to numeric values in [-1, 1] (FLOAT_IEEE Wave format), and
0 for not rescaling at all.}
}
\value{
The experimental data as a \code{exp} object.
}
\description{
A function that generates experiment data to produce several acoustic emission signal
             waveforms.
}
