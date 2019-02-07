#' @title Running experiment to produce acoustic emission signals due solid particle impacts onto a plate.
#' @description A function that runs experiment to produce acoustic emission signals
#'              due solid particle impacts onto a plate.
#' @param expName The experiment name, Default: 'exp1'
#' @param expData An object of class \code{exp} produced by \code{makeExp}, Default: makeExp()
#' @param nExp The number of experiments, Default: 2
#' @return A dataset with experiment run data embedded into the \code{pongaR} library.
# @details DETAILS
# @examples
# \dontrun{
# if(interactive()){
#  #EXAMPLE1
#  }
# }
#' @rdname runExp
#' @export
runExp <- function(expName = "exp1",
                   expData = makeExp(),
                   nExp = 2) {

  ## sanity check
  if (!is.character(expName)) {
    stop("Argument 'expName' must be a character string.")
  }
  if (class(expData) != "exp") {
    stop("Argument 'expData' must be an object of class exp")
  }
  if (!is.numeric(nExp) | (nExp <= 0)) {
    stop("Argument 'nExp' must be a positive number.")
  }

  aux <- list()
  aux$expData <- expData
  aux$runData <- list()
  n <- 1
  while (n <= nExp) {
    runData <- makeRun(expData = aux$expData)
    signal <- makeSignal(runData = runData)
    aux$runData[[n]] <- list()
    aux$runData[[n]]$expMatrix <- runData$expMatrix
    aux$runData[[n]]$out <- runData$out
    aux$runData[[n]]$signal <- signal
    n <- n + 1
  }
  assign(expName, aux)
  filePath <- paste("data/", expName, ".rda", sep = "")
  eval(parse(text = paste("save(", expName, ", file = filePath)", sep = "")))

}
