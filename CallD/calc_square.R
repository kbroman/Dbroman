calc.square <-
function(x)
{
  file <- "calc_square.so"
  if(!is.loaded("R_calc_square")) {
    dyn.load(file)
    cat(" -Loaded", file, "\n")
  }

 .C("R_calc_square",
    as.double(x),
    z = as.double(0))$z
}

calc.squares <-
function(x, method=c("A","A2","A2p", "A2pp", "B","direct"))
{
  method <- match.arg(method)
  Rfunc <- paste("R_calc_squares_", method, sep="")

  file <- "calc_squares.so"
  if(!is.loaded(Rfunc)) {
    dyn.load(file)
    cat(" -Loaded", file, "\n")
  }

 .C(Rfunc,
    as.integer(length(x)),
    as.double(x),
    z = as.double(x))$z
}
