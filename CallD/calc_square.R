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
function(x)
{
  file <- "calc_squares.so"
  if(!is.loaded("R_calc_squares")) {
    dyn.load(file)
    cat(" -Loaded", file, "\n")
  }

 .C("R_calc_squares",
    as.integer(length(x)),
    as.double(x),
    z = as.double(x))$z
}

calc.squares.alt <-
function(x)
{
  file <- "calc_squares.so"
  if(!is.loaded("R_calc_squares")) {
    dyn.load(file)
    cat(" -Loaded", file, "\n")
  }

 .C("R_calc_squares_alt",
    as.integer(length(x)),
    as.double(x),
    z = as.double(x))$z
}
