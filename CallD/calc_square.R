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
