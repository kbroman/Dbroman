Test for calling a D function from C

calc_square.R: R code
calc_square.d, calc_square_main.d: simple case (pass just one double)
calc_squares.d, calc_squares_main.d: pass a vector of doubles, 3 ways


From command-line, all methods work fine
$ calc_squares_main A
$ calc_squares_main B
$ calc_squares_main direct


From within R, methods B and direct work fine, but A gives a bus error

> source("calc_square.R")
> calc.squares(c(5.3, 7.5, 3.2, 1, 2), "A")

     *** caught bus error ***
    address 0x0, cause 'non-existent physical address'

    Traceback:
     1: .C(Rfunc, as.integer(length(x)), as.double(x), z = as.double(x))
     2: calc.squares(c(5.3, 7.5, 3.2, 1, 2), "A")
