all: calc_square.so calc_squares.so calc_square_main calc_squares_main

calc_square.so: calc_square.o
	gcc -m32 -shared calc_square.o -o calc_square.so -lphobos2

calc_square.o: calc_square.d
	dmd -fPIC -c calc_square.d -L-shared

calc_squares.so: calc_squares.o
	gcc -m32 -shared calc_squares.o -o calc_squares.so -lphobos2

calc_squares.o: calc_squares.d
	dmd -fPIC -c calc_squares.d -L-shared

calc_square_main: calc_square_main.d calc_square.d
	dmd -unittest calc_square_main.d calc_square.d -ofcalc_square_main 

calc_squares_main: calc_squares_main.d calc_squares.d
	dmd -unittest calc_squares_main.d calc_squares.d -ofcalc_squares_main

clean: 
	rm *.o *.so calc_square_main calc_squares_main