calc_square.so: calc_square.o
	gcc -m32 -shared calc_square.o -o calc_square.so

calc_square.o: calc_square.d
	dmd -c calc_square.d
