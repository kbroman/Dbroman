module calc_squares_mod;
import std.stdio;


extern(C) void Rprintf(const char *, ...);

void calc_squares(in double[] x, double[] y)
{
  Rprintf("x.length: %d\n", x.length);
  Rprintf("y.length: %d\n", y.length);
  foreach(i, xv; x) {
    y[i] = xv*xv;
  }
}

unittest {
  writeln(" --Unit test calc_squares");
  auto y = [1.0, 2.0, 4.0];
  auto z = new double[](y.length);
  calc_squares(y, z);
  assert(z == [1.0, 4.0, 16.0]);
}


extern(C) void R_calc_squares(int *n, double *x, double *output)
{
  double[] y = x[0..(*n)];
  double[] z = output[0..(*n)];
  calc_squares(y, z);
} 

extern(C) void R_calc_squares_alt(int *n, double *x, double *output)
{
  foreach(i; 0..(*n)) 
    output[i] = x[i]*x[i];
} 

