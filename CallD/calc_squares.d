module calc_squares_mod;
import std.stdio;

double[] calc_squares(in double[] x)
{
  auto y = x.dup;
  foreach(i, xv; x) {
    y[i] = xv*xv;
  }
  return(y);
}

unittest {
  writeln(" --Unit test calc_squares");
  assert(calc_squares([1.0, 2.0, 4.0]) == [1.0, 4.0, 16.0]);
}


extern(C) void R_calc_squares(int *n, double *x, double *output)
{
  double[] y = x[0..(*n)];
  output[0..(*n)] = calc_squares(y);
} 

