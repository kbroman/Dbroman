module calc_squares_mod;
import std.stdio;

double[] calc_squares(in double[] x)
{
  auto y = x.dup;
  foreach(i, xv; x) {
    y[i] = x[i]*x[i];
  }
  return(y);
}

extern(C) void R_calc_squares(int *n, double *x, double *output)
{
  double[] y = x[0..(*n)];
  output[0..(*n)] = calc_squares(y);
} 

