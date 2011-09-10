module calc_squares_mod;

double[] calc_squares(in double[] x)
{
  auto y = x.dup;
  foreach(i, xv; x) {
    y[i] = x[i]*x[i];
  }
  return(y);
}

extern(C) void R_calc_squares(int n, double *x, double *output)
{
  double[] y = x[0..n];
  auto result = calc_squares(y);
  output[0..n] = result;
} 

