module calc_squares_mod;
import std.stdio;

double[] calc_squares_A(in double[] x)
{
  auto y = x.dup;
  foreach(i, xv; x) {
    y[i] = xv*xv;
  }
  return(y);
}



void calc_squares_B(in double[] x, double[] y)
{
  foreach(i, xv; x) {
    y[i] = xv*xv;
  }
}

unittest {
  writeln(" --Unit test calc_squares_A");
  assert(calc_squares_A([1.0, 2.0, 4.0]) == [1.0, 4.0, 16.0]);

  writeln(" --Unit test calc_squares_B");
  auto y = [1.0, 2.0, 4.0];
  auto z = new double[](y.length);
  calc_squares_B(y, z);
  assert(z == [1.0, 4.0, 16.0]);
}


extern(C) void R_calc_squares_A(int *n, double *x, double *output)
{
  double[] y = x[0..(*n)];
  output[0..(*n)] = calc_squares_A(y);
} 

extern(C) void R_calc_squares_B(int *n, double *x, double *output)
{
  double[] y = x[0..(*n)];
  double[] z = output[0..(*n)];
  calc_squares_B(y, z);
} 

extern(C) void R_calc_squares_direct(int *n, double *x, double *output)
{
  foreach(i; 0..(*n)) 
    output[i] = x[i]*x[i];
} 

