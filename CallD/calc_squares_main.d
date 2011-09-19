
import std.stdio;
import calc_squares_mod;

void main()
{
  writeln(" --main");
  auto y = [5.3, 7.5, 3.2, 1, 2];
  auto z = y.dup;
  int n = y.length;
  writeln("    input = ", y);
  R_calc_squares(&n, &(y[0]), &(z[0]));
  writeln("   output = ", z);
}