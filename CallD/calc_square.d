module calc_square_mod;
import std.stdio;

double calc_square(in double x)
{
  return(x*x);
}

unittest {
  writeln(" --Unit test calc_square");
  assert(calc_square(5.0) == 25.0);
}

extern(C) void R_calc_square(ref double x, ref double output)
{
  output = calc_square(x);
} 

