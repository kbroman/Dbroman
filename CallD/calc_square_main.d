import std.stdio, std.conv, std.algorithm;

double calc_square(double x)
{
  return(x*x);
}

void R_calc_square(ref double x, ref double output)
{
  output = calc_square(x);
}

void main()
{
  auto y = 5.3;
  auto z = 2.0;
  R_calc_square(y, z);
  writeln(y, " ", z);
}