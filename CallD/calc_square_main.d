import std.stdio;
import calc_square_mod;


void main()
{
  writeln(" --main");
  auto y = 5.3;
  auto z = 2.0;
  R_calc_square(y, z);
  writeln("    input = ", y);
  writeln("   output = ", z);
}