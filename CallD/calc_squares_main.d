
import std.stdio;
import calc_squares_mod;

void main(string[] args)
{
  assert(args.length > 1, "Give argument \"A\", \"B\" or \"direct\".");
  auto method = args[1];

  writeln(" --main");

  auto y = [5.3, 7.5, 3.2, 1, 2];
  auto z = y.dup;
  int n = y.length;
  writeln("    input = ", y);

  switch(method) {
  case "A":
    writeln("   [Method A]");
    R_calc_squares_A(&n, &(y[0]), &(z[0])); break;
  case "B":
    writeln("   [Method B]");
    R_calc_squares_B(&n, &(y[0]), &(z[0])); break;
  case "direct":
    writeln("   [Method direct]");
    R_calc_squares_direct(&n, &(y[0]), &(z[0])); break;
  default: 
    throw new Exception("Give argument \"A\", \"B\" or \"direct\".");
  }
  writeln("   output = ", z);
}