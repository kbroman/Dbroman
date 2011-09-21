
import std.stdio;
import calc_squares_mod;

void main(string[] args)
{
  auto error_msg = "Give argument \"A\", \"A2\", \"A2p\", \"A2pp\", \"B\" or \"direct\".";
  assert(args.length > 1, error_msg);
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
  case "A2":
    writeln("   [Method A2]");
    R_calc_squares_A2(&n, &(y[0]), &(z[0])); break;
  case "A2p":
    writeln("   [Method A2p]");
    R_calc_squares_A2p(&n, &(y[0]), &(z[0])); break;
  case "A2pp":
    writeln("   [Method A2pp]");
    R_calc_squares_A2pp(&n, &(y[0]), &(z[0])); break;
  case "B":
    writeln("   [Method B]");
    R_calc_squares_B(&n, &(y[0]), &(z[0])); break;
  case "direct":
    writeln("   [Method direct]");
    R_calc_squares_direct(&n, &(y[0]), &(z[0])); break;
  default: 
    throw new Exception(error_msg);
  }
  writeln("   output = ", z);
}