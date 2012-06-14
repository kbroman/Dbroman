/**
 * Playing with classes
 **/

module classes;

import std.math, std.stdio, std.path;
import std.exception;
import std.conv;

class Cross {
  string cross_type;

  abstract double myfuncA();
  abstract double myfuncB();
}

class BC : Cross {
  this() { cross_type = "BC"; }

  double myfuncA() { return 1.0; }
  double myfuncB() { return myfuncA; }
}

class F2 : Cross {
  this() { cross_type = "F2"; }

  double myfuncA() { return 2.0; }
  double myfuncB() { return 3.0; }
}


Cross form_cross(string cross_type) {
  Cross cross;
  switch(cross_type) {
  case "BC": cross = new BC; break;
  case "F2": cross = new F2; break;
  default: throw new Exception("Cross type " ~ cross_type ~ " not supported.");
  }
  return cross;
}


double anotherfunc(Cross cross, string which)
{
  switch(which) {
  case "A": return cross.myfuncA();
  case "B": return cross.myfuncB();
  default: throw new Exception("which = \"" ~ which ~ "\" not supported.");
  }
}


unittest {
  writeln("Unit test " ~ __FILE__);

  auto bc = form_cross("BC");
  auto f2 = form_cross("F2");

  assert(bc.cross_type == "BC");
  assert(f2.cross_type == "F2");

  assert(bc.myfuncA() == 1);
  assert(bc.myfuncB() == 1);

  assert(f2.myfuncA() == 2);
  assert(f2.myfuncB() == 3);

  assert(anotherfunc(bc, "A") == 1);
  assert(anotherfunc(bc, "B") == 1);

  assert(anotherfunc(f2, "A") == 2);
  assert(anotherfunc(f2, "B") == 3);
}
