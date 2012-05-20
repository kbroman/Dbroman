
import std.stdio, std.random;

// simple test of uniform random number generator
unittest {
  Random gen;

  double a;
  foreach(i; 0..10) {
    a = uniform(0.0, 1.0, gen);
    writef("%8.6f ", a);
  }
  writeln();
  writeln();

  gen.seed(unpredictableSeed);

  foreach(i; 0..10) {
    a = uniform(0.0, 1.0, gen);
    writef("%8.6f ", a);
  }
  writeln();
  writeln();

  gen.seed(87513765);

  foreach(i; 0..10) {
    a = uniform(0.0, 1.0, gen);
    writef("%8.6f ", a);
  }
  writeln();

  gen.seed(unpredictableSeed);
  writeln("dice:");
  foreach(i; 0..30)
    writef("%d", dice(gen, 1, 1));
  writeln;
  foreach(i; 0..30)
    writef("%d", dice(gen, 1, 1, 1));
  writeln();

  gen.seed(87513765);
  writeln("dice:");
  foreach(i; 0..30)
    writef("%d", dice(gen, 1, 1));
  writeln;
  foreach(i; 0..30)
    writef("%d", dice(gen, 25, 25, 25));
  writeln();
}
