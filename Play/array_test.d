import std.stdio;

// writevec: print out a vector

unittest {
  // array test
  size_t y[] = [3, 7, 8];
  writeln("array test");
  writeln("y.length: ", y.length);
  
  writevec(y);

  foreach(i, ref yv; y)
    yv = i;

  writevec(y);
}


void writevec(T)(T[] input)
{
  foreach(x; input) {
    writef("%d ", x);
  }
  writeln();
}  