import std.stdio;

void main() {}

unittest {
  // array test
  double y[] = [12.0, 123.1, 15.3];
  writeln("array test");
  writeln("y.length: ", y.length);
  
  auto x = y;
  writeln("x is y? ", x is y);
  writeln("x.length: ", x.length);
  write("x: ");
  foreach (xv; x) 
    writef("%5.1f ", xv);
  writeln();

  writeln("Swapping last two values");
  auto temp = y[2];
  y[2] = y[1];
  y[1] = temp;
  writeln("x is y? ", x is y);


  writeln("y.length: ", y.length);
  write("y: ");
  foreach (xv; y) 
    writef("%5.1f ", xv);
  writeln();
  
  writeln("x.length: ", x.length);
  write("x: ");
  foreach (xv; x) 
    writef("%5.1f ", xv);
  writeln();


  writeln("&(x[0]) is &(y[0])? ", &(x[0]) is &(y[0]));
  
  writeln("add to y");
  y ~= -13.9;
  writeln("x is y? ", x is y);
  writeln("&(x[0]) is &(y[0])? ", &(x[0]) is &(y[0]));

  writeln("y.length: ", y.length);
  write("y: ");
  foreach (xv; y) 
    writef("%5.1f ", xv);
  writeln();
  
  writeln("x.length: ", x.length);
  write("x: ");
  foreach (xv; x) 
    writef("%5.1f ", xv);
  writeln();
  
  writeln("Swapping last two values");
  temp = y[3];
  y[3] = y[2];
  y[2] = temp;

  writeln("y.length: ", y.length);
  write("y: ");
  foreach (xv; y) 
    writef("%5.1f ", xv);
  writeln();
  
  writeln("x.length: ", x.length);
  write("x: ");
  foreach (xv; x) 
    writef("%5.1f ", xv);
  writeln();
}
