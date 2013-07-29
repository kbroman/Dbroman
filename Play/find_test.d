import std.stdio, std.algorithm;

// playing with countUntil and find for an array of strings
unittest {
  auto xvect = ["blah", "house", "car", "airplane"];

  writeln("countUntil asda:", countUntil!("a == b")(xvect, "asda"));
  writeln("countUntil blah:", countUntil!("a == b")(xvect, "blah"));
  writeln("countUntil house:", countUntil!("a == b")(xvect, "house"));
  writeln("countUntil car:", countUntil!("a == b")(xvect, "car"));
  writeln("countUntil airplane:", countUntil!("a == b")(xvect, "airplane"));

  writeln;

  writeln("find asda:", find!("a == b")(xvect, "asda"));
  writeln("find blah:", find!("a == b")(xvect, "blah"));
  writeln("find house:", find!("a == b")(xvect, "house"));
  writeln("find car:", find!("a == b")(xvect, "car"));
  writeln("find airplane:", find!("a == b")(xvect, "airplane"));
}