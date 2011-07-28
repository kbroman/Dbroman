// generic functions

module generics; 

import std.stdio;

immutable GENOTYPE_NA = -1;

enum F2  { NA = GENOTYPE_NA, A, H, B, HorB, HorA }; 
enum BC  { NA = GENOTYPE_NA, A, H };

F2[] possibleGeno(F2 one)
{
  auto possible_genotypes = [F2.A, F2.H, F2.B];
  
  return possible_genotypes;
}

BC[] possibleGeno(BC one)
{
  auto possible_genotypes = [BC.A, BC.H];

  return possible_genotypes;
}


unittest {
  assert(possibleGeno(BC.A) == [BC.A, BC.H]);
  assert(possibleGeno(F2.A) == [F2.A, F2.H, F2.B]);
}
  

void main()
{
  writeln("BC: ", possibleGeno(BC.A));

  writeln("F2: ", possibleGeno(F2.A));
}