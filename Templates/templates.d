// try template/mixin type things

module templates.d;

import std.stdio;

immutable GENOTYPE_NA = -1;
enum F2  { NA = GENOTYPE_NA, A, H, B, HorB, HorA }; 
enum BC  { NA = GENOTYPE_NA, A, H };

mixin template grabfifthCode(T)
{
  T grabfifth(T[] vector)
  {
    if(vector.length < 5) 
      throw new Exception("vector is too short");

    return vector[4];
  }
}


mixin grabfifthCode!F2;
mixin grabfifthCode!BC;


void main()
{
  auto bcvect = [BC.A, BC.H, BC.A, BC.H, BC.A];
  auto f2vect = [F2.A, F2.H, F2.A, F2.H, F2.B];

  writeln("F2: ", grabfifth(f2vect));
  
  writeln("BC: ", grabfifth(bcvect));
  
}