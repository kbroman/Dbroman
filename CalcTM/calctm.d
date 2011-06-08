/**********************************************************************
 * calctm.d
 * Karl Broman, 7 June 2011
 * 
 * calculate transition matrix for one locus, RIL by sibmating
 * 
 **********************************************************************/

import std.stdio, std.conv, std.algorithm, std.regex;

void main(string[] args) {

  assert(args.length > 1, "Give no. strains");

  int n_str = to!int(args[1]);

  auto strains = ["A", "B", "C", "D"];
  assert(n_str <= strains.length, "Can't use more than " ~ to!string(strains.length) ~ " strains");
  assert(n_str % 2 == 0, "No. strains must be multiple of 2");

  strains = strains[0..n_str];

  auto individuals = getPairs(strains, "");

  auto pairs = getPairs(individuals, "x");

  auto pairLookup = getPairLookup(pairs);

  auto prototypes = getPrototypes(pairLookup);
  prototypes = prototypes.sort;

  auto transitionMatrix = getTM(pairLookup, prototypes);

  if(args.length > 2 && args[2] == "maxima") {
    writeTMmaxima(transitionMatrix, prototypes, "P" ~ to!string(n_str));
  }
  else {
    writeTM(transitionMatrix, prototypes);
  }
}


string[] getPairs(string[] symbols, string sep)
{
  auto n_pair = symbols.length^^2;
  auto pairs = new string[n_pair];

  int k = 0;
  foreach(i; symbols) {
    foreach(j; symbols) {
      pairs[k] = i ~ sep ~ j;
      k++;
    }
  }

  return pairs;
}

unittest {
  assert(getPairs(["A"], "") == ["AA"]);
  assert(getPairs(["A", "B"], "") == ["AA", "AB", "BA", "BB"]);

  assert(getPairs(["A"], "x") == ["AxA"]);
  assert(getPairs(["A", "B"], "x") == ["AxA", "AxB", "BxA", "BxB"]);

  assert(getPairs(["AA"], "x") == ["AAxAA"]);
  assert(getPairs(["AA", "AB", "BB"], "x") == 
	 ["AAxAA", "AAxAB", "AAxBB", "ABxAA", "ABxAB", "ABxBB", "BBxAA", "BBxAB", "BBxBB"]);
}


string[] getOffspring(string parentpair)
{
  auto individuals = split(parentpair, regex("x"));

  auto offspring = new string[4];

  int k = 0;
  foreach(i; individuals[0]) {
    foreach(j; individuals[1]) {
      offspring[k] = [i] ~ [j];
      k++;
    }
  }
  return offspring;
}

unittest {
  assert(getOffspring("AAxAA") == ["AA", "AA", "AA", "AA"]);
  assert(getOffspring("AAxAB") == ["AA", "AB", "AA", "AB"]);
  assert(getOffspring("ABxCD") == ["AC", "AD", "BC", "BD"]);
}  

string switchPair(string parentpair) 
{
  auto individuals = split(parentpair, regex("x"));

  auto newpair = individuals[1] ~ "x" ~ individuals[0];

  return newpair;
}

unittest {
  assert(switchPair("AAxAA") == "AAxAA");
  assert(switchPair("AAxAB") == "ABxAA");
  assert(switchPair("ABxCD") == "CDxAB");
}

string switchOneInd(string parentpair, int which2switch) 
{
  auto individuals = split(parentpair, regex("x"));

  individuals[which2switch] = [individuals[which2switch][1]] ~ [individuals[which2switch][0]];
  auto newpair = individuals[0] ~ "x" ~ individuals[1];

  return newpair;
}
  
unittest {
  assert(switchOneInd("AAxAA", 0) == "AAxAA");
  assert(switchOneInd("AAxAB", 0) == "AAxAB");
  assert(switchOneInd("ABxCD", 0) == "BAxCD");

  assert(switchOneInd("AAxAA", 1) == "AAxAA");
  assert(switchOneInd("AAxAB", 1) == "AAxBA");
  assert(switchOneInd("ABxCD", 1) == "ABxDC");
}

char[] exchangeLetters(string parentpair, char[char] char2switch) 
{
  auto newpair = new char[parentpair.length];

  foreach(i; 0..parentpair.length) {
    if(parentpair[i] in char2switch) {
      newpair[i] = char2switch[parentpair[i]];
    }
    else {
      newpair[i] = parentpair[i];
    }
  }
  return newpair;
}

unittest {
  assert(exchangeLetters("AAxAA", ['A':'B', 'B':'A']) == "BBxBB");
  assert(exchangeLetters("AAxAB", ['A':'B', 'B':'A']) == "BBxBA");
  assert(exchangeLetters("AAxAC", ['A':'C', 'C':'A']) == "CCxCA");
  assert(exchangeLetters("AAxAC", ['B':'C', 'C':'B']) == "AAxAB");
}


string[] getEquivPairs(string parentpair)
{
  int[string] pairHash;

  pairHash[parentpair] = 1;

  // switch pair
   foreach(pair; pairHash.keys) {
     pairHash[switchPair(pair)] = 1;
   }

  // switch letters in first
  foreach(pair; pairHash.keys) {
    pairHash[switchOneInd(pair, 0)] = 1;
  }

  // switch letters in second
  foreach(pair; pairHash.keys) {
    pairHash[switchOneInd(pair, 1)] = 1;
  }

  // switch A<->B
  foreach(pair; pairHash.keys) {
    pairHash[to!string(exchangeLetters(pair, ['A':'B', 'B':'A']))] = 1;
  }

  // switch C<->D
  foreach(pair; pairHash.keys) {
    pairHash[to!string(exchangeLetters(pair, ['C':'D', 'D':'C']))] = 1;
  }

  // switch A<->C, B<->D
  foreach(pair; pairHash.keys) {
    pairHash[to!string(exchangeLetters(pair, ['A':'C', 'C':'A', 'B':'D', 'D':'B']))] = 1;
  }
  
  return pairHash.keys;
}

unittest {
  assert(getEquivPairs("AAxAA").sort == ["AAxAA", "BBxBB", "CCxCC", "DDxDD"]);
  assert(getEquivPairs("AAxAB").sort == ["AAxAB", "AAxBA", "ABxAA", "ABxBB", 
					 "BAxAA", "BAxBB", "BBxAB", "BBxBA",
					 "CCxCD", "CCxDC", "CDxCC", "CDxDD",
					 "DCxCC", "DCxDD", "DDxCD", "DDxDC"]);
}


string[string] getPairLookup(string[] parentpairs) 
{
  string[] equivPairs;
  string[string] seenPairs;

  foreach(pair; parentpairs) {
    if(pair in seenPairs) {
      foreach(anEquivPair; getEquivPairs(pair)) {
	seenPairs[anEquivPair] = seenPairs[pair];
      }
    }
    else {
      foreach(anEquivPair; getEquivPairs(pair)) {
	seenPairs[anEquivPair] = pair;
      }
    }
  }

  return seenPairs;
}




string[] getPrototypes(string[string] pairLookup)
{
  int[string] prototypes;

  foreach(pair; pairLookup.keys) {
    prototypes[pairLookup[pair]] = 1;
  }

  return prototypes.keys;
}

unittest {
  assert(getPrototypes(["A":"B", "B":"B", "C":"D", "D":"D", "E":"A"]).sort == ["A", "B", "D"]);
}


int[string][string] getTM(string[string] pairLookup, string[] prototypes)
{
  int[string][string] result;

  foreach(prototype; prototypes) {
    auto offspring = getOffspring(prototype);
    auto pairs = getPairs(offspring, "x");
    auto n_pairs = pairs.length;

    foreach(pair; pairs) {
      ++result[prototype][pairLookup[pair]];
    }
  }

  return result;
}

void writeTM(int[string][string] tm, string[] prototypes)
{
  write("       ");
  foreach(i; prototypes) {
    write(i, "  ");
  }
  writeln();

  int rowSum;
  foreach(i; tm[prototypes[0]].keys) {
    rowSum += tm[prototypes[0]][i];
  }

  foreach(i; prototypes) {
    write(i, "  ");
    foreach(j; prototypes) {
      if(j in tm[i]) {
	writef("%6.4f ", cast(double)tm[i][j]/rowSum);
      }
      else {
	writef("%6.4f ", 0.0);
      }
    }
    writeln();
  }
}

void writeTMmaxima(int[string][string] tm, string[] prototypes, string label="P")
{
  writeln(label, " :  matrix(");

  int rowSum;
  foreach(i; tm[prototypes[0]].keys) {
    rowSum += tm[prototypes[0]][i];
  }

  write("/*               ");
  foreach(i; prototypes) {
    writef("  %s   ", i);
  }
  writeln("*/");


  foreach(i; prototypes) {
    write("/* ", i, " */    [ ");
    foreach(j; prototypes) {
      if(j in tm[i]) {
	writef("%4d/%-2d ", tm[i][j], rowSum);
      }
      else {
	writef(" %4d   ", 0);
      }
      if(j !is prototypes[$-1]) {
	write(", ");
      }
    }
    if(i !is prototypes[$-1]) {
      writeln("],");
    }
    else {
      writeln("] )$");
    }
  }
}


/* end of calctm.d */
