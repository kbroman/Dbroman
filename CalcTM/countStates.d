/**********************************************************************
 * countStates.d
 * Karl Broman
 * first written 20 June 2011
 * last modified 20 June 2011
 * 
 * Count states one locus, RIL by sibmating
 * 
 **********************************************************************/

import std.stdio, std.conv, std.algorithm, std.regex;

void main(string[] args) {

  assert(args.length > 2, "Give no. strains (2 or 4) and chromosome type (A or X).");

  int n_strains = to!int(args[1]);
  char chr_type = args[2][0];

  auto strains = ["A", "B", "C", "D"];
  assert(n_strains <= strains.length, "Can't use more than " ~ to!string(strains.length) ~ " strains.");
  assert(n_strains % 2 == 0, "No. strains must be multiple of 2.");
  assert(chr_type == 'A' || chr_type == 'X', "Chromosome type must be 'A' or 'X'.");

  strains = strains[0..n_strains];

  string[] pairs;
  if(chr_type == 'A') {
    auto individuals = getPairs(strains, strains, "");
    pairs = getPairs(individuals, individuals, "x");
  }
  else {
    if(n_strains==4) {
      strains = strains[0 .. 3];
      n_strains = 3;
    }
    auto females = getPairs(strains, strains, "");
    auto males = strains.dup;
    pairs = getPairs(females, males, "x");
  }

  bool swapLetters;
  if(args.length > 3) {
    swapLetters = to!bool(args[3]);
  }
  else {
    swapLetters = true;
  }

  auto pairLookup = getPairLookup(pairs, chr_type, n_strains, true, swapLetters);

  // for 2 strains, autosome, allow starting point different than AAxBB
  //    then check that A and B can be swapped
  string start;
  if(n_strains==2 && chr_type == 'A' && swapLetters) {
    auto pairLookup2 = getPairLookup(pairs, chr_type, n_strains, false, swapLetters);
    if(args.length > 4) {
      start = args[4];
    }
    else {
        start = "AAxBB";
    }
    assert(start in pairLookup, "Incompatible starting point.");

    auto swapstart = exchangeLetters(start, ['A':'B', 'B':'A']);
    if(pairLookup2[swapstart] != pairLookup2[start]) {
      writeln("Can't swap A and B");
      pairLookup = pairLookup2;
    }
  }

  auto prototypes = getPrototypes(pairLookup);
  prototypes = prototypes.sort;

  /* count states */
  int[string] n_states;
  foreach(pair; pairs) {
    ++n_states[pairLookup[pair]];
  }

  foreach(prototype; prototypes) {
    writefln("%5s %2d", prototype, n_states[prototype]);
  }


}


string[] getPairs(string[] symbols_left, string[] symbols_right, string sep)
{
  auto n_pair = symbols_left.length*symbols_right.length;
  auto pairs = new string[n_pair];

  int k = 0;
  foreach(i; symbols_left) {
    foreach(j; symbols_right) {
      pairs[k] = i ~ sep ~ j;
      k++;
    }
  }

  return pairs;
}

unittest {
  assert(getPairs(["A"], ["A"], "") == ["AA"]);
  assert(getPairs(["A", "B"], ["A", "B"], "") == ["AA", "AB", "BA", "BB"]);

  assert(getPairs(["A"], ["A"], "x") == ["AxA"]);
  assert(getPairs(["A", "B"], ["A", "B"], "x") == ["AxA", "AxB", "BxA", "BxB"]);

  assert(getPairs(["AA"], ["AA"], "x") == ["AAxAA"]);
  assert(getPairs(["AA", "AB", "BB"], ["AA", "AB", "BB"], "x") == 
         ["AAxAA", "AAxAB", "AAxBB", "ABxAA", "ABxAB", "ABxBB", "BBxAA", "BBxAB", "BBxBB"]);
  assert(getPairs(["AA", "AB", "BB"], ["A", "B"], "x") == 
         ["AAxA", "AAxB", "ABxA", "ABxB", "BBxA", "BBxB"]);
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

  assert(switchOneInd("ABxA", 0) == "BAxA");
  assert(switchOneInd("ABxC", 0) == "BAxC");
  assert(switchOneInd("CDxB", 0) == "DCxB");

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

  assert(exchangeLetters("AAxC", ['A':'B', 'B':'A']) == "BBxC");
  assert(exchangeLetters("ACxB", ['A':'B', 'B':'A']) == "BCxA");
}


string[] getEquivPairs(string parentpair, char chr_type, int n_strains, bool swapAB, bool swapLetters)
{
  int[string] pairHash;

  pairHash[parentpair] = 1;

  if(chr_type == 'A') {
    // switch pair
    foreach(pair; pairHash.keys) {
      pairHash[switchPair(pair)] = 1;
    }

    // switch letters in second
    foreach(pair; pairHash.keys) {
      pairHash[switchOneInd(pair, 1)] = 1;
    }
  }

  // switch letters in first
  foreach(pair; pairHash.keys) {
    pairHash[switchOneInd(pair, 0)] = 1;
  }

  if(!(chr_type == 'X' && n_strains == 2) && swapAB && swapLetters) {
    // switch A<->B
    foreach(pair; pairHash.keys) {
      pairHash[to!string(exchangeLetters(pair, ['A':'B', 'B':'A']))] = 1;
    }
  }

  if(chr_type == 'A' && n_strains > 2 && swapLetters) {
    // switch C<->D
    foreach(pair; pairHash.keys) {
      pairHash[to!string(exchangeLetters(pair, ['C':'D', 'D':'C']))] = 1;
    }

    // switch A<->C, B<->D
    foreach(pair; pairHash.keys) {
      pairHash[to!string(exchangeLetters(pair, ['A':'C', 'C':'A', 'B':'D', 'D':'B']))] = 1;
    }
  }
  
  return pairHash.keys;
}

unittest {
  assert(getEquivPairs("AAxAA", 'A', 4, true, true).sort == ["AAxAA", "BBxBB", "CCxCC", "DDxDD"]);
  assert(getEquivPairs("AAxAB", 'A', 4, true, true).sort == ["AAxAB", "AAxBA", "ABxAA", "ABxBB", 
						       "BAxAA", "BAxBB", "BBxAB", "BBxBA",
						       "CCxCD", "CCxDC", "CDxCC", "CDxDD",
						       "DCxCC", "DCxDD", "DDxCD", "DDxDC"]);
  assert(getEquivPairs("AAxAB", 'A', 2, true, true).sort == ["AAxAB", "AAxBA", "ABxAA", "ABxBB", 
						       "BAxAA", "BAxBB", "BBxAB", "BBxBA"]);
  assert(getEquivPairs("AAxAB", 'A', 4, true, true).sort == ["AAxAB", "AAxBA", "ABxAA", "ABxBB", 
						       "BAxAA", "BAxBB", "BBxAB", "BBxBA",
						       "CCxCD", "CCxDC", "CDxCC", "CDxDD",
						       "DCxCC", "DCxDD", "DDxCD", "DDxDC"]);

  assert(getEquivPairs("AAxA", 'X', 2, true, true).sort == ["AAxA"]);
  assert(getEquivPairs("ABxA", 'X', 2, true, true).sort == ["ABxA","BAxA"]);
  assert(getEquivPairs("ABxA", 'X', 4, true, true).sort == ["ABxA","ABxB", "BAxA","BAxB"]);
  assert(getEquivPairs("ABxC", 'X', 4, true, true).sort == ["ABxC","BAxC"]);
  assert(getEquivPairs("ACxA", 'X', 4, true, true).sort == ["ACxA","BCxB", "CAxA", "CBxB"]);
}


string[string] getPairLookup(string[] parentpairs, char chr_type, int n_strains, bool swapAB, bool swapLetters) 
{
  string[] equivPairs;
  string[string] seenPairs;

  foreach(pair; parentpairs) {
    if(pair in seenPairs) {
      foreach(anEquivPair; getEquivPairs(pair, chr_type, n_strains, swapAB, swapLetters)) {
        seenPairs[anEquivPair] = seenPairs[pair];
      }
    }
    else {
      foreach(anEquivPair; getEquivPairs(pair, chr_type, n_strains, swapAB, swapLetters)) {
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


/* end of countStates.d */
