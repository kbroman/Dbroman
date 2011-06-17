/**********************************************************************
 * calcprob.d
 * Karl Broman
 * first written 8 June 2011
 * last modified 17 June 2011
 * 
 * calculate genotype probabilities at each generation for RIL by sib-mating
 * 
 **********************************************************************/

import std.stdio, std.conv, std.algorithm, std.regex;

void main(string[] args) {

  assert(args.length > 3, "Give no. strains (2 or 4), chromosome type (A or X) and no. generations.");

  int n_strains = to!int(args[1]);
  char chr_type = args[2][0];
  int n_gen = to!int(args[3]);

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

  auto pairLookup = getPairLookup(pairs, chr_type, n_strains, true);

  string start;
  if(args.length > 4) {
    start = args[4];
  }
  else {
    if(chr_type == 'A') {
      if(n_strains == 2) {
        start = "AAxBB";
      }
      else {
        start = "ABxCD";
      }
    }
    else {
      if(n_strains == 2) {
        start = "AAxB";
      }
      else {
        start = "ABxC";
      }
    }
  }
  assert(start in pairLookup, "Incompatible starting point.");

  // check whether to swap AB
  if(n_strains==2 && chr_type == 'A') {
    auto pairLookup2 = getPairLookup(pairs, chr_type, n_strains, false);

    auto swapstart = exchangeLetters(start, ['A':'B', 'B':'A']);
    if(pairLookup2[swapstart] != pairLookup2[start]) {
      writeln("Can't swap A and B");
      pairLookup = pairLookup2;
    }
  }

  auto prototypes = getPrototypes(pairLookup);
  prototypes = prototypes.sort;
  assert(xiny(start, prototypes), "Starting point not a prototype.");

  auto transitionMatrix = getTM(pairLookup, prototypes, chr_type);

  double[string][int] probs;

  foreach(prototype; prototypes) {
    probs[0][prototype] = 0.0;
  }
  probs[0][start] = 1.0;

  foreach(gen; 1..(n_gen+1)) {
    probs[gen] = getNextGen(probs[gen-1], transitionMatrix);
  }

  double[string][int] indprob;
  string whichInd;
  string[string] indLookup;
  if(args.length > 5) {
    whichInd = args[5];
    assert((whichInd=="ind" && chr_type=='A') ||
           ((whichInd=="female" || whichInd=="male") && chr_type=='X'),
           "chr_type and whichInd incompatible.");
    
    indLookup = getIndLookup(prototypes, chr_type, n_strains, whichInd);

    foreach(gen; 0..(n_gen+1)) {
      indprob[gen] = getIndProbs(probs[gen], whichInd, indLookup);
    }
    writeProbs(indprob, indprob[0].keys.sort);
  }
  else {
    writeProbs(probs, prototypes);
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




string[] getOffspringA(string parentpair) 
{
  auto individuals = split(parentpair, regex("x"));

  auto offspring = new string[individuals[0].length * individuals[1].length];

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
  assert(getOffspringA("AAxAA") == ["AA", "AA", "AA", "AA"]);
  assert(getOffspringA("AAxAB") == ["AA", "AB", "AA", "AB"]);
  assert(getOffspringA("ABxCD") == ["AC", "AD", "BC", "BD"]);
  assert(getOffspringA("AAxA") == ["AA", "AA"]);
  assert(getOffspringA("AAxB") == ["AB", "AB"]);
  assert(getOffspringA("ABxC") == ["AC", "BC"]);
}  


string[] getOffspringXmale(string parentpair)
{
  auto individuals = split(parentpair, regex("x"));

  auto offspring = new string[individuals[0].length];

  int k = 0;
  foreach(i; individuals[0]) {
    offspring[k] = [i];
    k++;
  }

  return offspring;
}

unittest {
  assert(getOffspringXmale("AAxA") == ["A", "A"]);
  assert(getOffspringXmale("AAxB") == ["A", "A"]);
  assert(getOffspringXmale("ABxC") == ["A", "B"]);
  assert(getOffspringXmale("ABxA") == ["A", "B"]);
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


string[] getEquivPairs(string parentpair, char chr_type, int n_strains, bool swapAB)
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

  if(!(chr_type == 'X' && n_strains == 2) && swapAB) {
    // switch A<->B
    foreach(pair; pairHash.keys) {
      pairHash[to!string(exchangeLetters(pair, ['A':'B', 'B':'A']))] = 1;
    }
  }

  if(chr_type == 'A' && n_strains > 2) {
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
  assert(getEquivPairs("AAxAA", 'A', 4, true).sort == ["AAxAA", "BBxBB", "CCxCC", "DDxDD"]);
  assert(getEquivPairs("AAxAB", 'A', 4, true).sort == ["AAxAB", "AAxBA", "ABxAA", "ABxBB", 
						       "BAxAA", "BAxBB", "BBxAB", "BBxBA",
						       "CCxCD", "CCxDC", "CDxCC", "CDxDD",
						       "DCxCC", "DCxDD", "DDxCD", "DDxDC"]);
  assert(getEquivPairs("AAxAB", 'A', 2, true).sort == ["AAxAB", "AAxBA", "ABxAA", "ABxBB", 
						       "BAxAA", "BAxBB", "BBxAB", "BBxBA"]);
  assert(getEquivPairs("AAxAB", 'A', 4, true).sort == ["AAxAB", "AAxBA", "ABxAA", "ABxBB", 
						       "BAxAA", "BAxBB", "BBxAB", "BBxBA",
						       "CCxCD", "CCxDC", "CDxCC", "CDxDD",
						       "DCxCC", "DCxDD", "DDxCD", "DDxDC"]);

  assert(getEquivPairs("AAxA", 'X', 2, true).sort == ["AAxA"]);
  assert(getEquivPairs("ABxA", 'X', 2, true).sort == ["ABxA","BAxA"]);
  assert(getEquivPairs("ABxA", 'X', 4, true).sort == ["ABxA","ABxB", "BAxA","BAxB"]);
  assert(getEquivPairs("ABxC", 'X', 4, true).sort == ["ABxC","BAxC"]);
  assert(getEquivPairs("ACxA", 'X', 4, true).sort == ["ACxA","BCxB", "CAxA", "CBxB"]);
}

string[] getEquivInd(string ind, char chr_type, int n_strains)
{
  int[string] indHash;

  indHash[ind] = 1;

  if(!(chr_type == 'X' && n_strains == 2)) {
    // switch A<->B
    foreach(i; indHash.keys) {
      indHash[to!string(exchangeLetters(i, ['A':'B', 'B':'A']))] = 1;
    }
  }

  if(chr_type == 'A' && n_strains > 2) {
    // switch C<->D
    foreach(i; indHash.keys) {
      indHash[to!string(exchangeLetters(i, ['C':'D', 'D':'C']))] = 1;
    }

    // switch A<->C, B<->D
    foreach(i; indHash.keys) {
      indHash[to!string(exchangeLetters(i, ['A':'C', 'C':'A', 'B':'D', 'D':'B']))] = 1;
    }
  }
  
  return indHash.keys;
}

unittest {
  assert(getEquivInd("AA", 'A', 4).sort == ["AA", "BB", "CC", "DD"]);
  assert(getEquivInd("AB", 'A', 4).sort == ["AB", "BA", "CD", "DC"]); 

  assert(getEquivInd("AA", 'X', 2).sort == ["AA"]);
  assert(getEquivInd("AB", 'X', 2).sort == ["AB"]);
  assert(getEquivInd("AB", 'X', 4).sort == ["AB","BA"]);
  assert(getEquivInd("AC", 'X', 4).sort == ["AC","BC"]);
}


string[string] getPairLookup(string[] parentpairs, char chr_type, int n_strains, bool swapAB) 
{
  string[string] seenPairs;

  foreach(pair; parentpairs) {
    if(pair in seenPairs) {
      foreach(anEquivPair; getEquivPairs(pair, chr_type, n_strains, swapAB)) {
        seenPairs[anEquivPair] = seenPairs[pair];
      }
    }
    else {
      foreach(anEquivPair; getEquivPairs(pair, chr_type, n_strains, swapAB)) {
        seenPairs[anEquivPair] = pair;
      }
    }
  }

  return seenPairs;
}


string[string] getIndLookup(string[] parentpairs, char chr_type, int n_strains, string whichInd) 
{
  string[string] seenInd;

  string[] inds;
  foreach(pair; parentpairs) {
    inds = split(pair, regex("x"));
    
    if(whichInd=="female") {
      inds = [inds[0]];
    }
    else if(whichInd == "male") {
      inds = [inds[1]];
    }

    foreach(ind; inds) {
      if(ind in seenInd) {
        foreach(anEquivInd; getEquivInd(ind, chr_type, n_strains)) {
          seenInd[anEquivInd] = seenInd[ind];
        }
      }
      else {
        foreach(anEquivInd; getEquivInd(ind, chr_type, n_strains)) {
          seenInd[anEquivInd] = ind;
        }
      }
    }
  }
    
  return seenInd;
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


double[string][string] getTM(string[string] pairLookup, string[] prototypes, char chr_type)
{
  double[string][string] result;

  string[] pairs;
  foreach(prototype; prototypes) {
    if(chr_type == 'A') {
      auto offspring = getOffspringA(prototype);
      pairs = getPairs(offspring, offspring, "x");
    }
    else {
      auto female_offspring = getOffspringA(prototype);
      auto male_offspring = getOffspringXmale(prototype);
      pairs = getPairs(female_offspring, male_offspring, "x");
    }

    foreach(pair; pairs) {
      ++result[prototype][pairLookup[pair]];
    }
  }

  int rowSum;
  foreach(i; result[prototypes[0]].keys) {
    rowSum += result[prototypes[0]][i];
  }

  foreach(i; result.keys) {
    foreach(j; result[i].keys) {
      result[i][j] /= rowSum;
    }
  }

  foreach(i; prototypes) {
    foreach(j; prototypes) {
      if(j !in result[i]) {
        result[i][j] = 0.0;
      }
    }
  }

  return result;
}

double[string] getNextGen(double[string]thisgen, double[string][string] transitionMatrix)
{
  auto prototypes = thisgen.keys;

  double[string] nextgen;
  foreach(i; prototypes) {
    nextgen[i] = 0.0;
    foreach(j; prototypes) {
      nextgen[i] += thisgen[j] * transitionMatrix[j][i];
    }
  }
  return nextgen;
}

double[string] getIndProbs(double[string] pairprobs, string whichInd, string[string] indLookup)
{
  string[] pairs = pairprobs.keys;
  double[string] indprob; 
  string[] inds;

  assert(whichInd == "ind" || whichInd == "female" || whichInd == "male");

  double factor=1.0;

  foreach(pair; pairs) {
    inds = split(pair, regex("x"));

    
    if(whichInd=="female") {
      inds = [inds[0]];
    }
    else if(whichInd == "male") {
      inds = [inds[1]];
    }
    else {
      factor = 0.5;
    }

    foreach(ind; inds) {
      if(indLookup[ind] in indprob) {
        indprob[indLookup[ind]] += pairprobs[pair]*factor;
      }
      else {
        indprob[indLookup[ind]] = pairprobs[pair]*factor;
      }
    }
  }

  return indprob;
}

unittest {
  auto x = getIndProbs(["AAxAA":0.5, "AAxAB":0.5], "ind", ["AA":"AA", "AB":"AB"]);
  assert(x.keys.sort == ["AA","AB"]);
  assert(x["AA"] == 0.75);
  assert(x["AB"] == 0.25);

  x = getIndProbs(["AAxA":0.5, "AAxB":0.5], "female", ["AA":"AA"]);
  assert(x.keys == ["AA"]);
  assert(x["AA"] == 1.0);

  x = getIndProbs(["AAxA":0.5, "AAxB":0.5], "male", ["A":"A", "B":"B"]);
  assert(x.keys.sort == ["A","B"]);
  assert(x["A"] == 0.5);
  assert(x["B"] == 0.5);
}


void writeProbs(double[string][int] probs, string[] prototypes) 
{
  write("    ");
  foreach(i; prototypes) {
    writef("%9s ", i);
  }
  writeln();

  foreach(gen; 0..probs.length) {
    writef("%2d  ", gen);
    foreach(i; prototypes) {
      writef("%9.7f ", probs[gen][i]);
    }
    writeln();
  }
}

bool xiny(string x, string[] y)
{
  foreach(yv; y) {
    if(x == yv) return(true);
  }
  return(false);
}

/* end of calcprob.d */
