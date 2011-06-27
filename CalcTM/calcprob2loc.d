/**********************************************************************
 * calcprob2loc.d
 * Karl Broman
 * first written 27 June 2011
 * last modified 27 June 2011
 * 
 * calculate genotype probabilities at each generation, two loci RIL by sibmating
 * 
 **********************************************************************/

import std.stdio, std.conv, std.algorithm, std.regex;

void main(string[] args) {

  assert(args.length > 4, "Give no. strains (2 or 4), chromosome type (A or X), rec'n fraction and no. generations.");

  int n_strains = to!int(args[1]);
  char chr_type = args[2][0];
  double recfrac = to!double(args[3]);
  int n_gen = to!int(args[4]);

  auto strains = ["A", "B", "C", "D"];
  assert(n_strains <= strains.length, "Can't use more than " ~ to!string(strains.length) ~ " strains.");
  assert(n_strains % 2 == 0, "No. strains must be multiple of 2.");
  assert(chr_type == 'A' || chr_type == 'X', "Chromosome type must be 'A' or 'X'.");

  strains = strains[0..n_strains];

  string[] pairs;
  writeln("Construct individuals and pairs");
  string[] haplotypes;
  if(chr_type == 'A') {
    haplotypes = getPairs(strains, strains, "");
    auto individuals = getPairs(haplotypes, haplotypes, "|");
    pairs = getPairs(individuals, individuals, "x");
  }
  else {
    if(n_strains==4) {
      strains = strains[0 .. 3];
      n_strains = 3;
    }
    haplotypes = getPairs(strains, strains, "");
    auto females = getPairs(haplotypes, haplotypes, "|");
    auto males = haplotypes.dup;
    pairs = getPairs(females, males, "x");
  }

  writeln("    No. pairs: ", pairs.length);

  writeln("Get pair lookup");
  auto pairLookup = getPairLookup(pairs, chr_type, n_strains, true);

  // for 2 strains, autosome, allow starting point different than AA|AAxBB|BB
  //    then check that A and B can be swapped
  string start;
  if(args.length > 5) {
    start = args[5];
  }
  else {
    if(chr_type == 'A') {
      if(n_strains == 2) {
        start = "AA|AAxBB|BB";
      }
      else {
        start = "AA|BBxCC|DD";
      }
    }
    else {
      if(n_strains == 2) {
        start = "AA|AAxBB";
      }
      else {
        start = "AA|BBxCC";
      }
    }
  }
  assert(start in pairLookup, "Incompatible starting point.");

  if(n_strains==2 && chr_type == 'A') {
    auto pairLookup2 = getPairLookup(pairs, chr_type, n_strains, false);

    auto swapstart = exchangeLetters(start, ['A':'B', 'B':'A']);
    if(pairLookup2[swapstart] != pairLookup2[start]) {
      writeln("    Note: Can't swap A and B");
      pairLookup = pairLookup2;
    }
  }

  writeln("Get prototypes");
  auto prototypes = getPrototypes(pairLookup);
  prototypes = prototypes.sort;
  writeln("    No. prototypes: ", prototypes.length);

  if(n_strains==2 && chr_type=='A') {
    assert(xiny(start, prototypes), "Starting point not a prototype.");
  }
    
  writeln("Get transition matrix");
  auto transitionMatrix = getTM(pairLookup, prototypes, chr_type, recfrac);

  writeln("Get probs");
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
  if(args.length > 6) {
    whichInd = args[6];
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
  writeln("[unittest getPairs]");
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


double[string] getMeioticProducts(string individual, double recfrac)
{
  double[string] result;

  auto haplotypes = split(individual, regex("\\|"));

  if(haplotypes.length==1) { // male on X chromosome 
    result[haplotypes[0]]=1.0;
    return(result);
  }

  foreach(hap; haplotypes) {
    result[hap] += 0.5*(1-recfrac);
  }
  if(haplotypes.length==2) {
    result[[haplotypes[0][0]] ~ [haplotypes[1][1]]] += 0.5*recfrac;
    result[[haplotypes[1][0]] ~ [haplotypes[0][1]]] += 0.5*recfrac;
  }

  return result;
}

unittest {
  writeln("[unittest getMeioticProducts]");
  assert(getMeioticProducts("AA|AB", 0.05) == ["AA":0.5, "AB":0.5]);
  assert(getMeioticProducts("AB|AB", 0.05) == ["AB":1.0]);
  assert(getMeioticProducts("AA|BB", 0.05) == ["AA":0.475, "BB":0.475, "AB":0.025, "BA":0.025]);

  assert(getMeioticProducts("AA|AB", 0.05) == ["AA":0.5, "AB":0.5]);
  assert(getMeioticProducts("AB|CD", 0.1) == ["AB":0.45, "CD":0.45, "AD":0.05, "CB":0.05]);
  assert(getMeioticProducts("AA|BB", 0.05) == ["AA":0.475, "BB":0.475, "AB":0.025, "BA":0.025]);

  assert(getMeioticProducts("AA", 0.05) == ["AA":1.0]);
  assert(getMeioticProducts("AD", 0.05) == ["AD":1.0]);
}


double[string] getOffspringA(string parentpair, double recfrac)
{
  auto individuals = split(parentpair, regex("x"));

  auto offspring1 = getMeioticProducts(individuals[0], recfrac);
  auto offspring2 = getMeioticProducts(individuals[1], recfrac);
  double[string] result;

  foreach(kid1; offspring1.keys) {
    foreach(kid2; offspring2.keys) {
      result[kid1 ~ "|" ~ kid2] += offspring1[kid1] * offspring2[kid2];
    }
  }
  
  return result;
}

unittest {
  writeln("[unittest getOffspringA]");
  assert(getOffspringA("AA|AAxBB|BB", 0.05) == ["AA|BB":1.0]);
  assert(getOffspringA("AA|BBxAA|AA", 0.2) == ["AA|AA":0.4, "BB|AA":0.4,
					       "AB|AA":0.1, "BA|AA":0.1]);
  assert(getOffspringA("AA|ABxCC|DD", 0.1) == ["AA|CC":0.225, "AA|DD":0.225,
					       "AB|CC":0.225, "AB|DD":0.225,
					       "AA|CD":0.025, "AA|DC":0.025,
					       "AB|CD":0.025, "AB|DC":0.025]);

  double recfrac = 0.1;
  double[] v = [((1-recfrac)/2)^^2, recfrac*(1-recfrac)/4, recfrac^^2/4];
  assert(getOffspringA("AA|BBxCC|DD", 0.1) == ["AA|CC":v[0], "AA|DD":v[0],
					       "BB|CC":v[0], "BB|DD":v[0],
					       "AB|CC":v[1], "AB|DD":v[1],
					       "BA|CC":v[1], "BA|DD":v[1],
					       "AA|CD":v[1], "AA|DC":v[1],
					       "BB|CD":v[1], "BB|DC":v[1],
					       "AB|CD":v[2], "AB|DC":v[2],
					       "BA|CD":v[2], "BA|DC":v[2]]);
					       

  recfrac = 0.01;
  v = [((1-recfrac)/2)^^2, recfrac*(1-recfrac)/4, recfrac^^2/4];
  assert(getOffspringA("AA|BBxCC|DD", 0.01) == ["AA|CC":v[0], "AA|DD":v[0],
					       "BB|CC":v[0], "BB|DD":v[0],
					       "AB|CC":v[1], "AB|DD":v[1],
					       "BA|CC":v[1], "BA|DD":v[1],
					       "AA|CD":v[1], "AA|DC":v[1],
					       "BB|CD":v[1], "BB|DC":v[1],
					       "AB|CD":v[2], "AB|DC":v[2],
					       "BA|CD":v[2], "BA|DC":v[2]]);
					       
}  


double[string] getOffspringXmale(string parentpair, double recfrac)
{
  auto individuals = split(parentpair, regex("x"));

  return getMeioticProducts(individuals[0], recfrac);
}

unittest {
  writeln("[unittest getOffspringXmale]");
  assert(getOffspringXmale("AA|AAxAA", 0.05) == ["AA":1.0]);
  assert(getOffspringXmale("AA|AAxBB", 0.05) == ["AA":1.0]);

  assert(getOffspringXmale("AA|ABxBB", 0.05) == ["AA":0.5, "AB":0.5]);
  assert(getOffspringXmale("AA|ABxCC", 0.05) == ["AA":0.5, "AB":0.5]);

  assert(getOffspringXmale("AA|BBxCC", 0.2) == ["AA":0.4, "BB":0.4,
						"AB":0.1, "BA":0.1]);
}  


string switchPair(string parentpair) 
{
  auto individuals = split(parentpair, regex("x"));

  auto newpair = individuals[1] ~ "x" ~ individuals[0];

  return newpair;
}

unittest {
  writeln("[unittest switchPair]");
  assert(switchPair("AA|AAxAA|AA") == "AA|AAxAA|AA");
  assert(switchPair("AA|ABxAA|BB") == "AA|BBxAA|AB");
  assert(switchPair("AA|BBxCC|DD") == "CC|DDxAA|BB");
}

string switchOneInd(string parentpair, int which2switch) 
{
  auto individuals = split(parentpair, regex("x"));
  auto ind2swap = split(individuals[which2switch], regex("\\|"));

  individuals[which2switch] = ind2swap[1] ~ "|" ~ ind2swap[0];

  return individuals[0] ~ "x" ~ individuals[1];
}
  
unittest {
  writeln("[unittest switchOneInd]");
  assert(switchOneInd("AA|AAxAA|AA", 0) == "AA|AAxAA|AA");
  assert(switchOneInd("AA|AAxAA|AB", 0) == "AA|AAxAA|AB");
  assert(switchOneInd("AA|ABxCC|DD", 0) == "AB|AAxCC|DD");

  assert(switchOneInd("AA|BBxAA", 0) == "BB|AAxAA");
  assert(switchOneInd("AA|BBxCC", 0) == "BB|AAxCC");
  assert(switchOneInd("CD|DDxBB", 0) == "DD|CDxBB");

  assert(switchOneInd("AA|BBxAA|BB", 1) == "AA|BBxBB|AA");
  assert(switchOneInd("AA|ABxCD|AB", 1) == "AA|ABxAB|CD");
  assert(switchOneInd("AB|CDxAB|CD", 1) == "AB|CDxCD|AB");
}



string switchHap(string individual)
{
  auto haplotypes = split(individual, regex("\\|"));

  if(haplotypes.length==2) {
    return haplotypes[1] ~ "|" ~ haplotypes[0];
  }
  else { 
    return individual;
  }
}
  
unittest {
  writeln("[unittest switchHap]");
  assert(switchHap("AA|AA") == "AA|AA");
  assert(switchHap("AA|AB") == "AB|AA");
  assert(switchHap("CC|DD") == "DD|CC");

  assert(switchHap("AA|BB") == "BB|AA");
  assert(switchHap("CD|DD") == "DD|CD");

  assert(switchHap("BA") == "BA");
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
  writeln("[unittest exchangeLetters]");
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

  // switch loci
  foreach(pair; pairHash.keys) {
    pairHash[switchLoci(pair)] = 1;
  }

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
  writeln("[unittest getEquivPairs]");
  assert(getEquivPairs("AA|AAxAA|AA", 'A', 4, true).sort == ["AA|AAxAA|AA", "BB|BBxBB|BB", 
							     "CC|CCxCC|CC", "DD|DDxDD|DD"]);
  assert(getEquivPairs("AA|AAxAB|AB", 'A', 2, true).sort == ["AA|AAxAB|AB", "AA|AAxBA|BA",
							     "AB|ABxAA|AA", "AB|ABxBB|BB",
							     "BA|BAxAA|AA", "BA|BAxBB|BB",
							     "BB|BBxAB|AB", "BB|BBxBA|BA"]);
}


string switchLoci(string parentpair)
{
  string[] haplotypes = split(parentpair, regex("[\\|x]"));

  foreach(i, ref hap; haplotypes) {
    hap = [hap[1]] ~ [hap[0]];
  }

  if(haplotypes.length==4) {
    return haplotypes[0] ~ "|" ~ haplotypes[1] ~ "x" ~
      haplotypes[2] ~ "|" ~ haplotypes[3];
  }
  else if(haplotypes.length==3) {
    return haplotypes[0] ~ "|" ~ haplotypes[1] ~ "x" ~
      haplotypes[2];
  }
  else if(haplotypes.length==2) {
    return haplotypes[0] ~ "|" ~ haplotypes[1];
  }
  else {
    return haplotypes[0];
  }
}

unittest {
  writeln("[unittest switchLoci]");
  assert(switchLoci("AA|AAxAA|AA") == "AA|AAxAA|AA");
  assert(switchLoci("AA|AAxAA|AB") == "AA|AAxAA|BA");
  assert(switchLoci("AA|ABxCC|DD") == "AA|BAxCC|DD");

  assert(switchLoci("AA|BBxAA") == "AA|BBxAA");
  assert(switchLoci("AA|BBxAC") == "AA|BBxCA");
  assert(switchLoci("CD|DDxBC") == "DC|DDxCB");

  assert(switchLoci("AA|BBxAA|BB") == "AA|BBxAA|BB");
  assert(switchLoci("AA|ABxCD|AB") == "AA|BAxDC|BA");
  assert(switchLoci("AB|CDxAB|CD") == "BA|DCxBA|DC");

  assert(switchLoci("AA|BB") == "AA|BB");
  assert(switchLoci("AA|AB") == "AA|BA");
  assert(switchLoci("AB|CD") == "BA|DC");

  assert(switchLoci("AA") == "AA");
  assert(switchLoci("CD") == "DC");
  assert(switchLoci("AB") == "BA");
}



string[] getEquivInd(string ind, char chr_type, int n_strains)
{
  int[string] indHash;

  indHash[ind] = 1;

  foreach(i; indHash.keys) {
    indHash[to!string(switchLoci(i))] = 1;
  }

  foreach(i; indHash.keys) {
    indHash[to!string(switchHap(i))] = 1;
  }

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
  assert(getEquivInd("AA|BB", 'A', 4).sort == ["AA|BB", "BB|AA", "CC|DD", "DD|CC"]);
					       
  assert(getEquivInd("AA|AB", 'A', 4).sort == ["AA|AB", "BB|BA", "CC|CD", "DD|DC",
					       "AA|BA", "BB|AB", "CC|DC", "DD|CD",
					       "AB|AA", "BA|BB", "CD|CC", "DC|DD",
					       "BA|AA", "AB|BB", "DC|CC", "CD|DD"].sort);

  assert(getEquivInd("AA|AB", 'X', 4).sort == ["AA|AB", "BB|BA", "AA|BA", "BB|AB",
					       "AB|AA", "BA|BB", "BA|AA", "AB|BB"].sort);

  assert(getEquivInd("AA", 'X', 4).sort == ["AA", "BB"]);
  assert(getEquivInd("AA", 'X', 2).sort == ["AA"]);
  assert(getEquivInd("AB", 'X', 2).sort == ["AB","BA"]);
}

string[string] getPairLookup(string[] parentpairs, char chr_type, int n_strains, bool swapAB) 
{
  string[] equivPairs;
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
  writeln("[unittest getPrototypes]");
  assert(getPrototypes(["A":"B", "B":"B", "C":"D", "D":"D", "E":"A"]).sort == ["A", "B", "D"]);
}


double[string] getPairProb(double[string] femaleProb, double[string] maleProb)
{
  double[string] result;

  foreach(female; femaleProb.keys) {
    foreach(male; maleProb.keys) {
      result[female ~ "x" ~ male] += femaleProb[female] * maleProb[male];
    }
  }

  return result;
}

unittest {
  writeln("[unittest getPairProb]");
  assert(getPairProb(["AA|AA":1.0], ["CC":1.0]) == ["AA|AAxCC":1.0]);
  assert(getPairProb(["AA|AA":1.0], ["CC":0.5, "DD":0.5]) == ["AA|AAxCC":0.5, "AA|AAxDD":0.5]);
  assert(getPairProb(["AA|AB":0.2, "AA|BA":0.8], ["CC":0.5, "DD":0.5]) == ["AA|ABxCC":0.1, "AA|ABxDD":0.1,
									   "AA|BAxCC":0.4, "AA|BAxDD":0.4]);
}
  

double[string][string] getTM(string[string] pairLookup, string[] prototypes, char chr_type, double rec_frac)
{
  double[string][string] result;

  double[string] pairProb;
  foreach(prototype; prototypes) {
    if(chr_type == 'A') {
      auto offspring = getOffspringA(prototype, rec_frac);
      pairProb = getPairProb(offspring, offspring);
    }
    else {
      auto female_offspring = getOffspringA(prototype, rec_frac);
      auto male_offspring = getOffspringXmale(prototype, rec_frac);
      pairProb = getPairProb(female_offspring, male_offspring);
    }

    foreach(pair; pairProb.keys) {
      result[prototype][pairLookup[pair]] += pairProb[pair];
    }

    foreach(pro2; prototypes) {
      if(pro2 !in result[prototype])
	result[prototype][pro2] = 0.0;
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
    writef("%11s ", i);
  }
  writeln();

  foreach(gen; 0..probs.length) {
    writef("%2d  ", gen);
    foreach(i; prototypes) {
      writef("%11.9f ", probs[gen][i]);
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

/* end of calcprob2loc.d */
