/**********************************************************************
 * countIndStates.d
 * Karl Broman
 * first written 8 July 2011
 * last modified 8 July 2011
 * 
 * Count states two loci, sib mating
 * 
 **********************************************************************/

import std.stdio, std.conv, std.algorithm, std.regex;

void main(string[] args) {

  assert(args.length > 2, "Give no. strains (2 or 4 or 8), chromosome type (A or X)");

  int n_strains = to!int(args[1]);
  char chr_type = args[2][0];

  writeln("no. strains = ", n_strains);
  writeln("chr_type = ", chr_type);

  auto strains = ["A", "B", "C", "D", "E", "F", "G", "H"];
  assert(n_strains <= strains.length, "Can't use more than " ~ to!string(strains.length) ~ " strains.");
  assert(n_strains % 2 == 0, "No. strains must be multiple of 2.");
  assert(chr_type == 'A' || chr_type == 'X', "Chromosome type must be 'A' or 'X'.");

  if(chr_type=='X') {
    switch(n_strains) {
    case 2: strains = strains[0..n_strains]; break;
    case 4: strains = strains[0..3]; break;
    case 8: strains = [strains[0], strains[1], strains[2], strains[4], strains[5]]; break;
    }
  }
  else {
    strains = strains[0..n_strains];
  }

  writeln("Construct individuals");
  auto haplotypes = getPairs(strains, strains, "");
  auto individuals = getPairs(haplotypes, haplotypes, "|");

  writeln("No. haplotypes = ", haplotypes.length);
  writeln("No. individuals = ", individuals.length);

  auto indLookup = getIndLookup(individuals, chr_type, n_strains); 

  if(args.length <= 3 || args[3] == "") {
    auto prototypes = getPrototypes(indLookup);
    writeln("No. prototypes = ", prototypes.length);
  }
  else if(to!int(args[3])==1) {
    auto prototypes = getPrototypeCounts(indLookup);
    writeln("No. prototypes = ", prototypes.keys.length);
    auto sorted_prototypes = prototypes.keys.sort;
    int total=0;
    foreach(prototype; sorted_prototypes) {
      writefln("%5s  %3d", prototype, prototypes[prototype]);
      total += prototypes[prototype];
    }
    writeln("Total = ", total);
  }
  else { 
    auto prototypes = getPrototypeCounts(indLookup);
    writeln("No. prototypes = ", prototypes.keys.length);
    auto sorted_prototypes = prototypes.keys.sort;

    foreach(i; individuals) {
      writefln("%5s  %5s", i, indLookup[i]);
    }
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
  assert(exchangeLetters("AA|AA", ['A':'B', 'B':'A']) == "BB|BB");
  assert(exchangeLetters("AA|AB", ['A':'B', 'B':'A']) == "BB|BA");
  assert(exchangeLetters("AA|AC", ['A':'C', 'C':'A']) == "CC|CA");
  assert(exchangeLetters("AA|AC", ['B':'C', 'C':'B']) == "AA|AB");

  assert(exchangeLetters("AAxC", ['A':'B', 'B':'A']) == "BBxC");
  assert(exchangeLetters("ACxB", ['A':'B', 'B':'A']) == "BCxA");
}



string switchLoci(string individual)
{
  auto haplotypes = split(individual, regex("[\\|]"));

  foreach(i, ref hap; haplotypes) {
    hap = [hap[1]] ~ [hap[0]];
  }

  return haplotypes[0] ~ "|" ~ haplotypes[1];
}

unittest {
  writeln("[unittest switchLoci]");
  assert(switchLoci("AA|AA") == "AA|AA");
  assert(switchLoci("AA|AB") == "AA|BA");
  assert(switchLoci("AA|BB") == "AA|BB");
  assert(switchLoci("CD|DD") == "DC|DD");
  assert(switchLoci("AB|CD") == "BA|DC");
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

  if(chr_type=='A') {
    // switch A<->B
    foreach(i; indHash.keys) {
      indHash[to!string(exchangeLetters(i, ['A':'B', 'B':'A']))] = 1;
    }
    if(n_strains > 2) {
      // switch C<->D
      foreach(i; indHash.keys) {
	indHash[to!string(exchangeLetters(i, ['C':'D', 'D':'C']))] = 1;
      }

      // switch A<->C, B<->D
      foreach(i; indHash.keys) {
        indHash[to!string(exchangeLetters(i, ['A':'C', 'C':'A', 'B':'D', 'D':'B']))] = 1;
      }
    }
    if(n_strains > 4) {
      // switch E<->F
      foreach(i; indHash.keys) {
	indHash[to!string(exchangeLetters(i, ['E':'F', 'F':'E']))] = 1;
      }

      // switch G<->H
      foreach(i; indHash.keys) {
	indHash[to!string(exchangeLetters(i, ['G':'H', 'H':'G']))] = 1;
      }

      // switch E<->G, F<->H
      foreach(i; indHash.keys) {
        indHash[to!string(exchangeLetters(i, ['E':'G', 'G':'E', 'F':'H', 'H':'F']))] = 1;
      }
      // switch A<->E, B<->F, C<->G, D<->H
      foreach(i; indHash.keys) {
        indHash[to!string(exchangeLetters(i, ['A':'E', 'E':'A', 'B':'F', 'F':'B',
					      'C':'G', 'G':'C', 'D':'H', 'H':'D']))] = 1;
      }
    }
  }
  else {
    if(n_strains > 2) {
      // switch A<->B
      foreach(i; indHash.keys) {
	indHash[to!string(exchangeLetters(i, ['A':'B', 'B':'A']))] = 1;
      }
    }
    if(n_strains > 4) {
      // switch E<->F
      foreach(i; indHash.keys) {
	indHash[to!string(exchangeLetters(i, ['E':'F', 'F':'E']))] = 1;
      }
    }
  }
  
  return indHash.keys;
}

/*
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
*/

string[string] getIndLookup(string[] inds, char chr_type, int n_strains) 
{
  string[string] seenInd;

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
    
  return seenInd;
}


string[] getPrototypes(string[string] indLookup)
{
  int[string] prototypes;

  foreach(ind; indLookup.keys) {
    prototypes[indLookup[ind]] = 1;
  }

  return prototypes.keys;
}

int[string] getPrototypeCounts(string[string] indLookup)
{
  int[string] prototypes;

  foreach(ind; indLookup.keys) {
    prototypes[indLookup[ind]] = prototypes[indLookup[ind]] + 1;
  }

  return prototypes;
}




/* end of countIndStates.d */
