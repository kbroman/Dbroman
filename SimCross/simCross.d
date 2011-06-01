#!/usr/bin/rdmd
/**********************************************************************
 * simCross.d
 *
 * Karl Broman, 31 May 2011
 *
 * code to simulate genotype data for a simple cross
 **********************************************************************/

import std.stdio, std.conv;

void main()
{
  readMap("map10.txt");
}

void readMap(in string filename)
{
  // prepare to read file
  auto f = File(filename, "r");

  // read no. markers (and no. columns)
  int nMarkers, ncol;
  f.readf("%d %d\n", &nMarkers, &ncol);

  // skip a line
  string junk;
  f.readf("%s\n", &junk);
  
  writeln("no. markers=", nMarkers);

  auto markers = new string [](nMarkers);
  char[] tempchr;
  double temppos;
  double[string] pos;
  string[string] chr;

  foreach(i; 0 .. nMarkers) {
    f.readf("%s %s %f\n", &(markers[i]), &tempchr, &temppos);
    pos[markers[i]] = temppos;
    chr[markers[i]] = to!string(tempchr);
  }

  foreach(marker; markers) {
    if(chr[marker] == "12") 
      writefln("%-10s %2s %6.2f", marker, chr[marker], pos[marker]);
  }

}

// end of simCross.d
