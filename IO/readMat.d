#!/usr/bin/rdmd
/****************************************
 * readMat.d
 * Karl Broman, 31 May 2011
 * 
 * read two matrices from text files and 
 * calculate their product
 ****************************************/

import std.stdio;

void main() {
  auto Amat = readMatrix("A.txt");
  auto Bmat = readMatrix("B.txt");

  auto Cmat = matmult(Amat, Bmat);

  writeln("\nresult:");
  printMatrix(Cmat);
}


// read a matrix from a text file
double [][]readMatrix(string filename)
{
  write("Reading from ", filename, ":");
  auto f = File(filename, "r");
  int ncol, nrow;

  f.readf("%d %d\n", &nrow, &ncol);
  writeln("    nrow=", nrow, " ncol=", ncol);
  auto result = new double [][](nrow, ncol);
  
  foreach(i; 0 .. nrow) {
    foreach(j; 0 .. ncol) {
      f.readf("%s", &(result[i][j]));
    }
  }
  return result;
}


// multiply two matrices
double [][]matmult(double [][]Amat, double [][]Bmat)
{
  auto Cmat = new double[][](Amat.length, Bmat[0].length);

  // should ensure that Bmat.length == Amat[0].length

  foreach(i; 0 .. Amat.length) {
    foreach(j; 0 .. Bmat[0].length) {
      Cmat[i][j] = 0.0;
      foreach(k; 0 .. Amat[0].length) {
        Cmat[i][j] += Amat[i][k] * Bmat[k][j];
      }
    }
  }
  return Cmat;
}


// print a matrix
void printMatrix(T)(T[][] matrix)
{
  foreach(i; 0 .. matrix.length) {
    foreach(j; 0 .. matrix[0].length) {
      writef("%s ", matrix[i][j]);
    }
    write("\n");
  }
}
 
// end of readMat.d
