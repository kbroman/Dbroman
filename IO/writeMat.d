#!/usr/bin/rdmd
/****************************************
 * writeMat.d
 * Karl Broman, 31 May 2011
 * 
 * read two matrices from text files
 * write them to binary files
 * read them in again
 * multiply them
 ****************************************/

import std.stdio;

void main() {
  string Atxtfile="A.txt", Abinfile="A.bin";
  string Btxtfile="B.txt", Bbinfile="B.bin";

  auto Amat = readMatrix(Atxtfile);
  writeBinMatrix(Amat, Abinfile);

  auto Bmat = readMatrix(Btxtfile);
  writeBinMatrix(Bmat, Bbinfile);

  auto Amat2 = readBinMatrix(Abinfile);
  auto Bmat2 = readBinMatrix(Bbinfile);
  
  writeln("Amat2:");
  printMatrix(Amat2);

  writeln("\nBmat2:");
  printMatrix(Bmat2);

  auto Product = matmult(Amat2, Bmat2);
  writeln("\nA * B:");
  printMatrix(Product);
}



// read matrix from text file
double [][]readMatrix(string filename)
{
  write("Reading from ", filename, ":");
  auto f = File(filename, "r");
  int ncol, nrow;

  f.readf("%d %d\n", &nrow, &ncol);
  writeln("    nrow=", nrow, " ncol=", ncol);
  auto result = new double [][](nrow, ncol);
  
  foreach(i; 0..nrow) {
    foreach(ref el; result[i]) {
      f.readf("%s", &el);
    }
  }
  return result;
}


// write matrix to binary file
void writeBinMatrix(double [][] matrix, string filename)
{
  auto f = File(filename, "w");
  f.rawWrite([matrix.length, matrix[0].length]);
  foreach(i; 0 .. matrix.length) {
    f.rawWrite(matrix[i]);
  }
}
  
// read matrix from binary file
double[][] readBinMatrix(string filename)
{
  write("Reading from ", filename, ":");
  auto f = File(filename, "r");

  auto dim = new uint[2];
  f.rawRead(dim);
  auto nrow = dim[0];
  auto ncol = dim[1];

  writeln("    nrow=", nrow, "  ncol=", ncol);

  auto result = new double[][](nrow, ncol);
  foreach(ref el; result) {
    f.rawRead(el);
  }
  return result;
}
  
// print matrix
void printMatrix(T)(T[][] matrix)
{
  foreach(i; 0 .. matrix.length) {
    foreach(j; 0 .. matrix[0].length) {
      writef("%s ", matrix[i][j]);
    }
    write("\n");
  }
}

// multiply two matrices
double [][]matmult(double [][]firstMatrix, double [][]secondMatrix)
{
  auto result = new double[][](firstMatrix.length, secondMatrix[0].length);

  // should check that secondMatrix.length == firstMatrix[0].length

  foreach(i; 0 .. firstMatrix.length) {
    foreach(j; 0 .. secondMatrix[0].length) {
      result[i][j] = 0.0;
      foreach(k; 0 .. firstMatrix[0].length) {
        result[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
      }
    }
  }
  return result;
}

// end of writeMat.d
