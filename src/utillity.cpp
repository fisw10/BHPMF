//
//  utillity.cpp
//
//
//  Created by Farideh Fazayeli on 7/18/13.
//
//

#include "utillity.h"

void mvrnorm(double *des, double *mu, double *cholCov, int dim, bool upper){
  
  int i;
  int inc = 1;
  double one = 1.0;
  //make some std norm draws
  for(i = 0; i < dim; i++)
    des[i] = rnorm(0.0, 1.0);
  
  //mult this vector by the lower triangle of the cholCov
  if(upper)
    F77_NAME(dtrmv)("U", "T", "N", &dim, cholCov, &dim, des, &inc);
  else
    F77_NAME(dtrmv)("L", "N", "N", &dim, cholCov, &dim, des, &inc);
  
  //add the mean to the result
  F77_NAME(daxpy)(&dim, &one, mu, &inc, des, &inc);
  
  
}

/*
Print an array, assumes array dimension is known
*/
void PrintArr(double *myArr, int nRow){
  for (int ii = 0; ii < nRow; ii++) {
    Rprintf("%5.3f ", myArr[ii]);
  }
  Rprintf("\n");
}

void PrintArr(int *myArr, int nRow){
  
  for (int ii = 0; ii < nRow; ii++) {
    Rprintf("%5.3f ", myArr[ii]);
    //cout << myArr[ii] << " ";
  }
  Rprintf("\n");
  //cout << endl;
}

/*
Print a matrix stored in an array, assumes matrix dimension is known
*/
void printMat(float *myArr, int nRow, int nCol){
  
  for (int ii = 0; ii < nRow; ii++) {
    for (int jj = 0; jj < nCol; jj++) {
      Rprintf("%5.3f ", myArr[ii*nCol+jj]);
      //cout << myArr[ii*nCol+jj] << " ";
    }
    Rprintf("\n");
    //cout << endl;
  }
}

void printMat(double *myArr, int nRow, int nCol){
  
  for (int ii = 0; ii < nRow; ii++) {
    for (int jj = 0; jj < nCol; jj++) {
      Rprintf("%5.3f ", myArr[ii*nCol+jj]);
    }
    Rprintf("\n");
  }
}

/*
Print a matrix from a given matrix, assumes matrix dimension is known
*/
void printMat(float **myMat, int nRow, int nCol){
  
  for (int ii = 0; ii < nRow; ii++) {
    for (int jj = 0; jj < nCol; jj++) {
      Rprintf("%5.3f ", myMat[ii][jj]);
      //cout << myMat[ii][jj] << " ";
    }
    Rprintf("\n");
    //cout << endl;
  }
}


void printMat(int **myMat, int nRow, int nCol){
  
  for (int ii = 0; ii < nRow; ii++) {
    for (int jj = 0; jj < nCol; jj++) {
      Rprintf("%5.3f ", myMat[ii][jj]);
      //cout << myMat[ii][jj] << " ";
    }
    Rprintf("\n");
    //cout << endl;
  }
}

SEXP getListElement (SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  
  for (i = 0; i < length(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }
  return elmt;
}

/*
void dlmreadVec(char *filePath, matVecFloat &myMat, int &nRow, int nCol, int verbose)

Read a text file into a matrix, a very simple code with assumption number of columns of matrix is known.
*/
void dlmreadVec(string filePath, matVecFloat &myMat, int &nRow, int nCol, int verbose){
  
  vector<float>  myArr;
  FILE * pFile;
  pFile = fopen (filePath.c_str(), "r");
  
  if (!pFile)
  {
    // cout << "Cannot open file "<< filePath << endl;
    //cerr << "Failed to open " << pFile << endl;
    //            exit(1);  //abort program
  }
  
  // cout << "ncol: " << nCol << endl;
  
  nRow = 0;
  float num;
  
  while(fscanf(pFile, "%f", &num) != EOF){
    myArr.push_back(num);
  }
  
  fclose(pFile);
  nRow = myArr.size()/nCol;
  
  myMat.resize(nRow);
  for(int ii = 0; ii < nRow; ii++){
    myMat[ii].resize(nCol);
    for (int jj = 0; jj < nCol; jj++)
      myMat[ii][jj] = myArr[ii*nCol+jj];
  }
  if (verbose) {
    Rprintf("There are exactly %i rows and %i columns in the matrix \n", nRow, nCol);
  }
}


/*
void dlmreadVec(char *filePath, matVecInt &myMat, int &nRow, int nCol, int verbose)

Read a text file into a matrix, a very simple code with assumption number of columns of matrix is known.
*/
void dlmreadVec(string filePath, matVecInt &myMat, int &nRow, int nCol, int verbose){
  
  vector<int> myArr;
  FILE * pFile;
  pFile = fopen (filePath.c_str(), "r");
  
  
  if (!pFile)
  {
    // cout << "Cannot open file "<< filePath << endl;
    //cerr << "Failed to open " << pFile << endl;
    //            exit(1);  //abort program
  }
  
  // cout <<"nCol: " <<  nCol << endl;
  
  nRow = 0;
  int num;
  
  while(fscanf(pFile, "%d", &num) != EOF){
    myArr.push_back(num);
  }
  
  fclose(pFile);
  
  nRow = myArr.size()/nCol;
  
  myMat.resize(nRow);
  for(int ii = 0; ii < nRow; ii++){
    myMat[ii].resize(nCol);
    for (int jj = 0; jj < nCol; jj++)
      myMat[ii][jj] = myArr[ii*nCol+jj];
  }
  if (verbose) {
    Rprintf("There are exactly %i rows and %i columns in the matrix \n", nRow, nCol);
  }
}



/*
void dlmreadVec(char *filePath, matVecInt &myMat, int &nRow, int nCol)

Read a text file into a matrix, a very simple code with assumption number of columns of matrix is known.
*/
void dlmreadVec(string filePath, vector<int> &myArr){
  
  FILE * pFile;
  pFile = fopen (filePath.c_str(), "r");
  
  
  if (!pFile)
  {
    // cout << "Cannot open file "<< filePath << endl;
    //cerr << "Failed to open " << pFile << endl;
    //            exit(1);  //abort program
  }
  
  int num;
  
  while(fscanf(pFile, "%d", &num) != EOF){
    myArr.push_back(num);
  }
  
  fclose(pFile);
  
}

/*
void dlmreadVec(string filePath, matVecInt &myMat, int nRow, vector<int> nCols, int verbose) {
  
  vector<int> myArr;
  FILE * pFile;
  pFile = fopen (filePath.c_str(), "r");
  
  
  if (!pFile)
  {
    // cout << "Cannot open file "<< filePath << endl;
    //cerr << "Failed to open " << pFile << endl;
    //            exit(1);  //abort program
  }
  int num;
  
  while(fscanf(pFile, "%d", &num) != EOF) {
    myArr.push_back(num);
  }
  
  fclose(pFile);
  
  myMat.resize(nRow);
  for(int ii = 0; ii < nRow; ii++) {
    myMat[ii].resize(nCols[ii]);
    for (int jj = 0; jj < nCols[ii]; jj++)
      myMat[ii][jj] = myArr[ii*nCols[ii]+jj];
  }
  
  if (verbose) {
    Rprintf("There are exactly %i rows and %i columns in the matrix \n", nRow, nCols);
  }
}
*/

void dlmwriteVec(const char *inFilePath, const char *outFilePath, int nRow,
                 int nCol, int gap, int burn, int nSams) {
  matVecFloat myMat;
  FILE * ipFile;
  ipFile = fopen (inFilePath, "r");
  
  if (!ipFile)
  {
    // cout << "Cannot open file "<< inFilePath << endl;
    //cerr << "Failed to open " << pFile << endl;
    //            exit(1);  //abort program
  }
  
  float num;
  
  myMat.resize(nRow);
  for(int ii = 0; ii < nRow; ii++) {
    myMat[ii].resize(nCol);
  }
  
  int count = 1;
  for(int sam = 0; sam < nSams; sam++){
    for(int ii = 0; ii < nRow; ii++){
      for(int jj = 0; jj < nCol; jj++){
        fscanf(ipFile, "%f", &num);
        if(sam == gap){
          myMat[ii][jj] = num;
        }
        else if ((sam > gap) & ((sam-gap) % burn == 0)){
          myMat[ii][jj] += num;
          if (!ii & !jj) count++;
        }
      }
    }
  }
  
  fclose(ipFile);
  
  FILE * opFile;
  opFile = fopen (outFilePath, "w");
  
  if (!opFile)
  {
    // cout << "Cannot open file "<< outFilePath << endl;
    //            exit(1);  //abort program
  }
  
  for(int ii = 0; ii < nRow; ii++){
    for (int jj = 0; jj < nCol; jj++){
      fprintf(opFile,"%5.4f \t", myMat[ii][jj]/count);
    }
    fprintf(opFile,"\n");
  }
  
  fclose(opFile);
  
}
