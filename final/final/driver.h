//
//  Filename:     driver.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the driver header for Assignment 6
//

#ifndef hw3_driver_h
#define hw3_driver_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "matrix.h"
#include "matrix_triangular_upper.h"
#include "matrix_triangular_lower.h"
#include "matrix_diagonal.h"
#include "matrix_tridiagonal.h"
#include "matrix_symmetrical.h"
#include "matrix_banded.h"
#include "matrix_poisson.h"
#include "cholesky.h"
#include "gauss_elim.h"
#include "gauss_seidel.h"
#include "pde_test.h"
#include "pde_final.h"


typedef enum
{
	kMatrixTypeDense = 0,
	kMatrixTypeTriangularUpper = 1,
	kMatrixTypeTriangularLower = 2,
	kMatrixTypeDiagonal = 3,
	kMatrixTypeTridiagonal = 4,
	kMatrixTypeSymmetrical = 5
} kMatrixType;


double actualSolution(double aX, double aY)
{
  return ( 1.0 / sinh(M_PI) ) * ( sin(aX) * sinh(M_PI - aY) + sin(aY) * sinh(M_PI - aX));
}


//Description: Find solution, b for Ax = b for the matrix and vector
//             contained in aFile
//Pre:         see openFile for format conditions
//Post:        results are diaplayed in console
void solveFile(const char * aFile, kMatrixType aType);

//Description: Open a file for matrix input
//Pre:         file should be in proper format
//
//                  n
//                  A(0,0) .. A(0,n)
//
//                  .   .     .
//
//                  .      .  .
//
//                  A(n,0) .. A(n,n)
//
//                  x[0] .. x[n]
//
//             where:
//                n       = the number of rows, columns and equations
//                A(i, j) = the element of aMatrix at the ith row and jth column
//                x[i]    = the ith element of aVector
//
//Post:        aMatrix, aVector are set up with the contents of a file called aFile
bool openFile(matrix_base<double>* aMatrix, vector<double>& aVector, const char * aFile);

void run(const size_t aN, pde_base<double>& aPDE);

template <class T_method>
bool solveMatrix(vector<double>& aX, const matrix_base<double>& aMatrix, const vector<double>& aB, T_method aMethod);

void displayMatrixTypes();

//Description: Run tests for matrix classes
//Pre:         none
//Post:        test results are displayed in the console
void testMatrices();

//Description: Run test
//Pre:         none
//Post:        aName is displayed and is aExpression is true PASS, otherwise FAIL
void test(bool aExpression, std::string aName);

#endif
