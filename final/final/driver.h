//
//  Filename:     driver.h
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the driver header for Assignment 7
//

#ifndef final_driver_h
#define final_driver_h

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
#include "runtime.h"

typedef enum
{
	kMatrixTypeDense = 0,
	kMatrixTypeTriangularUpper = 1,
	kMatrixTypeTriangularLower = 2,
	kMatrixTypeDiagonal = 3,
	kMatrixTypeTridiagonal = 4,
	kMatrixTypeSymmetrical = 5
} kMatrixType;


//Description: Calculates and returns the actual solution for the given problem
//Pre:
//Post:
double actualSolution(double aX, double aY)
{
  return ( 1.0 / sinh(M_PI) ) * ( sin(aX) * sinh(M_PI - aY) + sin(aY) * sinh(M_PI - aX));
}

//Description: Calculates and returns the actual solution for the 
//Pre:
//Post:
void run(const size_t aN, pde_base<double>& aPDE);

//Description: Calculates and returns the actual solution for the 
//Pre:
//Post:
void runSolvers(vector<double>& x, const matrix_base<double>& m, const vector<double>& b);

double checkError(const vector<double>& aX, const vector<point2d<double>>& aXMapping);

template <class T_method>
bool solveMatrix(vector<double>& aX, const matrix_base<double>& aMatrix, const vector<double>& aB, T_method aMethod);

//Description: Run tests for matrix classes
//Pre:         none
//Post:        test results are displayed in the console
void testMatrices();

//Description: Run test
//Pre:         none
//Post:        aName is displayed and is aExpression is true PASS, otherwise FAIL
void test(bool aExpression, std::string aName);

#endif
