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
#include <typeinfo>
#include <pthread.h>
#include <unistd.h>

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
#include "pde_notes.h"
#include "runtime.h"


typedef enum
{
	kMenuChoiceQuit = 0,
	kMenuChoiceChangeMeshDensity = 1,
  kMenuChoiceChangeBounds = 2,
	kMenuChoiceCompare = 3,
	kMenuChoiceSolve = 4,
	kMenuChoiceMatlab = 5,
	kMenuChoiceTest = 6
} kMenuChoice;


//Description:  print a message
//Pre:          none
//Post:         outputs the message
void printMessage(const std::string& aMessage);


//Description:  Calculates and returns the actual solution for the given problem
//Pre:          None
//Post:         returns the numerical solution at a certain (x,y) point
double actualSolution(double aX, double aY)
{
	return ( 1.0 / sinh(M_PI) ) * ( sin(aX) * sinh(M_PI - aY) + sin(aY) * sinh(M_PI - aX));
}


//Description: Starts menu
//Pre:         none
//Post:        outputs menu
void runMenu();

//Description: Prompts user to solve pde
//Pre:         none
//Post:        solves pde based on method chosen
unsigned int promptForSolve(pde_base<double>& aPDE);

bool loading = false;
void *startLoading(void *aMessage);

//Description: Creates the b and xMapping vectors for some N
//Pre:         none
//Post:        aB and aXMapping are populated with the approprate values for N
void createSystem(const pde_base<double>& aPDE, vector<double>& aB, vector<point2d<double>>& aXMapping);


//Description: Solves a aPDE
//Pre:         none
//Post:        aPDE is populated with the solution
template <class T_method>
void solvePDE(pde_base<double>& aPDE, T_method aMethod);


//Description: Runs each method on aPDE of mesh size aN
//Pre:         none
//Post:        returns each method's running time and error
std::string runSolvers(pde_base<double>& aPDE, const bool aShouldIterate = true, const bool aPrettyPrint = true);


//Description: Checks how close an approximation is the the actual solution
//Pre:         aX's size should ge greater than zero
//Post:        returns the error in percent
double checkError(const vector<double>& aX, const vector<point2d<double>>& aXMapping);


//Description: Calculates and returns the actual solution for the
//Pre: 				 matrix aA must be a (n x n) matrix
//             vector aB must be of size n
//						 T_method must be one of the defined solving methods
//Post:				 vector aX will contain the size n solution to Ax = b
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
