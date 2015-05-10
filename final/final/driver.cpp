//
//  Filename:			driver.cpp
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the driver implementation for Assignment 7
//                
//

#include "driver.h"

int main()
{
  pde_final<double> pde(0,M_PI);

//	MATLAB OUTPUT
//  solvePDE(25, pde);
//  std::cout << pde.matlabOutput() << std::endl;
  
  runSolvers(20, pde);
  
  testMatrices();
  
  return 0;
}


void createSystem(const size_t aN, const pde_base<double>& aPDE, vector<double>& aB, vector<point2d<double>>& aXMapping)
{
  const size_t SIZE = (aN-1)*(aN-1);
  double lowerBound = aPDE.lowerBound();
  double upperBound = aPDE.upperBound();
  double inc = fabs((upperBound - lowerBound)) / aN;
  size_t count = 0;
  double bValue = 0;
  
  aB.clear();
  aXMapping.clear();
  
  aB.reserve(SIZE, true);
  aXMapping.reserve(SIZE, true);
  count = 0;
  
  for (double i = lowerBound + inc; !equivalent(i, upperBound); i += inc)
  {
    for (double j = lowerBound + inc; !equivalent(j, upperBound); j += inc)
    {
      aXMapping[count].set(j, i);
      bValue = 0;
      
      for (int k = 0; k < 4; k++)
      {
        double x = j, y = i;
        
        switch (k)
        {
          case 0:
            x -= inc;
            break;
            
          case 1:
            y -= inc;
            break;
            
          case 2:
            x += inc;
            break;
            
          case 3:
            y += inc;
            break;
            
          default:
            break;
        }
        
        if (equivalent(x, lowerBound) || equivalent(y, lowerBound) || equivalent(x, upperBound) || equivalent(y, upperBound))
        {
          bValue += aPDE(x, y);
        }
      }
      aB[count] = bValue*0.25;
      count++;
    }
  }
}


void solvePDE(const size_t aN, pde_base<double>& aPDE)
{
  matrix_poisson<double> m(aN);
  vector<point2d<double>> xMapping;
  vector<double> b;
  vector<double> x;
  point2d<double> p;
  double iBound = 0;
  double jBound = 0;
  double lowerBound = aPDE.lowerBound();
  double upperBound = aPDE.upperBound();
  double inc = fabs((upperBound - lowerBound)) / aN;
  double tolerance = inc / 2.0;
  
  aPDE.clearPoints();
  
  createSystem(aN, aPDE, b, xMapping);
  
  solveMatrix(x, m, b, cholesky<double>());
  
  // changing the previous for loops with while loops prevents
  // round-off error
  
  iBound = lowerBound;
  
  while (std::abs(iBound-(upperBound+inc)) >= tolerance)
  {
    jBound = lowerBound;
    
    while (std::abs(jBound-(upperBound+inc)) >= tolerance)
    {
      if (equivalent(jBound, lowerBound) || equivalent(iBound, lowerBound) || equivalent(jBound, upperBound) || equivalent(iBound, upperBound))
      {
        p.set(jBound, iBound);
        aPDE.addKnownPoint(p);
      }
      
      jBound += inc;
    }
    
    iBound += inc;
  }
  
  for (size_t i = 0; i < x.size(); i++)
  {
    aPDE.addPoint(xMapping[i], x[i]);
  }
}


void runSolvers(const size_t aN, const pde_base<double>& aPDE, const bool aShouldIterate)
{
  matrix_poisson<double>* m;
  vector<point2d<double>> xMapping;
  vector<double> b;
  vector<double> x;
  vector<double> times;
  runtime timer;
  size_t start = (aShouldIterate ? 5 : aN);
  
  std::cout << std::setw(3) << "N" << std::setw(12) << "Cholesky" << std::setw(25) << "Gaussian Elimination"<< std::setw(15) << "Gauss-Seidel" << std::setw(12) << "Error" << std::endl << std::endl;
  
  for (size_t n = start; n <= aN; n++)
  { 
    m = new matrix_poisson<double>(n);
    createSystem(n, aPDE, b, xMapping);
    times.clear();
    
    timer.begin();
    solveMatrix(x, *m, b, cholesky<double>());
    timer.end();
    times.push_back(timer.elapsed());
    
    timer.begin();
    solveMatrix(x, *m, b, gauss_elim<double>());
    timer.end();
    times.push_back(timer.elapsed());
    
    timer.begin();
    solveMatrix(x, *m, b, gauss_seidel<double>());
    timer.end();
    times.push_back(timer.elapsed());
    
    delete m;
    
    std::cout << std::setw(3) << n
    << std::setw(12) << times[0]
    << std::setw(25) << times[1]
    << std::setw(15) << times[2]
    << std::setw(12) << checkError(x, xMapping)
    << std::endl;
  }
  std::cout << std::endl;
}


double checkError(const vector<double>& aX, const vector<point2d<double>>& aXMapping)
{
  double res = 0;
  
  for (size_t i = 0; i < aX.size(); i++)
  {
    res += fabs((aX[i] - actualSolution(aXMapping[i].x(), aXMapping[i].y())));
  }
  
  return (res / aX.size())*100.0;
}


template <class T_method>
bool solveMatrix(vector<double>& aX, const matrix_base<double>& aMatrix, const vector<double>& aB, T_method aMethod)
{
  return aMethod(aX, aMatrix, aB);
}


void testMatrices()
{ 
  // seed random
  srand(static_cast<unsigned int>(time(NULL)));
  
  // size of matrix to test (> 3)
  const size_t START_SIZE = 10;
  
  // vector of matrices
  vector<matrix_base<double>*> mats;
  
  // dummy variables for testing
  matrix<double> dummyDense;
  matrix<double> zero;
  matrix_triangular_lower<double> dummyLower;
  matrix_triangular_upper<double> dummyUpper;
  matrix_diagonal<double> dummyDiagonal;
  matrix_banded<double> dummyBanded;
  bool dummyBool;
  size_t dummySize;
  
  // matrices
  matrix<double> mDense;
  matrix_triangular_lower<double> mLower;
  matrix_triangular_upper<double> mUpper;
  matrix_banded<double> mBanded(10, 4);
  matrix_diagonal<double> mDiagonal;
  matrix_tridiagonal<double> mTridiagonal;
  matrix_symmetrical<double> mSymmetrical;
  matrix_poisson<double> mPoisson;
  matrix_poisson<double>* mPoissonPtr;
  
  // insert matrices into vector
  mats.push_back(&mDense);
  mats.push_back(&mLower);
  mats.push_back(&mUpper);
  mats.push_back(&mBanded);
  mats.push_back(&mDiagonal);
  mats.push_back(&mTridiagonal);
  mats.push_back(&mSymmetrical);
  mats.push_back(&mPoisson);
  
  //// TEST: tests
  
  ////  The following tests are tested for every matrix in the array
  ////  for both even and odd dimensions
  
  ////  this also tests polymorphism
  
  for (size_t size = START_SIZE; size > START_SIZE - 2; size--)
  { 
    for (size_t i = 0; i < mats.size(); i++)
    {
      matrix_base<double>* matrix = mats[i];
      
      matrix->resize(size, size);
      zero.resize(matrix->rows(), matrix->columns());
      matrix->randomize();
      
      
      //// TEST: Printing description
      
      std::cout << matrix->description() << ":\n" << std::endl;
      
      
      //// TEST: if all elements were copied over
      
      dummyDense = *matrix;
      test(dummyDense.isEqualTo(*matrix), "ELEMENTS EQUAL TO DENSE");
      
      
      //// TEST: if derived class does not equal dense (memorySize comparison)
      
      if (i != 0)
      {
        dummyDense = *matrix;
        test(!(dummyDense == *matrix), "DERIVED NOT EQUAL TO DENSE");
      }
      
      
      //// TEST: Addition
      
      dummyDense = *matrix + *matrix;
      test(dummyDense.isEqualTo(2 * *matrix), "ADDITION");
      
      
      //// TEST: Subtraction
      
      dummyDense = *matrix - *matrix;
      test(dummyDense.isEqualTo(zero), "SUBTRACTION");
      
      
      //// TEST: Scalar multiplication
      
      dummyDense = 0 * *matrix;
      test(dummyDense.isEqualTo(zero), "SCALAR MULTIPLICATION");
      
      
      //// TEST: Matrix multiplication
      
      dummyBool = true;
      dummyDense = *matrix * *matrix;
      
      for (size_t row = 0; row < dummyDense.rows(); row++)
      {
        for (size_t col = 0; col < dummyDense.columns(); col++)
        {
          double sum = 0;
          
          for (size_t j = 0; j < dummyDense.columns(); j++)
          {
            sum += matrix->operator()(row, j) * matrix->operator()(j, col);
          }
          
          if (sum != dummyDense(row, col))
          {
            dummyBool = false;
          }
        }
      }
      
      test(dummyBool, "MATRIX MULTIPLICATION");
      
      
      //// TEST: Clear
      
      matrix->randomize();
      matrix->clear();
      
      test(matrix->isEqualTo(zero), "CLEAR");
      
      std::cout << std::endl << std::endl;
    }
  }
  
  //// The following tests are tested for individual cases
  
  
  //// TEST: Lower matrix multiplication
  
  mLower.randomize();
  dummyLower.randomize();
  
  dummyLower = mLower * mLower;
  dummyDense = dummyLower;
  
  test(dummyLower.isEqualTo(dummyDense), "LOWER MATRIX MULTIPLICATION");
  
  
  //// TEST: Upper matrix multiplication
  
  mUpper.randomize();
  dummyUpper.randomize();
  
  dummyDense = mUpper * mUpper;
  dummyUpper = dummyDense;
  
  test(dummyUpper.isEqualTo(dummyDense), "UPPER MATRIX MULTIPLICATION");
  
  
  //// TEST: Diagonal matrix multiplication
  
  mDiagonal.randomize();
  dummyDiagonal.randomize();
  
  dummyDense = mDiagonal * mDiagonal;
  dummyDiagonal = dummyDense;
  
  test(dummyDiagonal.isEqualTo(dummyDense), "DIAGONAL MATRIX MULTIPLICATION");
  
  
  //// TEST: Upper and lower matrix addition
  
  dummyBool = true;
  mUpper.randomize();
  mLower.randomize();
  
  dummyDense = dummyUpper + dummyLower;
  
  for (size_t row = 0; row < dummyDense.rows() && dummyBool == true; row++)
  {
    for (size_t col = 0; col < dummyDense.columns() && dummyBool == true; col++)
    {
      if (row > col)
      {
        if (dummyDense(row, col) != dummyLower(row, col))
        {
          dummyBool = false;
        }
      }
      else if (row < col)
      {
        if (dummyDense(row, col) != dummyUpper(row, col))
        {
          dummyBool = false;
        }
      }
      else
      {
        if (dummyDense(row, col) != dummyUpper(row, col) + dummyLower(row, col))
        {
          dummyBool = false;
        }
      }
    }
  }
  
  test(dummyBool, "UPPER + LOWER MATRIX ADDITION");
  
  
  //// TEST: Dense matrix memory
  
  dummyBool = true;
  
  for (size_t i = 25; i > 1 && dummyBool == true; i--)
  {
    dummyDense.resize(i, i);
    
    if (dummyDense.memorySize() != (i * i))
    {
      dummyBool = false;
    }
  }
  
  test(dummyBool, "DENSE MATRIX MEMORY");
  
  //// TEST: Triangular matrix memory
  
  dummyBool = true;
  
  for (size_t i = 25; i > 1 && dummyBool == true; i--)
  {
    dummyUpper.resize(i, i);
    
    if (dummyUpper.memorySize() != (((1 + i) * i) / 2))
    {
      dummyBool = false;
    }
  }
  
  test(dummyBool, "TRIANGULAR MATRIX MEMORY");
  
  //// TEST: Banded matrix memory
  
  dummyBool = true;
  
  for (size_t i = 25; i > 1 && dummyBool == true; i--)
  {
    dummyBanded.resize(i, i);
    dummySize = 0;
    
    for (int i = -static_cast<int>(dummyBanded.band()); i <= static_cast<int>(dummyBanded.band()); i++)
    {
      dummySize += dummyBanded.size() - abs(i);
    }
    
    if (dummyBanded.memorySize() != dummySize)
    {
      dummyBool = false;
    }
  }
  
  test(dummyBool, "BANDED MATRIX MEMORY");
  
  dummyBool = true;
  
  for (size_t n = 5; n < 50 && dummyBool == true; n++)
  {
    mPoissonPtr = new matrix_poisson<double>(n);
    
    if (mPoissonPtr->memorySize() != 2)
    {
      dummyBool = false;
    }
    
    if (mPoissonPtr->rows() != (n-1)*(n-1) ||
        mPoissonPtr->columns() != (n-1)*(n-1) )
    {
      dummyBool = false;
    }
    
    for (size_t row = 0; row < mPoissonPtr->rows(); row++)
    {
      for (size_t col = 0; col < mPoissonPtr->columns(); col++)
      {
        if (row == col)
        {
          if (mPoissonPtr->operator()(row, col) != 1)
          {
            dummyBool = false;
          }
        }
        else if (row == col+1 && (row % (mPoissonPtr->slices() - 1) != 0))
        {
          
          for (size_t i = 1; i < row; i++)
          {
            if (i % (mPoissonPtr->slices() - 1) != 0)
            {
              if (mPoissonPtr->operator()(row, col) != -0.25)
              {
                dummyBool = false;
              }
            }
          }
        }
        else if (row == col+(mPoissonPtr->slices() - 1))
        {
          for (size_t i = (mPoissonPtr->slices() - 1); i < row; i++)
          {
            if (mPoissonPtr->operator()(row, col) != -0.25)
            {
              dummyBool = false;
            }
          }
        }
      }
    }
    
    delete mPoissonPtr;
  }
  
  test(dummyBool, "POISSON MATRIX CONSTRUCTION");
}


void test(bool aExpression, std::string aName)
{
  std::cout << std::setw(30) << aName << ": ";
  
  std::cout << (aExpression ? "PASS" : "**** FAIL") << std::endl;
}

