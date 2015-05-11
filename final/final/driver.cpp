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
  runMenu();
  
  return 0;
}


void runMenu()
{
  pde_final<double> pde(20, 0,M_PI);
  kMenuChoice choice = kMenuChoiceQuit;
  std::string input;
  float animFactor = 0.0;
  bool drawLines = false;
  bool outputToFile = false;
  bool shouldIterate = false;
  bool prettyPrint = false;
  point2d<double> bounds = pde.bounds();
  std::ofstream file;
  
  do
  {
    std::cout << std::string(2, '\n') << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    std::cout << "Object-Oriented Numerical Modeling Final" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    std::cout << "1. Change Mesh Density" << std::endl;
    std::cout << "2. Change Bounds" << std::endl;
    std::cout << "3. Compare Techniques" << std::endl;
    std::cout << "4. Solve" << std::endl;
    std::cout << "5. Output Matlab" << std::endl;
    std::cout << "6. Run Tests" << std::endl;
    std::cout << "0. Quit" << std::endl;
    std::cout << std::string(25, '-') << std::endl;
    
    std::cout << "\nMesh Density: " << pde.density() << "\nBounds: " << pde.bounds() << (pde.solved() ? "\n*** PDE Solution Stored ***" : "") << std::endl;
    
    std::cout << "\nChoice: ";
    
    std::cin >> input;
    choice = static_cast<kMenuChoice>(atoi(input.c_str()));
    
    
    switch (choice)
    {
      case kMenuChoiceChangeMeshDensity:
        std::cout << "Enter Mesh Density:" << std::endl;
        std::cin >> input;
        
        pde.setDensity(static_cast<size_t>(atoi(input.c_str())));
        
        break;
        
      case kMenuChoiceChangeBounds:
        std::cout << "Enter Bounds: 'X Y'" << std::endl;
        std::cin >> bounds;
        
        pde.setBounds(bounds);
        
        break;
        
      case kMenuChoiceCompare:
        
        std::cout << "Should Iterate? (1 or 0): ";
        std::cin >> shouldIterate;
        
        std::cout << "Print Pretty? (1 or 0): ";
        std::cin >> prettyPrint;
        
        printMessage("Comparing Techniques...");
        
        runSolvers(pde, shouldIterate, prettyPrint);
        
        break;
        
      case kMenuChoiceSolve:
        
        promptForSolve(pde);
        
        if (pde.solved())
        {
          std::cout << pde.pointsOutput() << std::endl;
        }
        
        break;
        
      case kMenuChoiceMatlab:
        
        promptForSolve(pde);
        
        if (pde.solved())
        {
          std::cout << "Enter Animation Factor: ";
          std::cin >> animFactor;
          
          std::cout << "Should Draw Lines? (1 or 0): ";
          std::cin >> drawLines;
          
          std::cout << "Output to file? (1 or 0): ";
          std::cin >> outputToFile;
          
          if (outputToFile)
          {
            file.open("matlab_output.txt");
            file.clear();
            file << pde.matlabOutput(animFactor, drawLines) << std::endl;
            std::cout << "\n*** Matlab output stored in 'matlab_output.txt' ***" << std::endl;
          }
          else
          {
            std::cout << pde.matlabOutput(animFactor, drawLines) << std::endl;
          }
        }
        
        break;
        
      case kMenuChoiceTest:
        
        printMessage("Testing...");
        
        testMatrices();
        
        std::cin.ignore();
        
        break;
        
      default:
        break;
    }
    
  } while (choice != kMenuChoiceQuit);
}


void promptForSolve(pde_base<double>& aPDE)
{
  int method = 0;
  
  std::cout << "\nChoose Technique:" << std::endl;
  std::cout << std::string(40, '-') << std::endl;
  std::cout << "1. Cholesky" << std::endl;
  std::cout << "2. Gaussian Elimination" << std::endl;
  std::cout << "3. Gauss-Seidel Iterative Method" << std::endl;
  
  if (aPDE.solved())
  {
    std::cout << "4. Use Current Solution" << std::endl;
  }
  
  std::cout << "0. Back" << std::endl;
  std::cout << std::string(40, '-') << std::endl;
  
  std::cin >> method;
  
  if (method != 0)
  {
    printMessage("Solving...");
  }
  
  switch (method)
  {
    case 1:
      solvePDE(aPDE, cholesky<double>());
      break;
      
    case 2:
      solvePDE(aPDE, gauss_elim<double>());
      break;
      
    case 3:
      solvePDE(aPDE, gauss_seidel<double>());
      break;
      
    default:
      break;
  }
}


void printMessage(const std::string& aMessage)
{
  std::string line = std::string(aMessage.length(), '-');
  
  std::cout << line << std::endl;
  std::cout << aMessage << std::endl;
  std::cout << line << "\n\n" << std::endl;
}


void createSystem(const pde_base<double>& aPDE, vector<double>& aB, vector<point2d<double>>& aXMapping)
{
  const size_t SIZE = (aPDE.density()-1)*(aPDE.density()-1);
  double lowerBound = aPDE.lowerBound();
  double upperBound = aPDE.upperBound();
  double inc = fabs((upperBound - lowerBound)) / aPDE.density();
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


template <class T_method>
void solvePDE(pde_base<double>& aPDE, T_method aMethod)
{
  matrix_poisson<double> m(aPDE.density());
  vector<point2d<double>> xMapping;
  vector<double> b;
  vector<double> x;
  point2d<double> p;
  
  aPDE.clearPoints();
  
  createSystem(aPDE, b, xMapping);
  
  solveMatrix(x, m, b, aMethod);
  
  for (size_t i = 0; i < x.size(); i++)
  {
    aPDE.addPoint(xMapping[i], x[i]);
  }
}


void runSolvers(pde_base<double>& aPDE, const bool aShouldIterate, const bool aPrettyPrint)
{
  // NOTE:  a lot of redundant output code
  //        a possible improvement could be a table class to easily output tables
  
  const size_t DENSITY = aPDE.density();
  const int TABLE_1_SIZE = 5;
  const int TABLE_2_SIZE = 3;
  const int TABLE_1_WIDTHS[TABLE_1_SIZE] = {5,10,22,15,16};
  const int TABLE_2_WIDTHS[TABLE_2_SIZE] = {16,10,16};
  const char* TABLE_1_HEADERS[TABLE_1_SIZE] =
  {"N", "Cholesky", "Gaussian Elimination", "Gauss-Seidel", "Solution Error"};
  const char* TABLE_2_HEADERS[TABLE_2_SIZE] =
  {"Error Tolerance", "Time", "Solution Error"};
  const int SEIDEL_ITERATIONS = 12;
  
  matrix_poisson<double>* m;
  vector<point2d<double>> xMapping;
  vector<double> b;
  vector<double> x;
  vector<double> times;
  runtime timer;
  double seidel_tol = 0;
  size_t start = (aShouldIterate ? 2 : DENSITY);
  std::stringstream ss;
  std::string ssLine;
  
  std::cout << "Variable Mesh Density (N)" << std::endl;
  
  if (aPrettyPrint)
  {
    ss << "|";
    
    for (int i = 0; i < TABLE_1_SIZE; i++)
    {
      ss << std::setw(TABLE_1_WIDTHS[i]) << TABLE_1_HEADERS[i] << " |";
    }
    
    ssLine = std::string(ss.str().length(), '-');
    
    std::cout << ssLine << " \n" << ss.str() << "\n" << ssLine << std::endl;
  }
  else
  {
    for (int i = 0; i < TABLE_1_SIZE; i++)
    {
      std::cout << TABLE_1_HEADERS[i] << "\t";
    }
    
    std::cout << std::endl << std::endl;
  }
  
  for (size_t n = start; n <= DENSITY; n++)
  {
    m = new matrix_poisson<double>(n);
    aPDE.setDensity(n);
    createSystem(aPDE, b, xMapping);
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
    
    if (aPrettyPrint)
    {
      std::cout << "|" << std::setw(TABLE_1_WIDTHS[0]) << n << " |"
      << std::setw(TABLE_1_WIDTHS[1]) << times[0] << " |"
      << std::setw(TABLE_1_WIDTHS[2]) << times[1] << " |"
      << std::setw(TABLE_1_WIDTHS[3]) << times[2] << " |"
      << std::setw(TABLE_1_WIDTHS[4]) << checkError(x, xMapping) << " |"
      << std::endl;
    }
    else
    {
      std::cout << n
      << "\t" << times[0]
      << "\t" << times[1]
      << "\t" << times[2]
      << "\t" << checkError(x, xMapping)
      << std::endl;
    }
  }
  
  std::cout << ssLine << std::endl;
  
  
  ss.clear();
  ss.str(std::string());
  
  std::cout << "\n\nVariable Error Tolerance (N = " << DENSITY << "):" << std::endl;
  
  if (aPrettyPrint)
  {
    ss << "|";
    
    for (int i = 0; i < TABLE_2_SIZE; i++)
    {
      ss << std::setw(TABLE_2_WIDTHS[i]) << TABLE_2_HEADERS[i] << " |";
    }
    
    ssLine = std::string(ss.str().length(), '-');
    
    std::cout << ssLine << " \n" << ss.str() << "\n" << ssLine << std::endl;
  }
  else
  {
    for (int i = 0; i < TABLE_2_SIZE; i++)
    {
      std::cout << TABLE_2_HEADERS[i] << "\t";
    }
    
    std::cout << std::endl << std::endl;
  }
  
  for (double i = 1; i <= SEIDEL_ITERATIONS; i++)
  {
    m = new matrix_poisson<double>(DENSITY);
    aPDE.setDensity(DENSITY);
    seidel_tol = pow(10, -i);
    createSystem(aPDE, b, xMapping);
    
    timer.begin();
    solveMatrix(x, *m, b, gauss_seidel<double>(seidel_tol));
    timer.end();
    
    delete m;
    
    if (aPrettyPrint)
    {
      std::cout << "|" << std::setw(TABLE_2_WIDTHS[0]) << seidel_tol << " |"
      << std::setw(TABLE_2_WIDTHS[1]) << timer.elapsed() << " |"
      << std::setw(TABLE_2_WIDTHS[2]) << checkError(x, xMapping) << " |"
      << std::endl;
    }
    else
    {
      std::cout << seidel_tol
      << "\t" << timer.elapsed()
      << "\t" << checkError(x, xMapping)
      << std::endl;
    }
  }
  
  std::cout << ssLine << std::endl;
  
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

