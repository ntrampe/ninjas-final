//
//  Filename:
//  Programmer:   Nicholas Trampe, James Kellerman
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - Solving Poisson's Equation
//
//  Description:  This is the driver implementation for Assignment 6
//                This program reads in a matrix, vector and optionally a matrix type
//                then solves the system
//

#include "driver.h"

int main(int argc, const char * argv[])
{ 
  matrix_poisson<double> m(4);
  vector<point<double>> x;
  point<double> p;
  size_t s = m.slices();
  
  
  for (size_t i = 1; i < s; i++)
  {
    for (size_t j = 1; j < s; j++)
    {
      p.set(j, i);
      x.push_back(p);
    }
  }
  
  //m.order(1);
  
  std::cout << m << std::endl;
  m.printMemory();
  std::cout << x << std::endl;
  
  return 0;
  
  // file containing matrix data
  std::string fileName;
  std::string typeStr;
  kMatrixType type;
  
  seedRandom();
  
  if (argc < 2)
  {
    std::cout << "Please specify a file:" << std::endl;
    std::cin >> fileName;
  }
  else
  {
    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0)
    {
      std::cout << "\nusage: ./driver <file> <matrix type>\n" << std::endl;
      std::cout << "<file> = file containing matrix and vector\n" << std::endl;
      std::cout << "<matrix type> = number representing matrix type as follows:\n" << std::endl;
      displayMatrixTypes();
      std::cout << std::endl;
      
      return 0;
    }
    else
    {
      fileName = argv[1];
    }
  }
  
  // the user has the option to enter a matrix type
  // if the user doesn't enter anything the matrix will be
  // treated as dense
  
  if (argc < 3)
  {
    type = kMatrixTypeDense;
  }
  else
  {
    type = static_cast<kMatrixType>(atoi(argv[2]));
  }
  
  std::cout << "\nFinding solution for input...\n" << std::endl;
  
  solveFile(fileName.c_str(), type);
  
  std::cout << "\nTesting...\n" << std::endl;
  
  testMatrices();
  
  return 0;
}


void solveFile(const char * aFile, kMatrixType aType)
{
  // coefficient matrix
  matrix_base<double>* a;
  
  // presentation matrix
  // for displaying to user
  // includes b vector
  matrix<double> p;
  
  // b vector
  vector<double> b;
  
  // solution vector
  vector<double> x;
  
  switch (aType)
  {
    case kMatrixTypeDense:
      a = new matrix<double>;
      break;
      
    case kMatrixTypeTriangularUpper:
      a = new matrix_triangular_upper<double>;
      break;
      
    case kMatrixTypeTriangularLower:
      a = new matrix_triangular_lower<double>;
      break;
      
    case kMatrixTypeDiagonal:
      a = new matrix_diagonal<double>;
      break;
      
    case kMatrixTypeTridiagonal:
      a = new matrix_tridiagonal<double>;
      break;
      
    case kMatrixTypeSymmetrical:
      a = new matrix_symmetrical<double>;
      break;
      
    default:
      break;
  }
  
  // open file
  if (openFile(a, b, aFile))
  {
    // copy a over to b
    p = *a;
    
    // create augmented matrix
    p.insertVectorAtColumn(b, p.columns());
    
    // display augmented matrix to user
    std::cout << "Augmented " << a->name() << " Input:" << std::endl << std::endl;
    std::cout << p << std::endl;
    std::cout << std::endl;
    
    // perform Gaussian Elimination
    if (a->solveMatrix(b, x))
    {
      std::cout << "Solution: x = " << x << std::endl;
    }
    else
    {
      std::cout << "No Solution :(" << std::endl;
    }
  }
  else
  {
    std::cout << "The file '" << aFile << "' could not be found!" << std::endl;
  }
  
  delete a;
}


bool openFile(matrix_base<double>* aMatrix, vector<double>& aVector, const char * aFile)
{
  std::ifstream file(aFile);
  size_t size = 0;
  
  aVector.clear();
  
  if (file.is_open())
  {
    file >> size;
    
    aMatrix->resize(size, size);
    aVector.reserve(size, true);
    
    file >> *aMatrix;
    file >> aVector;
    
    file.close();
    
    return true;
  }
  
  return false;
}


void displayMatrixTypes()
{
  std::cout << "0 = Dense" << std::endl;
  std::cout << "1 = Upper Triangular" << std::endl;
  std::cout << "2 = Lower Triangular" << std::endl;
  std::cout << "3 = Diagonal" << std::endl;
  std::cout << "4 = Tridiagonal" << std::endl;
  std::cout << "5 = Symmetrical" << std::endl;
}


void testMatrices()
{
  //TODO: use gtest or something because this is horrible
  
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
  
  // insert matrices into vector
  mats.push_back(&mDense);
  mats.push_back(&mLower);
  mats.push_back(&mUpper);
  mats.push_back(&mBanded);
  mats.push_back(&mDiagonal);
  mats.push_back(&mTridiagonal);
  mats.push_back(&mSymmetrical);
  
  //// TEST: tests
  
  ////  The following tests are tested for every matrix in the array
  ////  for both even and odd dimensions
  
  ////  this also tests polymorphism
  
  for (size_t size = START_SIZE; size > START_SIZE - 2; size--)
  {
    zero.resize(size, size);
    
    for (size_t i = 0; i < mats.size(); i++)
    {
      matrix_base<double>* matrix = mats[i];
      
      matrix->resize(size, size);
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
}


void test(bool aExpression, std::string aName)
{
  std::cout << std::setw(30) << aName << ": ";
  
  std::cout << (aExpression ? "PASS" : "**** FAIL") << std::endl;
}








