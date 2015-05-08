//
//  Filename:     cholesky.hpp
//  Programmer:   Nicholas Trampe
//  Class:        CS 5201 - Clayton Price
//  Assignment:   Final - A Parameterized Matrix Class and Gaussian Elimination
//
//  Description:  Here are the definitions of the functions declared in the cholesky class
//
//


template <class T>
bool cholesky<T>::operator()(vector<T>& aX, const matrix_base<T>& aA, const vector<T>& aB)
{
  matrix_symmetrical<T> aL(aA);
  matrix_symmetrical<T> aAug;
  vector<T> y;
  double sum = 0;
  int size = static_cast<int>(aL.size());
  
  y = aB;
  
  for (int k = 0; k < size; k++)
  {
    for (int i = 0; i <= k-1; i++)
    {
      sum = 0;
      for (int j = 0; j <= i-1; j++)
      {
        sum = sum + aL(i,j) * aL(k,j);
      }
      aL(k,i) = (aL(k,i) - sum)/aL(i,i);
    }
    
    sum = 0;
    
    for (int j = 0; j <= k-1; j++)
    {
      sum = sum + aL(k,j) * aL(k,j);
    }
    
    aL(k,k) = sqrt(aL(k,k) - sum);
  }
  
  aAug = aL;
  
  for (size_t col = 0; col < aAug.size(); col++)
  {
    if (aAug(col, col) != 0 && aAug(col, col) != 1)
    {
      y[col] = y[col] / aAug(col, col);
      
      aAug(col, col) = 1;
    }
    
    for (size_t row = col + 1; row < aAug.size(); row++)
    {
      y[row] = -aAug(row, col) * y[col] + y[row];
      
      aAug(row, col) = 0;
    }
  }
  
  double csum = 0;
  
  // create solution array
  aX.clear();
  aX.reserve(aL.columns(), true);
  
  // backward substitution
  // solve for the unknown variables from the bottom up
  for (int row = (int)aL.columns() - 1; row >= 0; row--)
  {
    csum = 0;
    
    // compute sum of the coefficients times the old
    // x values
    for (size_t col = row; col < aL.columns(); col++)
    {
      csum += aL(row, col) * aX[col];
    }
    
    // compute the new x by subtracting the sum from the
    // right-hand side of the equation and dividing by the
    // coefficient corresponding to the current x
    //
    // xi = (an - (a1 + a2 + ... + an-1 (not including ai))) / ai
    //
    // where i = row
    aX[row] = (y[row] - csum) / aL(row, row);
    
    // fix rounding errors (e.g. 0 != 0.000000000001)
    // for display purposes
    if (equivalent(aX[row], 0))
    {
      aX[row] = 0;
    }
    
    if (aX[row] != aX[row])
    {
      // aX[i] is not a number
      // no solution
      return false;
    }
  }
  
  // aX now contains the solution
  return true;
}

